module RtmFingerprint

!-----------------------------------------------------------------------
! !DESCRIPTION:
! TEMPORARY DIAGNOSTIC MODULE -- NOT FOR MERGE.
!
! Purpose: localize an exact-restart (ERS) failure in MOSART by emitting a
! bitwise fingerprint of selected state arrays at chosen points in the run.
!
! Method: reinterpret each real(r8) as a 64-bit integer and XOR-reduce over
! all local cells, then XOR-reduce across MPI ranks. The result is a single
! integer that is *bit-exactly* reproducible for a given set of field values
! and is completely independent of:
!   - domain decomposition (XOR is commutative and associative)
!   - loop / thread / rank ordering
!   - print formatting or precision loss
!
! Two runs whose fingerprints agree have bit-identical fields. Two runs whose
! fingerprints differ have at least one differing bit somewhere in the field.
! Because the model has been shown to be run-to-run deterministic, the FIRST
! call at which the continuous ("base") leg and the restarted ("rest") leg
! disagree identifies the first divergent quantity unambiguously.
!
! A companion count of cells with a nonzero value is also emitted, purely as a
! sanity check that the array is actually populated (a fingerprint of 0 could
! otherwise mean "all zeros" or "never touched").
!
! Usage: call rtm_fp(label, array) or rtm_fp2(label, array2d, itracer).
! Output lines are tagged 'RTMFP' for easy grep/diff:
!    grep RTMFP rof.log.base > base.fp
!    grep RTMFP rof.log.rest > rest.fp
!    diff base.fp rest.fp | head
!-----------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8, i8 => shr_kind_i8
  use RtmVar        , only : iulog
  use RtmSpmd       , only : masterproc, mpicom_rof, npes, iam
  use RunoffMod     , only : rtmCTL
  use RtmTimeManager, only : get_nstep, get_curr_date
  use mpi

  implicit none
  private

  public :: rtm_fp      ! fingerprint a 1-d real(r8) array
  public :: rtm_fp2     ! fingerprint one tracer slice of a 2-d real(r8) array
  public :: rtm_fp_tag  ! emit a bare marker line (no field)
  public :: rtm_fp_dam  ! fingerprint a dam-dimensioned array (bypasses cell guard)
  public :: rtm_fp_rep  ! fingerprint an array already replicated on every rank
  public :: rtm_fp_part ! fingerprint a rank-LOCAL PARTIAL array over a global index
  public :: rtm_fp_int  ! fingerprint a 1-d INTEGER array
  public :: rtm_fp_int2 ! fingerprint a 2-d INTEGER array, e.g. WRMUnit%myDam
  public :: rtm_fp_dump ! dump per-cell values to a side file, to NAME the cell
  public :: rtm_fp_gather ! gather a rank-local PARTIAL array from ALL ranks raw

  logical, public :: rtm_fp_on = .true.   ! master switch

  ! --- per-cell dump control (rtm_fp_dump) ---
  ! The dump is expensive and verbose, so it is restricted to a narrow nstep
  ! window bracketing the first known divergence. Widen only if needed.
  !
  ! Sizing: ~60k nonzero cells x 28 bytes x 4 labels x 3 Euler calls is
  ! ~20 MB per nstep per leg, so a 3-nstep window costs ~120 MB for both
  ! legs. Disk quota has already broken one build in this case directory --
  ! check quota before widening this window.
  !
  ! The window must bracket the FIRST divergence, which differs by build:
  !   debug     build: first divergence at nstep 54 -> window 53..55
  !   optimized build: first divergence at nstep 51 -> window 50..52
  integer, public :: rtm_fp_dump_lo = 48  ! first nstep to dump (inclusive)
  integer, public :: rtm_fp_dump_hi = 50  ! last  nstep to dump (inclusive)

  ! --- raw cross-rank gather control (rtm_fp_gather) ---
  ! Every fingerprint is a checksum, and every checksum has a blind spot. The
  ! one that matters here is a cross-rank index swap: rank A moves a value from
  ! dam index i to j while rank B moves an equal value from j to i. The raw
  ! XOR, the index-mixed hash and the rank-mixed hash are ALL unchanged by
  ! that, which is exactly the "attribution changes while the total is
  ! conserved" scenario under investigation. Only the raw bits from every rank
  ! settle it.
  !
  ! Sizing: npes x NDam x 8 bytes = 128 x 4246 x 8 ~ 4.3 MB per call, so this
  ! is restricted to ONE (nstep, iter) event, not a window.
  integer, public :: rtm_fp_gather_nstep = 49  ! the only nstep to gather at
  integer, public :: rtm_fp_gather_iter  = 1   ! the only iteration to gather at

  ! The two ERS legs share one RUNDIR, so dump filenames must not collide.
  ! The legs are distinguished by the first nstep each one ever sees (base
  ! starts at 0, rest starts at the restart nstep), captured on first use.
  integer, save :: fp_nstep0 = -1

contains

!-----------------------------------------------------------------------

  subroutine rtm_fp_hdr(sub, cycle)
    ! Write the common "where am I in the run" prefix.
    character(len=*), intent(in) :: sub
    integer, intent(in) :: cycle
    integer :: yr, mon, day, tod

    call get_curr_date(yr, mon, day, tod)
    write(iulog,'(a,i8,1x,i4.4,a,i2.2,a,i2.2,1x,i6,1x,a,1x,i4)') &
         'RTMFP nstep=', get_nstep(), yr, '-', mon, '-', day, tod, &
         trim(sub), cycle
  end subroutine rtm_fp_hdr

!-----------------------------------------------------------------------

  subroutine rtm_fp_reduce(arr, gidx, xorloc, nnzloc)
    ! XOR-reduce the bit patterns of arr, and count nonzeros.
    !
    ! Two reductions are combined:
    !   (a) plain XOR of the bit patterns -- catches any changed value;
    !   (b) XOR of (bits mixed with the cell's GLOBAL index) -- additionally
    !       catches the case where the same multiset of values ends up on
    !       different cells, which plain XOR is blind to (XOR is permutation
    !       insensitive, and a pair of equal values cancels to zero).
    ! Using the *global* index keeps (b) decomposition independent.
    real(r8), intent(in)  :: arr(:)
    integer,  intent(in)  :: gidx(:)
    integer(i8), intent(out) :: xorloc(2)
    integer(i8), intent(out) :: nnzloc(1)

    integer     :: i
    integer(i8) :: bits, mixed

    xorloc = 0_i8
    nnzloc = 0_i8
    do i = 1, size(arr)
       ! TRANSFER gives the raw IEEE-754 bit pattern; no arithmetic, so no
       ! rounding and no compiler reassociation can perturb it.
       bits      = transfer(arr(i), bits)
       xorloc(1) = ieor(xorloc(1), bits)
       ! Multiply the value bits by an odd function of the global index. The
       ! product is integer arithmetic (exact, wrapping) so it is still fully
       ! reproducible, but it now depends on WHERE the value sits.
       mixed     = bits * (2_i8 * int(gidx(i), i8) + 1_i8)
       xorloc(2) = ieor(xorloc(2), mixed)
       if (arr(i) /= 0._r8) nnzloc(1) = nnzloc(1) + 1_i8
    end do
  end subroutine rtm_fp_reduce

!-----------------------------------------------------------------------

  subroutine rtm_fp(label, arr, cycle)
    ! Emit a global bitwise fingerprint of a 1-d array.
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: arr(:)
    integer, optional, intent(in) :: cycle

    ! All reduction buffers are rank-1 so the MPI choice-buffer arguments are
    ! unambiguous for every compiler's mpi module interface.
    integer(i8) :: xorloc(2), nnzloc(1), xorglb(2), nnzglb(1)
    integer     :: ier, mycycle

    if (.not. rtm_fp_on) return

    mycycle = 0
    if (present(cycle)) mycycle = cycle

    if (size(arr) /= size(rtmCTL%gindex)) then
       ! Guard: the index-mixed hash requires a cell-dimensioned array.
       if (masterproc) write(iulog,'(a,a)') &
            'RTMFP WARNING: size mismatch, skipping fld=', trim(label)
       return
    endif

    call rtm_fp_reduce(arr, rtmCTL%gindex, xorloc, nnzloc)

    ! MPI_BXOR over ranks: also order-independent, so the answer does not
    ! depend on rank count or which cells live where.
    call mpi_allreduce(xorloc, xorglb, 2, MPI_INTEGER8, MPI_BXOR, mpicom_rof, ier)
    call mpi_allreduce(nnzloc, nnzglb, 1, MPI_INTEGER8, MPI_SUM,  mpicom_rof, ier)

    if (masterproc) then
       call rtm_fp_hdr(label, mycycle)
       write(iulog,'(a,a24,a,z16.16,a,z16.16,a,i10)') &
            'RTMFP   fld=', adjustl(label), ' xor=0x', xorglb(1), &
            ' hash=0x', xorglb(2), ' nnz=', nnzglb(1)
    endif
  end subroutine rtm_fp

!-----------------------------------------------------------------------

  subroutine rtm_fp2(label, arr, nt, cycle)
    ! Emit a global bitwise fingerprint of one tracer slice of a 2-d array.
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: arr(:,:)
    integer, intent(in) :: nt
    integer, optional, intent(in) :: cycle

    if (.not. rtm_fp_on) return
    call rtm_fp(label, arr(:,nt), cycle)
  end subroutine rtm_fp2

!-----------------------------------------------------------------------

  subroutine rtm_fp_dam(label, arr, gidx, cycle)
    ! Emit a global bitwise fingerprint of a DAM-dimensioned array.
    !
    ! rtm_fp() deliberately refuses arrays that are not cell-dimensioned, so
    ! the entire WRM dam state (LocalNumDam) has no coverage. This variant
    ! takes the global dam ID array explicitly (WRMUnit%damID) for the
    ! index-mixed hash, which keeps the result decomposition independent
    ! without RtmFingerprint having to depend on WRM_type_mod.
    !
    ! Note LocalNumDam may be 0 on some ranks; the local reduction then
    ! contributes identity (0) to the XOR, which is correct.
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: arr(:)
    integer,  intent(in) :: gidx(:)
    integer, optional, intent(in) :: cycle

    integer(i8) :: xorloc(2), nnzloc(1), xorglb(2), nnzglb(1)
    integer(i8) :: nloc(1), nglb(1)
    integer     :: ier, mycycle

    if (.not. rtm_fp_on) return

    mycycle = 0
    if (present(cycle)) mycycle = cycle

    if (size(arr) /= size(gidx)) then
       if (masterproc) write(iulog,'(a,a)') &
            'RTMFP WARNING: dam size mismatch, skipping fld=', trim(label)
       return
    endif

    call rtm_fp_reduce(arr, gidx, xorloc, nnzloc)
    nloc(1) = int(size(arr), i8)

    call mpi_allreduce(xorloc, xorglb, 2, MPI_INTEGER8, MPI_BXOR, mpicom_rof, ier)
    call mpi_allreduce(nnzloc, nnzglb, 1, MPI_INTEGER8, MPI_SUM,  mpicom_rof, ier)
    call mpi_allreduce(nloc,   nglb,   1, MPI_INTEGER8, MPI_SUM,  mpicom_rof, ier)

    if (masterproc) then
       call rtm_fp_hdr(label, mycycle)
       write(iulog,'(a,a24,a,z16.16,a,z16.16,a,i10,a,i10)') &
            'RTMFP   fld=', adjustl(label), ' xor=0x', xorglb(1), &
            ' hash=0x', xorglb(2), ' nnz=', nnzglb(1), ' ndam=', nglb(1)
    endif
  end subroutine rtm_fp_dam

!-----------------------------------------------------------------------

  subroutine rtm_fp_rep(label, arr, cycle)
    ! Emit a bitwise fingerprint of an array that is ALREADY REPLICATED
    ! identically on every rank -- e.g. the output of shr_mpi_sum(...,all=.true.)
    ! or of mct_aVect_bcast. Such arrays are global in the dam index, not
    ! decomposed, so there is nothing to reduce across ranks: masterproc simply
    ! hashes its own copy.
    !
    ! This closes the last blind spot in the diagnostic. rtm_fp/rtm_fp_dam are
    ! reductions over LOCAL cells/dams, so a quantity that is global-indexed and
    ! replicated (dam_uptake_sum, aVect_wdG%rAttr) is invisible to them: two runs
    ! can deliver an identical total to identical gridcells while ATTRIBUTING it
    ! to different dams, and every cell-dimensioned fingerprint stays equal.
    ! That is exactly the failure mode under investigation.
    !
    ! The index used for the mixed hash is the array's own position, which is
    ! the GLOBAL dam index for these arrays -- already decomposition
    ! independent, precisely because the array is replicated.
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: arr(:)
    integer, optional, intent(in) :: cycle

    integer     :: i, mycycle
    integer(i8) :: bits, mixed, xr(2), nnz

    if (.not. rtm_fp_on) return

    mycycle = 0
    if (present(cycle)) mycycle = cycle

    ! No MPI here: every rank holds the same data, so a reduction would only
    ! XOR npes identical copies together (and silently vanish for even npes).
    if (.not. masterproc) return

    xr  = 0_i8
    nnz = 0_i8
    do i = 1, size(arr)
       bits  = transfer(arr(i), bits)
       xr(1) = ieor(xr(1), bits)
       mixed = bits * (2_i8 * int(i, i8) + 1_i8)
       xr(2) = ieor(xr(2), mixed)
       if (arr(i) /= 0._r8) nnz = nnz + 1_i8
    end do

    call rtm_fp_hdr(label, mycycle)
    write(iulog,'(a,a24,a,z16.16,a,z16.16,a,i10,a,i10)') &
         'RTMFP   fld=', adjustl(label), ' xor=0x', xr(1), &
         ' hash=0x', xr(2), ' nnz=', nnz, ' nrep=', int(size(arr), i8)
  end subroutine rtm_fp_rep

!-----------------------------------------------------------------------

  subroutine rtm_fp_part(label, arr, cycle)
    ! Emit a bitwise fingerprint of a rank-LOCAL PARTIAL array laid out over a
    ! GLOBAL index -- the canonical case being dam_uptake(NDam) BEFORE its
    ! shr_mpi_sum. Every rank allocates the full global extent but writes only
    ! the entries its own gridcells contribute to, leaving the rest zero.
    !
    ! This is the probe that rtm_fp_rep must NOT be used for, and was: rtm_fp_rep
    ! hashes masterproc's copy alone, which for a partial array measures one
    ! rank's contribution and silently discards the other 127. That is how
    ! G_it_uptakeloc came out nnz=0 on all 792 records and carried no
    ! information.
    !
    ! Four accumulators are emitted because plain XOR has blind spots that all
    ! bite here:
    !   xor   -- raw XOR of the value bits. Blind to permutation, and two ranks
    !            contributing the SAME value to the SAME dam cancel to zero.
    !   hash  -- bits mixed with the global slot index. Fixes permutation.
    !   rhash -- bits mixed with an odd function of the RANK. Fixes the equal-
    !            contributions-on-different-ranks cancellation, which is not
    !            hypothetical: identical demand at mirrored gridcells is exactly
    !            the kind of symmetry this model produces.
    !   xrhash-- bits mixed with slot AND rank TOGETHER. This closes a blind spot
    !            the first three share: if rank A moves a value from dam index
    !            i to j while rank B moves an equal value from j to i, then the
    !            per-index totals, the per-rank totals, and the raw XOR are all
    !            unchanged, and only a joint (slot,rank) mix detects it. That
    !            cross-rank index swap is exactly the "attribution changes while
    !            the total is conserved" scenario under investigation, so a
    !            clean uptakeloc without this accumulator does NOT exclude a
    !            permuted input.
    ! nnz is SUMMED across ranks, so it counts total nonzero contributions
    ! rather than distinct dams -- a direct check that the array is populated.
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: arr(:)
    integer, optional, intent(in) :: cycle

    integer     :: i, ier, mycycle
    integer(i8) :: bits, xorloc(4), nnzloc(1), xorglb(4), nnzglb(1)

    if (.not. rtm_fp_on) return

    mycycle = 0
    if (present(cycle)) mycycle = cycle

    xorloc = 0_i8
    nnzloc = 0_i8
    do i = 1, size(arr)
       bits      = transfer(arr(i), bits)
       xorloc(1) = ieor(xorloc(1), bits)
       xorloc(2) = ieor(xorloc(2), bits * (2_i8 * int(i, i8) + 1_i8))
       xorloc(3) = ieor(xorloc(3), bits * (2_i8 * int(iam, i8) + 1_i8))
       xorloc(4) = ieor(xorloc(4), bits * (2_i8 * int(i,   i8) + 1_i8) &
                                        * (2_i8 * int(iam, i8) + 1_i8))
       if (arr(i) /= 0._r8) nnzloc(1) = nnzloc(1) + 1_i8
    end do

    ! XOR across ranks: order independent, and zero contributes identity, so
    ! ranks that touch no dam drop out cleanly.
    call mpi_allreduce(xorloc, xorglb, 4, MPI_INTEGER8, MPI_BXOR, mpicom_rof, ier)
    call mpi_allreduce(nnzloc, nnzglb, 1, MPI_INTEGER8, MPI_SUM,  mpicom_rof, ier)

    if (masterproc) then
       call rtm_fp_hdr(label, mycycle)
       write(iulog,'(a,a24,a,z16.16,a,z16.16,a,z16.16,a,z16.16,a,i10)') &
            'RTMFP   fld=', adjustl(label), ' xor=0x', xorglb(1), &
            ' hash=0x', xorglb(2), ' rhash=0x', xorglb(3), &
            ' xrhash=0x', xorglb(4), ' nnz=', nnzglb(1)
    endif
  end subroutine rtm_fp_part

!-----------------------------------------------------------------------

  subroutine rtm_fp_int(label, iarr, gidx, cycle)
    ! Emit a bitwise fingerprint of a 1-d INTEGER array.
    !
    ! Every other probe in this module takes real(r8), so the WRM integer
    ! dependency arrays -- WRMUnit%myDamNum, %icell, %damID -- have never been
    ! fingerprinted at all. They are the only inputs to the dam-uptake
    ! attribution that no probe has touched. A permuted or differently-built
    ! dependency list changes which dam receives which gridcell's uptake while
    ! leaving every real-valued field bit-identical, which is precisely the
    ! signature under investigation.
    !
    ! gidx supplies the GLOBAL index for the mixed hash, so the result is
    ! decomposition independent. Pass rtmCTL%gindex for cell-dimensioned arrays
    ! and WRMUnit%damID for dam-dimensioned ones. For arrays whose VALUES are
    ! local indices (WRMUnit%icell), map the values through rtmCTL%gindex at the
    ! call site -- a raw local index is not comparable across decompositions.
    !
    ! No TRANSFER: integers are already exact, so the value is used directly.
    ! Zero values contribute identity, which is what we want for the zero
    ! padding in myDam.
    character(len=*), intent(in) :: label
    integer, intent(in) :: iarr(:)
    integer, intent(in) :: gidx(:)
    integer, optional, intent(in) :: cycle

    integer     :: i, ier, mycycle
    integer(i8) :: v, xorloc(2), nnzloc(1), xorglb(2), nnzglb(1)

    if (.not. rtm_fp_on) return

    mycycle = 0
    if (present(cycle)) mycycle = cycle

    if (size(iarr) /= size(gidx)) then
       if (masterproc) write(iulog,'(a,a)') &
            'RTMFP WARNING: int size mismatch, skipping fld=', trim(label)
       return
    endif

    xorloc = 0_i8
    nnzloc = 0_i8
    do i = 1, size(iarr)
       v         = int(iarr(i), i8)
       xorloc(1) = ieor(xorloc(1), v)
       xorloc(2) = ieor(xorloc(2), v * (2_i8 * int(gidx(i), i8) + 1_i8))
       if (iarr(i) /= 0) nnzloc(1) = nnzloc(1) + 1_i8
    end do

    call mpi_allreduce(xorloc, xorglb, 2, MPI_INTEGER8, MPI_BXOR, mpicom_rof, ier)
    call mpi_allreduce(nnzloc, nnzglb, 1, MPI_INTEGER8, MPI_SUM,  mpicom_rof, ier)

    if (masterproc) then
       call rtm_fp_hdr(label, mycycle)
       write(iulog,'(a,a24,a,z16.16,a,z16.16,a,i10)') &
            'RTMFP   fld=', adjustl(label), ' xor=0x', xorglb(1), &
            ' hash=0x', xorglb(2), ' nnz=', nnzglb(1)
    endif
  end subroutine rtm_fp_int

!-----------------------------------------------------------------------

  subroutine rtm_fp_int2(label, iarr, gidx, cycle)
    ! Emit a bitwise fingerprint of a 2-d INTEGER array laid out as
    ! (slot, cell) -- specifically WRMUnit%myDam(mdn, begr:endr), the list of
    ! dams each gridcell depends on.
    !
    ! This is the single most interesting uninstrumented quantity in the
    ! investigation. myDam is built in WRM_subw_IO_mod.F90:307-335 by scanning
    ! gridID_from_Dam and matching against rtmCTL%gindex, appending each hit in
    ! discovery order. The ORDER of dams within a gridcell's list determines the
    ! order in which uptake is attributed in all three cases of
    ! ExtractionRegulatedFlow, and the array is resized in place as it grows.
    !
    ! The hash therefore mixes BOTH the slot position k and the global cell
    ! index, so a list that holds the same dams in a different order produces a
    ! different hash. Plain XOR alone would be blind to exactly that.
    character(len=*), intent(in) :: label
    integer, intent(in) :: iarr(:,:)
    integer, intent(in) :: gidx(:)
    integer, optional, intent(in) :: cycle

    integer     :: i, k, ier, mycycle
    integer(i8) :: v, xorloc(2), nnzloc(1), xorglb(2), nnzglb(1)

    if (.not. rtm_fp_on) return

    mycycle = 0
    if (present(cycle)) mycycle = cycle

    if (size(iarr,2) /= size(gidx)) then
       if (masterproc) write(iulog,'(a,a)') &
            'RTMFP WARNING: int2 size mismatch, skipping fld=', trim(label)
       return
    endif

    xorloc = 0_i8
    nnzloc = 0_i8
    do i = 1, size(iarr,2)
       do k = 1, size(iarr,1)
          v         = int(iarr(k,i), i8)
          xorloc(1) = ieor(xorloc(1), v)
          ! Mix slot AND global cell index: catches reordering within a cell's
          ! list, which is the whole point of this probe.
          xorloc(2) = ieor(xorloc(2), v * (2_i8 * int(k, i8) + 1_i8) &
                                        * (2_i8 * int(gidx(i), i8) + 1_i8))
          if (iarr(k,i) /= 0) nnzloc(1) = nnzloc(1) + 1_i8
       end do
    end do

    call mpi_allreduce(xorloc, xorglb, 2, MPI_INTEGER8, MPI_BXOR, mpicom_rof, ier)
    call mpi_allreduce(nnzloc, nnzglb, 1, MPI_INTEGER8, MPI_SUM,  mpicom_rof, ier)

    if (masterproc) then
       call rtm_fp_hdr(label, mycycle)
       write(iulog,'(a,a24,a,z16.16,a,z16.16,a,i10,a,i10)') &
            'RTMFP   fld=', adjustl(label), ' xor=0x', xorglb(1), &
            ' hash=0x', xorglb(2), ' nnz=', nnzglb(1), ' nslot=', &
            int(size(iarr,1), i8)
    endif
  end subroutine rtm_fp_int2

!-----------------------------------------------------------------------

  subroutine rtm_fp_dump(label, arr)
    ! Dump per-cell (global index, raw bit pattern, value) to a side file so
    ! the divergent CELL can be named, not just the divergent array.
    !
    ! The global fingerprints answer "which array, which timestep". They
    ! cannot answer "which cell" or "how many cells", because they are
    ! aggregates. This does.
    !
    ! Each rank writes its own file (no MPI-IO, no collective ordering), and
    ! only within [rtm_fp_dump_lo, rtm_fp_dump_hi]. Cells whose value is
    ! exactly zero are skipped -- they are the large majority and cannot
    ! differ. Post-process with:
    !
    !   cat fpdump.<leg>.<nstep>.<label>.*.txt | sort -k1,1n > leg.txt
    !   diff base.txt rest.txt
    !
    ! which yields the offending global indices directly.
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: arr(:)

    integer            :: i, nstep, unitn
    integer(i8)        :: bits
    character(len=256) :: fname
    character(len=8)   :: leg

    if (.not. rtm_fp_on) return

    nstep = get_nstep()
    if (fp_nstep0 < 0) fp_nstep0 = nstep     ! first nstep this leg ever sees

    if (nstep < rtm_fp_dump_lo .or. nstep > rtm_fp_dump_hi) return
    if (size(arr) /= size(rtmCTL%gindex)) then
       if (masterproc) write(iulog,'(a,a)') &
            'RTMFP WARNING: dump size mismatch, skipping fld=', trim(label)
       return
    endif

    ! Both ERS legs share one RUNDIR; tag by the leg's starting nstep so the
    ! base (starts at 0) and rest (starts at the restart nstep) never collide.
    write(leg,'(a,i5.5)') 'n', fp_nstep0
    write(fname,'(a,a,a,i5.5,a,a,a,i5.5,a)') &
         'fpdump.', trim(leg), '.', nstep, '.', trim(label), '.', iam, '.txt'

    open(newunit=unitn, file=trim(fname), status='unknown', &
         position='append', action='write')
    do i = 1, size(arr)
       if (arr(i) == 0._r8) cycle
       bits = transfer(arr(i), bits)
       ! gindex first so a plain numeric sort collates across ranks. Only the
       ! raw bit pattern is written -- it is the exact value, and a decimal
       ! copy alongside it would nearly double the file size for nothing.
       ! Recover the value afterwards only for the few cells that differ.
       write(unitn,'(i10,1x,z16.16)') rtmCTL%gindex(i), bits
    end do
    close(unitn)
  end subroutine rtm_fp_dump

!-----------------------------------------------------------------------

  subroutine rtm_fp_gather(label, arr, iter)
    ! Gather the RAW bits of a rank-local PARTIAL global-indexed array from
    ! every rank onto masterproc and write them to a side file, so the two ERS
    ! legs can be compared exactly rather than by checksum.
    !
    ! Why this exists: rtm_fp_part's three accumulators (plain XOR, index-mixed,
    ! rank-mixed) share a blind spot -- a cross-rank index swap leaves all three
    ! unchanged. Since the observed signature IS a conserved total with moved
    ! attribution, a clean rtm_fp_part does not exclude a permuted input. Raw
    ! bits from all ranks do.
    !
    ! Fires at exactly ONE (nstep, iter) event, the FIRST time it is reached,
    ! because the cost is npes x size(arr) x 8 bytes (~4.3 MB at 128 ranks and
    ! 4246 dams) and there are 3 Euler calls per nstep.
    !
    ! File format, one line per (rank, index) with a nonzero value:
    !     <rank> <index> <raw bits, hex>
    ! Post-process with:
    !     sort -k1,1n -k2,2n fpgather.<leg>.<label>.txt > leg.txt
    !     diff base.txt rest.txt
    ! which names the rank AND the dam index directly.
    character(len=*), intent(in) :: label
    real(r8), intent(in) :: arr(:)
    integer, intent(in)  :: iter

    integer     :: i, ip, n, ier, nstep, unitn
    integer(i8) :: bits
    integer(i8), allocatable :: sendbuf(:), recvbuf(:)
    logical, save :: done = .false.
    character(len=256) :: fname
    character(len=8)   :: leg

    if (.not. rtm_fp_on) return
    if (done) return

    nstep = get_nstep()
    if (fp_nstep0 < 0) fp_nstep0 = nstep     ! first nstep this leg ever sees

    if (nstep /= rtm_fp_gather_nstep) return
    if (iter  /= rtm_fp_gather_iter)  return
    done = .true.                            ! one event only, ever

    n = size(arr)
    allocate(sendbuf(n))
    do i = 1, n
       sendbuf(i) = transfer(arr(i), bits)
    end do

    ! Only masterproc needs a receive buffer, but allocating it everywhere
    ! keeps the call simple and costs nothing on the ranks that pass a
    ! zero-length slice.
    if (masterproc) then
       allocate(recvbuf(n*npes))
    else
       allocate(recvbuf(1))
    endif

    call mpi_gather(sendbuf, n, MPI_INTEGER8, recvbuf, n, MPI_INTEGER8, &
                    0, mpicom_rof, ier)

    if (masterproc) then
       write(leg,'(a,i5.5)') 'n', fp_nstep0
       write(fname,'(a,a,a,a,a,i2.2,a,i5.5,a)') &
            'fpgather.', trim(leg), '.', trim(label), '.iter', iter, &
            '.nstep', nstep, '.txt'
       open(newunit=unitn, file=trim(fname), status='unknown', &
            action='write')
       do ip = 0, npes-1
          do i = 1, n
             if (recvbuf(ip*n + i) == 0_i8) cycle   ! zeros are the vast majority
             write(unitn,'(i6,1x,i8,1x,z16.16)') ip, i, recvbuf(ip*n + i)
          end do
       end do
       close(unitn)
       write(iulog,'(a,a,a,a)') 'RTMFP   gather fld=', trim(label), &
            ' written to ', trim(fname)
    endif

    deallocate(sendbuf)
    deallocate(recvbuf)
  end subroutine rtm_fp_gather

!-----------------------------------------------------------------------

  subroutine rtm_fp_tag(label, cycle)
    ! Emit a bare marker so the two logs can be aligned even where no field
    ! is printed (e.g. "restart file just read").
    character(len=*), intent(in) :: label
    integer, optional, intent(in) :: cycle
    integer :: mycycle

    if (.not. rtm_fp_on) return
    mycycle = 0
    if (present(cycle)) mycycle = cycle
    if (masterproc) then
       call rtm_fp_hdr(label, mycycle)
       write(iulog,'(a,a)') 'RTMFP   mark=', trim(label)
    endif
  end subroutine rtm_fp_tag

end module RtmFingerprint
