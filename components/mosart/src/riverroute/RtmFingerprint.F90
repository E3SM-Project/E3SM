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
  public :: rtm_fp_dump ! dump per-cell values to a side file, to NAME the cell

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
  integer, public :: rtm_fp_dump_lo = 53  ! first nstep to dump (inclusive)
  integer, public :: rtm_fp_dump_hi = 55  ! last  nstep to dump (inclusive)

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
