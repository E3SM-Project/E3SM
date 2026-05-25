module seq_nlmap_mod

  !-----------------------------------------------------------------------------
  !
  ! Purpose: Apply a high-order linear map followed by a property restoring
  ! nonlinear global fixer.  The objective is to improve
  ! coarse-to-fine, conservative, monotone maps. Area-averaged (aave) maps
  ! produce grid imprint on the target fine grid because each coarse source cell
  ! has a constant field. This module applies a high-order linear map that
  ! produces smooth data on the fine target grid, then applies a nonlinear
  ! global fixer to restore conservation and bounds.
  !
  ! Usage:
  !   For any existing (area-averaged; see the next paragraph) mapfile
  ! SRC2TGT_TMAPFILE, where T is 'F' or 'S', optionally provide a second,
  ! high-order map by specifying SRC2TGT_TMAPFILE_NONLINEAR. If T is 'F', then
  ! the fixer restores global mass, as defined by the SRC2TGT_TMAPFILE map, and
  ! local bounds. If T is 'S', then the fixer restores the bounds.
  !   In practice, SRC2TGT_TMAPFILE should be an area-averaged map ('aave', aka
  ! 'mono'). This map provides the reference global mass on the target grid and
  ! the non-0 target cells. SRC2TGT_TMAPFILE_NONLINEAR can be anything, but its
  ! This requirement assures that SRC2TGT_TMAPFILE_NONLINEAR provides
  ! data in any target cell that SRC2TGT_TMAPFILE does. In the opposite
  ! direction, at runtime, any target cell that SRC2TGT_TMAPFILE does not affect
  ! is zeroed after SRC2TGT_TMAPFILE_NONLINEAR is applied.
  !
  ! Optionally, there is an algorithm variant to obtains exact
  ! mass conservation for atmosphere-to-surface flux maps in the case of
  ! overlapping surface grids with a nontrivial common refinement. (A unified
  ! surface grid permits conservation without the variant.) This variant is
  ! controlled by the boolean option NLMAPS_ATM2SRF_CONSERVE in the coupler
  ! options. If this option is true, both the atm2lnd and atm2ocn maps must be
  ! nonlinear; if they aren't, the module calls abort.  NOT YET IMPLMENTED IN MOAB
  !
  ! Algorithm is now implementd inside MOAB in src/Remapping/TempestOnlineMap.cpp
  !
  ! Following is a description of the algorithm and its properties,
  ! specializations of Alg. 3.1 of the following reference for this use case:
  !   Bradley, A. M., Bosler, P. A., Guba, O., Taylor, M. A., and Barnett,
  !   G. A. (2019). Communication-efficient property preservation in tracer
  !   transport. SIAM Journal on Scientific Computing, 41(3), C161-C193.
  !
  ! Notation:
  !   A x is the matrix-vector product A times x.
  !   a x is a scalar times a vector.
  !   x'y is the dot product of x and y.
  !   non-zero-pattern(A) is a matrix of size A with 1 where A(i,j) /= 0
  !     and 0 otherwise.
  !   (In)equalities and / are applied element-wise.
  !   In practice:
  !     Am is the area-average conservative and monotone linear map.
  !     A  is the high-order linear map.
  !     f  is a vector of area fractions.
  ! On input:
  !   1. Am is a monotone and conservative linear operator, where s is the
  !      vector of src areas; t, tgt areas.
  !   2. non-zero-pattern(Am) in non-zero-pattern(A).
  !   3. 0 <= f <= 1.
  ! Definition of conservative:
  !   t'Am x = s'x
  ! Definition of monotone:
  !   l, u := bounds(Am, x)
  !     => l(i) := min(x where Am(i,:) /= 0)
  !        u(i) similarly
  !   l <= Am x <= u
  !   => 4. 0 <= Am(i,j) <= 1
  ! With fractions:
  !   g := Am f
  !   ym := (Am (f x)) / g
  !   l, u := bounds(Am, x)
  !   Conservative:
  !     5. t'(g y) = t'(Am (f x)) = s'(f x)
  !   Monotone:
  !     c(j) := Am(i,j) f(j)
  !     0 <= c(j) <= 1 by 3 and 4.
  !     Let j vary over non-0s of Am(i,:).
  !     ym(i) = sum_j Am(i,j) f(j) x(j) / sum_j Am(i,j) f(j)
  !           = sum_j c(j) x(j) / sum_j c(j)
  !          <= max_j x(j) sum_j c(j) / sum_j c(j)
  !           = max_j x(j)
  !           = u(j)
  !     and similarly for l(j).
  !     => 6. l <= ym <= u
  ! Bounds:
  !   2 => 7. bounds(Am, x) in bounds(A, x)
  ! CAAS (Alg. 3.1 in doi:10.1137/18M1165414) applied to this case:
  !   CAAS(Am, A, f, x)
  !     g := Am f
  !     gym := Am (f x)
  !     M := t'gym
  !     y0 := (A (f x)) / g
  !     l, u := bounds(A, x)
  !     zero y0, l, u in any cell in which gym is 0
  !     y1 := max(l, min(u, y0))
  !     dM := M - t'(g y1)
  !     if dM >= 0
  !       w := u - y1
  !     else
  !       w := y1 - l
  !     y2 := y1 + (dM / (t g)'w) w
  !     return y2
  ! This algorithm has a nonempty constraint set because, by 5 and 6,
  !   ym := (Am (f x)) / g is conservative and in bounds.
  ! CAAS finds a solution if the constraint set is nonempty, proving y2 is
  ! conservative and in bounds. We can prove this in detail, as follows:
  !   y2 is conservative:
  !     t'(g y2) = t'(g y1 + (dM / (t g)'w) g w)
  !              = M - dM + (dM / (t g)'w) t'(g w)
  !              = M
  !   y2 is in bounds:
  !     w >= 0 by construction.
  !     Assume dM >= 0; the other case is similar.
  !     8. M <= t'(g u):
  !       M = t'(Am (f x)) = t'(g ym) <= t'(g u) by 6 and 7
  !     dM <= (t g)'w:
  !       (t g)'w = (t g)'u - (t g)'y1 = (t g)'u + dM - M
  !         = dM + [(t g)'u - M] >= dM by 8
  !     => y2 = y1 + (dM / (t g)'w) w <= y1 + w = u => y2 <= u.
  ! In this proof, the crucial part is point 8; it is equivalent to the
  ! statement that the constraint set is nonempty. This in turn makes
  !   0 <= (dM / (t g)'w) <= 1,
  ! permitting the final line to hold.
  !
  ! Conservative variant:  NOT YET IMPLEMENTED IN MOAB
  !   Some modifications are made for this case.
  !   f, the vector of area fractions, is not used. That is because it is atm
  ! fractions, which are uniformly and globally 1.
  !   Instead, we need nontrivial fraction fields. We use the notation of the
  ! documentation in seq_frac_mct.F90.
  !   For atm2lnd, on the atmosphere we use frac_s = fractions_a(lfrac); on the
  ! land we set frac_d to the 'frac' field in dom_cx_d.
  !   For atm2ocn, on the atmosphere we use frac_s = 1 - fractions_a(lfrac); on
  ! the ocn grid we again set frac_d to the 'frac' field in dom_cx_d, which is
  ! the sum of the ocn and ice fractions in the fractions_ox structure.
  !   Then in the above statement of the CAAS algorithm, these lines are used
  !  instead of the original ones:
  !     (1) M := sum( frac_s * area_s * x )
  !     (2) zero y0, l, u in any cell in which frac_d is 0
  !     (3) dM := M - sum( frac_d * area_d * y1 )
  !
  ! Author: A.M. Bradley, Mar,Apr-2023
  ! Update: A.M. Bradley, Mar-2025. a2s_cons feature.
  ! Update: MOAB Team, May 2026. MOAB version modifed from MCT original
  !
  !-----------------------------------------------------------------------------
  
  use shr_kind_mod     , only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN, &
                               SHR_KIND_CS, CX => SHR_KIND_CX
  use shr_sys_mod
  use shr_const_mod
  use mct_mod
  use seq_comm_mct
  use component_type_mod
  use seq_map_type_mod
  use seq_infodata_mod , only: nlmaps_exclude_max_number, nlmaps_exclude_nchar
  use shr_moab_mod     , only: mbGetnCells
  use iMOAB            , only: iMOAB_GetDoubleTagStorage
  use iso_c_binding    , only: C_NULL_CHAR

  implicit none
  save
  private
#include <mpif.h>

  ! All routines are intended to be called on the cpl pes.
  public :: seq_nlmap_setopts
  public :: seq_nlmap_init_a2oi_cons
  public :: seq_nlmap_init_a2l_cons
  public :: seq_nlmap_field_is_excluded

  ! List of fields that are to be excluded from the set that are nonlinearly
  ! mapped. For these excluded fields, the low-order linear map is used,
  ! instead.
  character(nlmaps_exclude_nchar) :: nlmaps_exclude_fields(nlmaps_exclude_max_number)
  integer :: nlmaps_exclude_n_fields
  logical :: atm2srf_conserve

contains

  subroutine seq_nlmap_setopts(nlmaps_exclude_fields_in, &
       nlmaps_atm2srf_conserve_in)
    ! Options in drv_in.

    character(nlmaps_exclude_nchar), optional, intent(in) :: &
         nlmaps_exclude_fields_in(nlmaps_exclude_max_number)
    logical, optional, intent(in) :: nlmaps_atm2srf_conserve_in

    integer :: i, n

    nlmaps_exclude_n_fields = 0
    if (present(nlmaps_exclude_fields_in)) then
       do i = 1, nlmaps_exclude_max_number
          n = len(trim(nlmaps_exclude_fields_in(i)))
          if (n > 0) then
             nlmaps_exclude_n_fields = nlmaps_exclude_n_fields + 1
             nlmaps_exclude_fields(nlmaps_exclude_n_fields) = nlmaps_exclude_fields_in(i)
          end if
       end do
    end if

    if (present(nlmaps_atm2srf_conserve_in)) then
       atm2srf_conserve = nlmaps_atm2srf_conserve_in
    end if

    if (atm2srf_conserve .and. nlmaps_exclude_n_fields > 0) then
       if (seq_comm_iamroot(CPLID)) then
          write(logunit,'(a)') 'nlmap> WARNING: When nlmaps_atm2srf_conserve is ON,&
               & the field exclusion list is ignored.'
       end if
    end if

  end subroutine seq_nlmap_setopts

  function seq_nlmap_field_is_excluded(fldname) result(is_excluded)
    ! True if fldname appears in the namelist nlmaps_exclude_fields list, in
    ! which case the nonlinear map must NOT be applied (the LOW-order map
    ! result is kept). This mirrors the loop in MCT's seq_nlmap_avNormArr at
    ! driver-mct/main/seq_nlmap_mod.F90 lines 691-698: for excluded fields
    ! avp_o keeps the mct_sMat_avMult low-order value rather than being
    ! overwritten with the nonlinear-fixer result.
    !
    ! When nlmaps_atm2srf_conserve is on, the exclude list is ignored
    ! (matches MCT's warning at the end of seq_nlmap_setopts).
    character(len=*), intent(in) :: fldname
    logical :: is_excluded

    integer :: i, n

    is_excluded = .false.
    if (atm2srf_conserve) return
    if (nlmaps_exclude_n_fields <= 0) return

    n = len_trim(fldname)
    if (n <= 0) return

    do i = 1, nlmaps_exclude_n_fields
       if ( trim(fldname) == trim(nlmaps_exclude_fields(i)) ) then
          is_excluded = .true.
          return
       end if
    end do
  end function seq_nlmap_field_is_excluded

  subroutine seq_nlmap_init_a2oi_cons(mapper, fractions_ax)
    ! Initialize frac_s, frac_d for the atm2ocn map.
    
    type(seq_map), pointer, intent(inout) :: mapper ! Fa2o
    type(mct_aVect)       , intent(in)    :: fractions_ax(:)

    integer(IN) :: k_sfrac, lsize_s, lsize_d, j
    integer :: ierr, ent_type, arrsize
    character(len=128) :: tagname

    if (.not. atm2srf_conserve) return

    if (.not. mapper%nl_available) then
       call shr_sys_abort('seq_nlmap_init_a2oi_cons ERROR: Nonlinear map not set.')
    end if

    k_sfrac = mct_aVect_indexRA(fractions_ax(1), 'lfrac')
    lsize_s = mct_aVect_lsize(fractions_ax(1))
    lsize_d = mbGetnCells(mapper%tgt_mbid)

    ! Guard against double-allocate on mapper re-init (would abort otherwise).
    ! mapper points into the static seq_maps(:) array (seq_map_type_mod), so
    ! it has program lifetime; releasing/re-claiming here also avoids a leak
    ! if a caller re-initializes the same mapper with a different source size.
    if (allocated(mapper%frac_s)) deallocate(mapper%frac_s)
    if (allocated(mapper%frac_d)) deallocate(mapper%frac_d)
    allocate(mapper%frac_s(lsize_s), mapper%frac_d(lsize_d))

    do j = 1,lsize_s
       mapper%frac_s(j) = 1 - fractions_ax(1)%rAttr(k_sfrac,j)
    end do
    ent_type = mapper%tag_entity_type
    arrsize = lsize_d
    tagname = 'frac'//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mapper%tgt_mbid, tagname, arrsize, ent_type, mapper%frac_d)
    if (ierr /= 0) then
       call shr_sys_abort('seq_nlmap_init_a2oi_cons ERROR: Failed to get frac tag from target MOAB mesh.')
    end if

  end subroutine seq_nlmap_init_a2oi_cons

  subroutine seq_nlmap_init_a2l_cons(mapper, fractions_ax)
    ! Initialize frac_s, frac_d for the atm2lnd map.

    type(seq_map), pointer, intent(inout) :: mapper ! Fa2l
    type(mct_aVect)       , intent(in)    :: fractions_ax(:)

    integer(IN) :: k_sfrac, lsize_s, lsize_d, j
    integer :: ierr, ent_type, arrsize
    character(len=128) :: tagname

    if (.not. atm2srf_conserve) return

    if (.not. mapper%nl_available) then
       call shr_sys_abort('seq_nlmap_init_a2l_cons ERROR: Nonlinear map not set.')
    end if

    k_sfrac = mct_aVect_indexRA(fractions_ax(1), 'lfrac')
    lsize_s = mct_aVect_lsize(fractions_ax(1))
    lsize_d = mbGetnCells(mapper%tgt_mbid)

    ! Guard against double-allocate on mapper re-init; see comment in
    ! seq_nlmap_init_a2oi_cons.
    if (allocated(mapper%frac_s)) deallocate(mapper%frac_s)
    if (allocated(mapper%frac_d)) deallocate(mapper%frac_d)
    allocate(mapper%frac_s(lsize_s), mapper%frac_d(lsize_d))

    do j = 1,lsize_s
       mapper%frac_s(j) = fractions_ax(1)%rAttr(k_sfrac,j)
    end do
    ent_type = mapper%tag_entity_type
    arrsize = lsize_d
    tagname = 'frac'//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mapper%tgt_mbid, tagname, arrsize, ent_type, mapper%frac_d)
    if (ierr /= 0) then
       call shr_sys_abort('seq_nlmap_init_a2l_cons ERROR: Failed to get frac tag from target MOAB mesh.')
    end if

  end subroutine seq_nlmap_init_a2l_cons

end module seq_nlmap_mod
