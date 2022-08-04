module scream_cpl_indices

  use iso_c_binding, only: c_int, c_char, c_bool

  implicit none
  private

  ! Focus only on the ones that scream imports/exports (subsets of x2a and a2x)
  integer, parameter, public :: num_scream_imports       = 9
  integer, parameter, public :: num_scream_exports       = 16
  integer, public :: num_cpl_imports, num_cpl_exports

  ! cpl indices for import/export fields
  integer(kind=c_int), public, allocatable, target :: index_x2a(:)
  integer(kind=c_int), public, allocatable, target :: index_a2x(:)

  ! Names used by scream for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_x2a(:)

  ! Vector component of import/export fields. If not a vector field, set to -1.
  integer(kind=c_int), public, allocatable, target :: vec_comp_a2x(:)
  integer(kind=c_int), public, allocatable, target :: vec_comp_x2a(:)

  ! Stores whether or not field can be exported during initialization
  logical(kind=c_bool), public, allocatable, target :: can_be_exported_during_init(:)

  public :: scream_set_cpl_indices

  contains

  subroutine scream_set_cpl_indices (x2a, a2x)
    use iso_c_binding,  only: C_NULL_CHAR
    use mct_mod,        only: mct_aVect, mct_avect_indexra

    !
    ! Input(s)
    !
    type(mct_avect), intent(in) :: x2a, a2x
    !
    ! Local(s)
    !
    integer :: i,idx

    ! total number of import/export var in data
    num_cpl_imports = size(x2a%rAttr,1)
    num_cpl_exports = size(a2x%rAttr,1)

    allocate (index_x2a(num_cpl_imports))
    allocate (index_a2x(num_cpl_exports))

    allocate (scr_names_x2a(num_cpl_imports))
    allocate (scr_names_a2x(num_cpl_exports))

    allocate (vec_comp_x2a(num_cpl_imports))
    allocate (vec_comp_a2x(num_cpl_exports))

    allocate (can_be_exported_during_init(num_cpl_exports))

    ! Initialize arrays
    do i=1,num_cpl_imports
      index_x2a(i) = i
      scr_names_x2a(i) = 'unused'
      vec_comp_x2a(i) = -1
    enddo
    do i=1,num_cpl_exports
      index_a2x(i) = i
      scr_names_a2x(i) = 'set_zero'
      vec_comp_a2x(i) = -1
      can_be_exported_during_init(i) = .true.
    enddo

    ! The following are imported to SCREAM
    scr_names_x2a(mct_avect_indexra(x2a,'Sx_avsdr'))  = 'sfc_alb_dir_vis'
    scr_names_x2a(mct_avect_indexra(x2a,'Sx_anidr'))  = 'sfc_alb_dir_nir'
    scr_names_x2a(mct_avect_indexra(x2a,'Sx_avsdf'))  = 'sfc_alb_dif_vis'
    scr_names_x2a(mct_avect_indexra(x2a,'Sx_anidf'))  = 'sfc_alb_dif_nir'
    scr_names_x2a(mct_avect_indexra(x2a,'Faxx_taux')) = 'surf_mom_flux'
    scr_names_x2a(mct_avect_indexra(x2a,'Faxx_tauy')) = 'surf_mom_flux'
    scr_names_x2a(mct_avect_indexra(x2a,'Faxx_sen'))  = 'surf_sens_flux'
    scr_names_x2a(mct_avect_indexra(x2a,'Faxx_evap')) = 'surf_evap'
    scr_names_x2a(mct_avect_indexra(x2a,'Faxx_lwup')) = 'surf_lw_flux_up'

    vec_comp_x2a(mct_avect_indexra(x2a,'Faxx_taux')) = 0
    vec_comp_x2a(mct_avect_indexra(x2a,'Faxx_tauy')) = 1

    ! The following are exported from SCREAM
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_z'))       = 'Sa_z'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_u'))       = 'horiz_winds'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_v'))       = 'horiz_winds'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_tbot'))    = 'T_mid'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_ptem'))    = 'Sa_ptem'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_pbot'))    = 'p_mid'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_shum'))    = 'qv'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_dens'))    = 'Sa_dens'
    scr_names_a2x(mct_avect_indexra(a2x,'Sa_pslv'))    = 'Sa_pslv'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_rainl')) = 'Faxa_rainl'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_snowl')) = 'Faxa_snowl'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_swndr')) = 'sfc_flux_dir_nir'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_swvdr')) = 'sfc_flux_dir_vis'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_swndf')) = 'sfc_flux_dif_nir'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_swvdf')) = 'sfc_flux_dif_vis'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_swnet')) = 'sfc_flux_sw_net'
    scr_names_a2x(mct_avect_indexra(a2x,'Faxa_lwdn'))  = 'sfc_flux_lw_dn'

    vec_comp_a2x(mct_avect_indexra(a2x,'Sa_u')) = 0
    vec_comp_a2x(mct_avect_indexra(a2x,'Sa_v')) = 1

    ! SCREAM needs to run before this field can be exported as they
    ! rely on fields computed by SCREAM
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_rainl')) = .false.
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_snowl')) = .false.
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_swndr')) = .false.
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_swvdr')) = .false.
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_swndf')) = .false.
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_swvdf')) = .false.
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_swnet')) = .false.
    can_be_exported_during_init(mct_avect_indexra(a2x,'Faxa_lwdn'))  = .false.

    ! Trim names
    do i=1,num_cpl_imports
      scr_names_x2a(i) = TRIM(scr_names_x2a(i)) // C_NULL_CHAR
    enddo

    do i=1,num_cpl_exports
      scr_names_a2x(i) = TRIM(scr_names_a2x(i)) // C_NULL_CHAR
    enddo

  end subroutine scream_set_cpl_indices

end module scream_cpl_indices
