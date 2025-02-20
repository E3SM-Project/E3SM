module eamxx_cpl_indices

  use iso_c_binding, only: c_int, c_double, c_char, c_bool

  implicit none
  private

  ! Focus only on the ones that scream imports/exports (subsets of x2a and a2x)
  integer, parameter, public :: num_scream_imports = 24
  integer, parameter, public :: num_scream_exports = 17
  integer, public :: num_cpl_imports, num_cpl_exports, import_field_size, export_field_size

  ! Names used by scream for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: import_field_names(:)
  character(len=32,kind=c_char), public, allocatable, target :: export_field_names(:)

  ! cpl indices for import/export fields
  integer(kind=c_int), public, allocatable, target :: import_cpl_indices(:)
  integer(kind=c_int), public, allocatable, target :: export_cpl_indices(:)

  ! Vector component of import/export fields. If not a vector field, set to -1.
  integer(kind=c_int), public, allocatable, target :: import_vector_components(:)
  integer(kind=c_int), public, allocatable, target :: export_vector_components(:)

  ! Constant multiple to apply to import/export fields. Currently used for fluxes
  ! which have different interpretation of direction between cpl and scream
  real(kind=c_double), public, allocatable, target :: import_constant_multiple(:)
  real(kind=c_double), public, allocatable, target :: export_constant_multiple(:)

  ! Stores whether or not field should/can be imported/exported during initialization
  logical(kind=c_bool), public, allocatable, target :: do_import_during_init(:)
  logical(kind=c_bool), public, allocatable, target :: do_export_during_init(:)

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
    integer :: i

    ! total number of import/export fields in data
    num_cpl_imports = size(x2a%rAttr,1)
    num_cpl_exports = size(a2x%rAttr,1)

    ! Size of import and export fields data
    import_field_size = size(x2a%rAttr,2)
    export_field_size = size(a2x%rAttr,2)

    ! IMPORT

    ! Import info is of size num_scream_imports. Any additional imports are not touched.
    allocate (import_field_names(num_scream_imports))
    allocate (import_cpl_indices(num_scream_imports))
    allocate (import_vector_components(num_scream_imports))
    allocate (import_constant_multiple(num_scream_imports))
    allocate (do_import_during_init(num_scream_imports))

    ! Initialize arrays
    do i=1,num_scream_imports
      import_vector_components(i) = -1
      import_constant_multiple(i) = 1
      do_import_during_init(i) = .false.
    enddo

    ! SCREAM names
    import_field_names(1)  = 'sfc_alb_dir_vis'
    import_field_names(2)  = 'sfc_alb_dir_nir'
    import_field_names(3)  = 'sfc_alb_dif_vis'
    import_field_names(4)  = 'sfc_alb_dif_nir'
    import_field_names(5)  = 'surf_radiative_T'
    import_field_names(6)  = 'T_2m'
    import_field_names(7)  = 'qv_2m'
    import_field_names(8)  = 'wind_speed_10m'
    import_field_names(9)  = 'snow_depth_land'
    import_field_names(10) = 'surf_lw_flux_up'
    import_field_names(11) = 'surf_mom_flux'
    import_field_names(12) = 'surf_mom_flux'
    import_field_names(13) = 'surf_sens_flux'
    import_field_names(14) = 'surf_evap'
    import_field_names(15) = 'ocnfrac'
    import_field_names(16) = 'landfrac'
    import_field_names(17) = 'icefrac'
    import_field_names(18) = 'fv'
    import_field_names(19) = 'ram1'
    import_field_names(20) = 'sst'
    import_field_names(21) = 'dstflx'
    import_field_names(22) = 'dstflx'
    import_field_names(23) = 'dstflx'
    import_field_names(24) = 'dstflx'

    ! CPL indices
    import_cpl_indices(1)  = mct_avect_indexra(x2a,'Sx_avsdr')
    import_cpl_indices(2)  = mct_avect_indexra(x2a,'Sx_anidr')
    import_cpl_indices(3)  = mct_avect_indexra(x2a,'Sx_avsdf')
    import_cpl_indices(4)  = mct_avect_indexra(x2a,'Sx_anidf')
    import_cpl_indices(5)  = mct_avect_indexra(x2a,'Sx_t')
    import_cpl_indices(6)  = mct_avect_indexra(x2a,'Sx_tref')
    import_cpl_indices(7)  = mct_avect_indexra(x2a,'Sx_qref')
    import_cpl_indices(8)  = mct_avect_indexra(x2a,'Sx_u10')
    import_cpl_indices(9)  = mct_avect_indexra(x2a,'Sl_snowh')
    import_cpl_indices(10) = mct_avect_indexra(x2a,'Faxx_lwup')
    import_cpl_indices(11) = mct_avect_indexra(x2a,'Faxx_taux')
    import_cpl_indices(12) = mct_avect_indexra(x2a,'Faxx_tauy')
    import_cpl_indices(13) = mct_avect_indexra(x2a,'Faxx_sen')
    import_cpl_indices(14) = mct_avect_indexra(x2a,'Faxx_evap')
    import_cpl_indices(15) = mct_avect_indexra(x2a,'Sf_ofrac')
    import_cpl_indices(16) = mct_avect_indexra(x2a,'Sf_lfrac')
    import_cpl_indices(17) = mct_avect_indexra(x2a,'Sf_ifrac') 
    import_cpl_indices(18) = mct_avect_indexra(x2a,'Sl_fv')
    import_cpl_indices(19) = mct_avect_indexra(x2a,'Sl_ram1')
    !sst
    import_cpl_indices(20) = mct_avect_indexra(x2a,'So_t')
    !dust fluxes
    import_cpl_indices(21) = mct_avect_indexra(x2a,'Fall_flxdst1')
    import_cpl_indices(22) = mct_avect_indexra(x2a,'Fall_flxdst2')
    import_cpl_indices(23) = mct_avect_indexra(x2a,'Fall_flxdst3')
    import_cpl_indices(24) = mct_avect_indexra(x2a,'Fall_flxdst4')

    ! Vector components
    !(Faxx_taux and Faxx_tauy)
    import_vector_components(11) = 0
    import_vector_components(12) = 1
    !(dust fluxes)
    import_vector_components(21) = 0
    import_vector_components(22) = 1
    import_vector_components(23) = 2
    import_vector_components(24) = 3

    ! Constant multiples
    import_constant_multiple(10) = -1
    import_constant_multiple(11) = -1
    import_constant_multiple(12) = -1
    import_constant_multiple(13) = -1
    import_constant_multiple(14) = -1

    ! Does this field need to be imported during intialization
    do_import_during_init(11) = .true.
    do_import_during_init(12) = .true.
    do_import_during_init(13) = .true.
    do_import_during_init(14) = .true.

    ! EXPORT

    ! Export info is of size num_cpl_imports. If SCREAM does not export a field, it is set to zero.
    allocate (export_field_names(num_scream_exports))
    allocate (export_cpl_indices(num_scream_exports))
    allocate (export_vector_components(num_scream_exports))
    allocate (export_constant_multiple(num_scream_exports))
    allocate (do_export_during_init(num_scream_exports))

    ! Initialize arrays
    do i=1,num_scream_exports
      export_vector_components(i) = -1
      export_constant_multiple(i) = 1
      do_export_during_init(i) = .false.
    enddo

    ! SCREAM names
    export_field_names(1)  = 'Sa_z'
    export_field_names(2)  = 'Sa_u'
    export_field_names(3)  = 'Sa_v'
    export_field_names(4)  = 'Sa_tbot' 
    export_field_names(5)  = 'Sa_ptem'
    export_field_names(6)  = 'Sa_pbot'
    export_field_names(7)  = 'Sa_shum'
    export_field_names(8)  = 'Sa_dens'
    export_field_names(9)  = 'Sa_pslv'
    export_field_names(10) = 'Faxa_rainl'
    export_field_names(11) = 'Faxa_snowl'
    export_field_names(12) = 'Faxa_swndr' 
    export_field_names(13) = 'Faxa_swvdr' 
    export_field_names(14) = 'Faxa_swndf' 
    export_field_names(15) = 'Faxa_swvdf' 
    export_field_names(16) = 'Faxa_swnet' 
    export_field_names(17) = 'Faxa_lwdn'  

    ! CPL indices
    export_cpl_indices(1)  = mct_avect_indexra(a2x,'Sa_z')
    export_cpl_indices(2)  = mct_avect_indexra(a2x,'Sa_u')
    export_cpl_indices(3)  = mct_avect_indexra(a2x,'Sa_v')
    export_cpl_indices(4)  = mct_avect_indexra(a2x,'Sa_tbot')
    export_cpl_indices(5)  = mct_avect_indexra(a2x,'Sa_ptem')
    export_cpl_indices(6)  = mct_avect_indexra(a2x,'Sa_pbot')
    export_cpl_indices(7)  = mct_avect_indexra(a2x,'Sa_shum')
    export_cpl_indices(8)  = mct_avect_indexra(a2x,'Sa_dens')
    export_cpl_indices(9)  = mct_avect_indexra(a2x,'Sa_pslv')
    export_cpl_indices(10) = mct_avect_indexra(a2x,'Faxa_rainl')
    export_cpl_indices(11) = mct_avect_indexra(a2x,'Faxa_snowl')
    export_cpl_indices(12) = mct_avect_indexra(a2x,'Faxa_swndr')
    export_cpl_indices(13) = mct_avect_indexra(a2x,'Faxa_swvdr')
    export_cpl_indices(14) = mct_avect_indexra(a2x,'Faxa_swndf')
    export_cpl_indices(15) = mct_avect_indexra(a2x,'Faxa_swvdf')
    export_cpl_indices(16) = mct_avect_indexra(a2x,'Faxa_swnet')
    export_cpl_indices(17) = mct_avect_indexra(a2x,'Faxa_lwdn')

    ! Does this field need to be imported during intialization
    do_export_during_init(1) = .true.
    do_export_during_init(2) = .true.
    do_export_during_init(3) = .true.
    do_export_during_init(4) = .true.
    do_export_during_init(5) = .true.
    do_export_during_init(6) = .true.
    do_export_during_init(7) = .true.
    do_export_during_init(8) = .true.
    do_export_during_init(9) = .true.

    ! Trim names
    do i=1,num_scream_imports
      import_field_names(i) = TRIM(import_field_names(i)) // C_NULL_CHAR
    enddo

    do i=1,num_scream_exports
      export_field_names(i) = TRIM(export_field_names(i)) // C_NULL_CHAR
    enddo

  end subroutine scream_set_cpl_indices

end module eamxx_cpl_indices
