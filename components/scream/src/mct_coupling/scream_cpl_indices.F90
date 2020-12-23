module scream_cpl_indices

  use iso_c_binding, only: c_int, c_char

  implicit none
  private

  ! Focus only on the ones that scream imports/exports (subsets of x2a and a2x)
  integer, parameter, public :: num_required_imports = 21
  integer, parameter, public :: num_required_exports = 12
  integer, parameter, public :: num_optional_imports = 0
  integer, parameter, public :: num_optional_exports = 1
  integer, parameter, public :: num_imports = num_required_imports + num_optional_imports
  integer, parameter, public :: num_exports = num_required_exports + num_optional_exports

  integer(kind=c_int), public, allocatable, target :: index_x2a(:)
  integer(kind=c_int), public, allocatable, target :: index_a2x(:)

  ! Names used by the component coupler for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: cpl_names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: cpl_names_x2a(:)

  ! Names used by scream for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_x2a(:)

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

    allocate (index_x2a(num_imports))
    allocate (index_a2x(num_exports))

    allocate (cpl_names_x2a(num_imports))
    allocate (cpl_names_a2x(num_exports))

    allocate (scr_names_x2a(num_imports))
    allocate (scr_names_a2x(num_exports))

    ! Determine attribute vector indices

    ! List of cpl names of inputs that scream cares about
    cpl_names_x2a(1)  = 'Faxx_evap' ! SHOC input
    cpl_names_x2a(2)  = 'Faxx_sen'  ! SHOC input
    cpl_names_x2a(3)  = 'Faxx_lat'  ! Needed for energy fixer (?)
    cpl_names_x2a(4)  = 'Faxx_taux' ! SHOC input
    cpl_names_x2a(5)  = 'Faxx_tauy' ! SHOC input
    cpl_names_x2a(6)  = 'Faxx_lwup' ! RRTMGP input
    cpl_names_x2a(7)  = 'Sx_avsdr'  ! RRTMGP input
    cpl_names_x2a(8)  = 'Sx_anidr'  ! RRTMGP input
    cpl_names_x2a(9)  = 'Sx_avsdf'  ! RRTMGP input
    cpl_names_x2a(10) = 'Sx_anidf'  ! RRTMGP input
    cpl_names_x2a(11) = 'Sx_t'      ! SHOC input (?)
    cpl_names_x2a(12) = 'Sl_snowh'  ! SHOC input (?)
    cpl_names_x2a(13) = 'Si_snowh'  ! UNUSED
    cpl_names_x2a(14) = 'Sx_tref'   ! UNUSED
    cpl_names_x2a(15) = 'Sx_qref'   ! UNUSED
    cpl_names_x2a(16) = 'Sx_u10'    ! UNUSED
    cpl_names_x2a(17) = 'Sf_ifrac'  ! RRTMGP input
    cpl_names_x2a(18) = 'Sf_ofrac'  ! SHOCK input (?)
    cpl_names_x2a(19) = 'Sf_lfrac'  ! RRTMGP input
    cpl_names_x2a(20) = 'So_ustar'  ! SHOC (?): doesn't seem to be cam_in%ustar though
    cpl_names_x2a(21) = 'So_re'     ! UNUSED

    ! Names used by scream for the input fields above
    scr_names_x2a(1)  = 'surface_water_evaporation_flux'
    scr_names_x2a(2)  = 'Faxx_sen'
    scr_names_x2a(3)  = 'Faxx_lat'
    scr_names_x2a(4)  = 'Faxx_taux'
    scr_names_x2a(5)  = 'Faxx_tauy'
    scr_names_x2a(6)  = 'Faxx_lwup'
    scr_names_x2a(7)  = 'Sx_avsdr'
    scr_names_x2a(8)  = 'Sx_anidr'
    scr_names_x2a(9)  = 'Sx_avsdf'
    scr_names_x2a(10) = 'Sx_anidf'
    scr_names_x2a(11) = 'Sx_t'
    scr_names_x2a(12) = 'Sl_snowh'
    scr_names_x2a(13) = 'Si_snowh'
    scr_names_x2a(14) = 'Sx_tref'
    scr_names_x2a(15) = 'Sx_qref'
    scr_names_x2a(16) = 'Sx_u10'
    scr_names_x2a(17) = 'Sf_ifrac'
    scr_names_x2a(18) = 'Sf_ofrac'
    scr_names_x2a(19) = 'Sf_lfrac'
    scr_names_x2a(20) = 'So_ustar'
    scr_names_x2a(21) = 'So_re'

    ! List of cpl names of outputs that scream needs to pass back to cpl
    cpl_names_a2x(1)  = 'Sa_tbot'     ! Temperature (HOMME)
    cpl_names_a2x(2)  = 'Sa_ptem'     ! Potential temperature (HOMME)
    cpl_names_a2x(3)  = 'Sa_z'        ! Z coord of bottom layers (in m? Pa?)
    cpl_names_a2x(4)  = 'Sa_u'        ! Horiz velocity u-comp (HOMME)
    cpl_names_a2x(5)  = 'Sa_v'        ! Horiz velocity v-comp (HOMME)
    cpl_names_a2x(6)  = 'Sa_pbot'     ! Pressure (HOMME)
    cpl_names_a2x(7)  = 'Sa_dens'     ! Density (derived from P, T)
    cpl_names_a2x(8)  = 'Sa_shum'     ! Specific humidity (Computed by ...?)
    cpl_names_a2x(9)  = 'Faxa_rainc'  ! Liquid convective precip. P3 (?)
    cpl_names_a2x(10) = 'Faxa_rainl'  ! Liquire large-scale precip. P3 (?)
    cpl_names_a2x(11) = 'Faxa_snowc'  ! Frozen convective precip. P3 (?)
    cpl_names_a2x(12) = 'Faxa_snowl'  ! Frozen large scale precip. P3 (?)
    cpl_names_a2x(13)  = 'Sa_co2prog' ! NOT COMPUTED by SCREAM (?)

    ! Names used by scream for the output fields above
    scr_names_a2x(1)  = 'surface_temperature'
    scr_names_a2x(2)  = 'Sa_ptem'
    scr_names_a2x(3)  = 'Sa_z'
    scr_names_a2x(4)  = 'Sa_u'
    scr_names_a2x(5)  = 'Sa_v'
    scr_names_a2x(6)  = 'Sa_pbot'
    scr_names_a2x(7)  = 'Sa_dens'
    scr_names_a2x(8)  = 'Sa_shum'
    scr_names_a2x(9)  = 'Faxa_rainc'
    scr_names_a2x(10) = 'Faxa_rainl'
    scr_names_a2x(11) = 'Faxa_snowc'
    scr_names_a2x(12) = 'Faxa_snowl'
    scr_names_a2x(13)  = 'Sa_co2prog'

    do i=1,num_required_imports
      index_x2a(i) = mct_avect_indexra(x2a,TRIM(cpl_names_x2a(i)))
      scr_names_x2a(i) = TRIM(scr_names_x2a(i)) // C_NULL_CHAR
    enddo
    do i=num_required_imports+1,num_imports
      index_x2a(i) = mct_avect_indexra(x2a,TRIM(cpl_names_x2a(i)),perrWith='quiet')
      scr_names_x2a(i) = TRIM(scr_names_x2a(i)) // C_NULL_CHAR
    enddo

    do i=1,num_required_exports
      index_a2x(i) = mct_avect_indexra(a2x,TRIM(cpl_names_a2x(i)))
      scr_names_a2x(i) = TRIM(scr_names_a2x(i)) // C_NULL_CHAR
    enddo
    do i=num_required_exports+1,num_exports
      index_a2x(i) = mct_avect_indexra(a2x,TRIM(cpl_names_a2x(i)),perrWith='quiet')
      scr_names_a2x(i) = TRIM(scr_names_a2x(i)) // C_NULL_CHAR
    enddo

    ! We no longer need the cpl names
    deallocate(cpl_names_a2x)
    deallocate(cpl_names_x2a)

  end subroutine scream_set_cpl_indices

end module scream_cpl_indices
