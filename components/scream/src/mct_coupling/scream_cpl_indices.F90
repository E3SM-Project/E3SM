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

  character(len=32,kind=c_char), public, allocatable, target :: names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: names_x2a(:)

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
    allocate (names_x2a(num_imports))
    allocate (index_a2x(num_exports))
    allocate (names_a2x(num_exports))

    ! Determine attribute vector indices

    ! List of inputs that scream cares about
    names_x2a(1)  = 'Faxx_evap'
    names_x2a(2)  = 'Faxx_sen'
    names_x2a(3)  = 'Faxx_lat'
    names_x2a(4)  = 'Faxx_taux'
    names_x2a(5)  = 'Faxx_tauy'
    names_x2a(6)  = 'Faxx_lwup'
    names_x2a(7)  = 'Sx_avsdr'
    names_x2a(8)  = 'Sx_anidr'
    names_x2a(9)  = 'Sx_avsdf'
    names_x2a(10) = 'Sx_anidf'
    names_x2a(11) = 'Sx_t'
    names_x2a(12) = 'Sl_snowh'
    names_x2a(13) = 'Si_snowh'
    names_x2a(14) = 'Sx_tref'
    names_x2a(15) = 'Sx_qref'
    names_x2a(16) = 'Sx_u10'
    names_x2a(17) = 'Sf_ifrac'
    names_x2a(18) = 'Sf_ofrac'
    names_x2a(19) = 'Sf_lfrac'
    names_x2a(20) = 'So_ustar'
    names_x2a(21) = 'So_re'

    ! Outputs that scream will pass back to cpl
    names_a2x(1)  = 'Sa_tbot'
    names_a2x(2)  = 'Sa_ptem'
    names_a2x(3)  = 'Sa_z'
    names_a2x(4)  = 'Sa_u'
    names_a2x(5)  = 'Sa_v'
    names_a2x(6)  = 'Sa_pbot'
    names_a2x(7)  = 'Sa_dens'
    names_a2x(8)  = 'Sa_shum'
    names_a2x(9)  = 'Faxa_rainc'
    names_a2x(10) = 'Faxa_rainl'
    names_a2x(11) = 'Faxa_snowc'
    names_a2x(12) = 'Faxa_snowl'
    names_a2x(13)  = 'Sa_co2prog'

    do i=1,num_required_imports
      index_x2a(i) = mct_avect_indexra(x2a,TRIM(names_x2a(i)))
      names_x2a(i) = TRIM(names_x2a(i)) // C_NULL_CHAR
    enddo
    do i=num_required_imports+1,num_imports
      index_x2a(i) = mct_avect_indexra(x2a,TRIM(names_x2a(i)),perrWith='quiet')
      names_x2a(i) = TRIM(names_x2a(i)) // C_NULL_CHAR
    enddo

    do i=1,num_required_exports
      index_a2x(i) = mct_avect_indexra(a2x,TRIM(names_a2x(i)))
      names_a2x(i) = TRIM(names_a2x(i)) // C_NULL_CHAR
    enddo
    do i=num_required_exports+1,num_exports
      index_a2x(i) = mct_avect_indexra(a2x,TRIM(names_a2x(i)),perrWith='quiet')
      names_a2x(i) = TRIM(names_a2x(i)) // C_NULL_CHAR
    enddo
  end subroutine scream_set_cpl_indices

end module scream_cpl_indices
