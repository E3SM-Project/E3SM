module glc_zocnclass_mod

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module contains data and routines for operating on GLC ocean z-level classes.

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod
  use seq_comm_mct, only : logunit
  use shr_log_mod, only : errMsg => shr_log_errMsg

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: glc_zocnclass_init            ! initialize GLC z-ocean class data
  public :: glc_zocnclass_clean           ! deallocate memory allocated here
  public :: glc_get_num_zocn_classes      ! get the number of z-ocean classes
  public :: glc_get_zlevels               ! get an array of the z-ocean levels
  public :: glc_get_zocnclass_bounds      ! get the boundaries of all z-ocean classes
  public :: glc_zocnclass_as_string       ! returns a string corresponding to a given z-ocean class
  public :: glc_all_zocnclass_strings     ! returns an array of strings for all z-ocean classes
  public :: glc_zocn_errcode_to_string    ! convert an error code into a string describing the error

  interface glc_zocnclass_init
     module procedure glc_zocnclass_init_default
     module procedure glc_zocnclass_init_override
  end interface glc_zocnclass_init


  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  ! Possible error code values
  integer, parameter, public :: GLC_ZOCNCLASS_ERR_NONE = 0      ! err_code indicating no error
  integer, parameter, public :: GLC_ZOCNCLASS_ERR_UNDEFINED = 1 ! err_code indicating z-ocean classes have not been defined
  integer, parameter, public :: GLC_ZOCNCLASS_ERR_TOO_LOW = 2   ! err_code indicating z-level below lowest z-ocean class
  integer, parameter, public :: GLC_ZOCNCLASS_ERR_TOO_HIGH = 3  ! err_code indicating z-level above highest z-ocean class

  ! String length for glc z-ocean classes represented as strings
  integer, parameter, public :: GLC_ZOCNCLASS_STRLEN = 2

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! number of elevation classes
  integer :: glc_nzoc  ! number of z-ocean classes

  ! z-level of each class.  Units are meters above sea level, so values should be <0
  ! indexing goes from shallowest to deepest levels
  real(r8), allocatable :: zocn_levels(:)
  ! upper and lower z-level limit for each class (m)
  ! first dimension: indexing goes from shallowest to deepest levels
  ! second dimension: index 1 is upper limit, index 2 is lower limit
  real(r8), allocatable :: zocn_bnds(:,:)


contains

  !-----------------------------------------------------------------------
  subroutine glc_zocnclass_init_default(my_glc_nzoc)
    !
    ! !DESCRIPTION:
    ! Initialize GLC z-ocean class data to default values, based on given glc_nzoc
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: my_glc_nzoc  ! number of GLC z-ocean classes
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_zocnclass_init'
    integer :: i
    !-----------------------------------------------------------------------

    glc_nzoc = my_glc_nzoc
    allocate(zocn_levels(glc_nzoc))

    select case (glc_nzoc)
    case(0)
       ! do nothing
    case(4)
       zocn_levels = [-250._r8,  -750._r8,  -1250._r8, -1750._r8]
    case(30)
       zocn_levels = [ -30._r8,   -90._r8,  -150._r8,  -210._r8,  -270._r8, &
                      -330._r8,  -390._r8,  -450._r8,  -510._r8,  -570._r8, &
                      -630._r8,  -690._r8,  -750._r8,  -810._r8,  -870._r8, &
                      -930._r8,  -990._r8, -1050._r8, -1110._r8, -1170._r8, &
                     -1230._r8, -1290._r8, -1350._r8, -1410._r8, -1470._r8, &
                     -1530._r8, -1590._r8, -1650._r8, -1710._r8, -1770._r8]
    case default
       write(logunit,*) subname,' ERROR: unknown glc_nzoc: ', glc_nzoc
       call shr_sys_abort(subname//' ERROR: unknown glc_nzoc')
    end select

    call glc_zocnclass_init_bnds()

  end subroutine glc_zocnclass_init_default

  !-----------------------------------------------------------------------
  subroutine glc_zocnclass_init_bnds()
    integer :: i

    allocate(zocn_bnds(2,glc_nzoc))
    zocn_bnds(:,:) = 0._r8
    if (glc_nzoc >= 2) then
       zocn_bnds(1,1) = 0._r8
       zocn_bnds(2,1) = 0.5_r8 * (zocn_levels(1) + zocn_levels(2))
       do i = 2, glc_nzoc - 1
          zocn_bnds(1,i) = 0.5_r8 * (zocn_levels(i-1) + zocn_levels(i))
          zocn_bnds(2,i) = 0.5_r8 * (zocn_levels(i) + zocn_levels(i+1))
       enddo
       zocn_bnds(1,glc_nzoc) = 0.5_r8 * (zocn_levels(glc_nzoc-1) + zocn_levels(glc_nzoc))
       zocn_bnds(2,glc_nzoc) = zocn_levels(glc_nzoc) + (zocn_levels(glc_nzoc) - zocn_bnds(1,glc_nzoc))
    endif
  end subroutine glc_zocnclass_init_bnds

  !-----------------------------------------------------------------------
  subroutine glc_zocnclass_init_override(my_glc_nzoc, my_zocn_levels)
    !
    ! !DESCRIPTION:
    ! Initialize GLC zocn class data to the given z-values
    !
    ! The input, my_zocn_levels, should have my_glc_nzoc elements.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: my_glc_nzoc     ! number of GLC z-ocean classes
    real(r8), intent(in) :: my_zocn_levels(:) ! z-ocean values (m)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_zocnlass_init_override'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(my_zocn_levels) == (/my_glc_nzoc/)), __FILE__, __LINE__)

    glc_nzoc = my_glc_nzoc
    allocate(zocn_levels(glc_nzoc))
    zocn_levels = my_zocn_levels
    allocate(zocn_bnds(2,glc_nzoc))

    call glc_zocnclass_init_bnds()

  end subroutine glc_zocnclass_init_override

  !-----------------------------------------------------------------------
  subroutine glc_zocnclass_clean()
    !
    ! !DESCRIPTION:
    ! Deallocate memory allocated in this module

    character(len=*), parameter :: subname = 'glc_zocnclass_clean'
    !-----------------------------------------------------------------------

    if (allocated(zocn_levels)) then
       deallocate(zocn_levels)
    end if
    if (allocated(zocn_bnds)) then
       deallocate(zocn_bnds)
    end if
    glc_nzoc = 0

  end subroutine glc_zocnclass_clean

  !-----------------------------------------------------------------------
  function glc_get_num_zocn_classes() result(num_zocn_classes)
    !
    ! !DESCRIPTION:
    ! Get the number of GLC z-ocean classes
    !
    ! !ARGUMENTS:
    integer :: num_zocn_classes  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_get_num_zocn_classes'
    !-----------------------------------------------------------------------

    num_zocn_classes = glc_nzoc

  end function glc_get_num_zocn_classes

  !-----------------------------------------------------------------------
  function glc_get_zlevels() result(zlevs)
    !
    ! !DESCRIPTION:
    ! Get all z-levels
    !
    ! This returns an array of size (glc_nzoc)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: zlevs(glc_nzoc)  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_get_zlevels'
    !-----------------------------------------------------------------------

    zlevs(:) = zocn_levels(:)

  end function glc_get_zlevels

  !-----------------------------------------------------------------------
  function glc_get_zocnclass_bounds() result(zocnclass_bounds)
    !
    ! !DESCRIPTION:
    ! Get the boundaries of all z-ocean classes.
    !
    ! This returns an array of size (glc_nzoc,2)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: zocnclass_bounds(2,glc_nzoc)  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_get_zocnclass_bounds'
    !-----------------------------------------------------------------------

    zocnclass_bounds(:,:) = zocn_bnds(:,:)

  end function glc_get_zocnclass_bounds

  !-----------------------------------------------------------------------
  function glc_zocnclass_as_string(zocn_class) result(zc_string)
    !
    ! !DESCRIPTION:
    ! Returns a string corresponding to a given elevation class.
    !
    ! This string can be used as a suffix for fields in MCT attribute vectors.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=GLC_ZOCNCLASS_STRLEN) :: zc_string  ! function result
    integer, intent(in) :: zocn_class
    !
    ! !LOCAL VARIABLES:
    character(len=16) :: format_string

    character(len=*), parameter :: subname = 'glc_zocnclass_as_string'
    !-----------------------------------------------------------------------

    ! e.g., for GLC_ZOCNCLASS_STRLEN = 2, format_string will be '(i2.2)'
    write(format_string,'(a,i0,a,i0,a)') '(i', GLC_ZOCNCLASS_STRLEN, '.', GLC_ZOCNCLASS_STRLEN, ')'

    write(zc_string,trim(format_string)) zocn_class
  end function glc_zocnclass_as_string

  !-----------------------------------------------------------------------
  function glc_all_zocnclass_strings(include_zero) result(zc_strings)
    !
    ! !DESCRIPTION:
    ! Returns an array of strings corresponding to all z-ocean classes from 1 to glc_nzoc
    !
    ! If include_zero is present and true, then includes z-ocean class 0 - so goes from
    ! 0 to glc_nzoc
    !
    ! These strings can be used as suffixes for fields in MCT attribute vectors.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=GLC_ZOCNCLASS_STRLEN), allocatable :: zc_strings(:)  ! function result
    logical, intent(in), optional :: include_zero   ! if present and true, include elevation class 0 (default is false)
    !
    ! !LOCAL VARIABLES:
    logical :: l_include_zero  ! local version of optional include_zero argument
    integer :: lower_bound
    integer :: i

    character(len=*), parameter :: subname = 'glc_all_zocnclass_strings'
    !-----------------------------------------------------------------------

    if (present(include_zero)) then
       l_include_zero = include_zero
    else
       l_include_zero = .false.
    end if

    if (l_include_zero) then
       lower_bound = 0
    else
       lower_bound = 1
    end if

    allocate(zc_strings(lower_bound:glc_nzoc))
    do i = lower_bound, glc_nzoc
       zc_strings(i) = glc_zocnclass_as_string(i)
    end do

  end function glc_all_zocnclass_strings


  !-----------------------------------------------------------------------
  function glc_zocn_errcode_to_string(err_code) result(err_string)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=256) :: err_string  ! function result
    integer, intent(in) :: err_code   ! error code (one of the GLC_ZOCNCLASS_ERR* values)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_errcode_to_string'
    !-----------------------------------------------------------------------

    select case (err_code)
    case (GLC_ZOCNCLASS_ERR_NONE)
       err_string = '(no error)'
    case (GLC_ZOCNCLASS_ERR_UNDEFINED)
       err_string = 'Z-ocean classes have not yet been defined'
    case (GLC_ZOCNCLASS_ERR_TOO_LOW)
       err_string = 'Z-level below the lower bound of the lowest z-ocean class'
    case (GLC_ZOCNCLASS_ERR_TOO_HIGH)
       err_string = 'Z-level above the upper bound of the highest z-ocean class'
    case default
       err_string = 'UNKNOWN ERROR'
    end select

  end function glc_zocn_errcode_to_string


end module glc_zocnclass_mod
