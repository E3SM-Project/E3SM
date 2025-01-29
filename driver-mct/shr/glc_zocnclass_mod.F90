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
  public :: glc_get_zocn_class            ! get the z-ocean class index for a given z-level
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

  ! upper z limit of each class (m)
  ! indexing starts at 0, with zocnmax(0) giving the lower elevation limit of z-ocean class 1
  ! indexing goes from ocean surface to deeper levels
  real(r8), allocatable :: zocnmax(:)


contains

  !-----------------------------------------------------------------------
  subroutine glc_zocnclass_init_default(my_glc_nzoc)
    !
    ! !DESCRIPTION:
    ! Initialize GLC -ocean class data to default boundaries, based on given glc_nzoc
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: my_glc_nzoc  ! number of GLC z-ocean classes
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_zocnclass_init'
    !-----------------------------------------------------------------------

    glc_nzoc = my_glc_nzoc
    allocate(zocnmax(0:glc_nzoc))

    select case (glc_nzoc)
    case(0)
       ! do nothing
    case(4)
       zocnmax = [0._r8,  -500._r8,  -1000._r8, -1500._r8, -2000._r8]
    case(30)
       zocnmax = [   0._r8,   -60._r8,  -120._r8,  -180._r8,  -240._r8, &
                  -300._r8,  -360._r8,  -420._r8,  -480._r8,  -540._r8, &
                  -600._r8,  -660._r8,  -720._r8,  -780._r8,  -840._r8, &
                  -900._r8,  -960._r8, -1020._r8, -1080._r8, -1140._r8, &
                 -1200._r8, -1260._r8, -1320._r8, -1380._r8, -1440._r8, &
                 -1500._r8, -1560._r8, -1620._r8, -1680._r8, -1740._r8]
    case default
       write(logunit,*) subname,' ERROR: unknown glc_nzoc: ', glc_nzoc
       call shr_sys_abort(subname//' ERROR: unknown glc_nzoc')
    end select

  end subroutine glc_zocnclass_init_default

  !-----------------------------------------------------------------------
  subroutine glc_zocnclass_init_override(my_glc_nzoc, my_zocnmax)
    !
    ! !DESCRIPTION:
    ! Initialize GLC zocn class data to the given z-ocean class boundaries.
    !
    ! The input, my_zocnmax, should have (my_glc_nzoc + 1) elements.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: my_glc_nzoc     ! number of GLC z-ocean classes
    real(r8), intent(in) :: my_zocnmax(0:) ! z-ocean class boundaries (m)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_zocnlass_init_override'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(my_zocnmax) == (/my_glc_nzoc/)), __FILE__, __LINE__)

    glc_nzoc = my_glc_nzoc
    allocate(zocnmax(0:glc_nzoc))
    zocnmax = my_zocnmax

  end subroutine glc_zocnclass_init_override

  !-----------------------------------------------------------------------
  subroutine glc_zocnclass_clean()
    !
    ! !DESCRIPTION:
    ! Deallocate memory allocated in this module

    character(len=*), parameter :: subname = 'glc_zocnclass_clean'
    !-----------------------------------------------------------------------

    if (allocated(zocnmax)) then
       deallocate(zocnmax)
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
  subroutine glc_get_zocn_class(zlev, zocn_class, err_code)
    !
    ! !DESCRIPTION:
    ! Get the zocn class index associated with a given ocean z-level (depth).
    !
    ! The returned zocn_class will be between 1 and num_zocn_classes, if this
    ! z-level is contained in a z-ocean class. In this case, err_code will
    ! be GLC_ZOCNCLASS_ERR_NONE (no error).
    !
    ! If there are no z-ocean classes defined, the returned value will be 0, and
    ! err_code will be GLC_ZOCNCLASS_ERR_UNDEFINED
    !
    ! If this z-level is below the lowest zocean class, the returned value
    ! will be 1, and err_code will be GLC_ZOCNCLASS_ERR_TOO_LOW.
    !
    ! If this z-level is above the highest z-ocean class, the returned value
    ! will be (num_zocn_classes), and err_code will be GLC_ZOCNCLASS_ERR_TOO_HIGH.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: zlev ! z-level in ocean (depth) (m)
    integer, intent(out) :: zocn_class  ! z-ocean class index
    integer, intent(out) :: err_code ! error code (see above for possible codes)
    !
    ! !LOCAL VARIABLES:
    integer :: zc  ! temporary z-ocean class

    character(len=*), parameter :: subname = 'glc_get_zocn_class'
    !-----------------------------------------------------------------------

    if (glc_nzoc < 1) then
       zocn_class = 0
       err_code = GLC_ZOCNCLASS_ERR_UNDEFINED
    else if (zlev < zocnmax(0)) then
       zocn_class = 1
       err_code = GLC_ZOCNCLASS_ERR_TOO_LOW
    else if (zlev >= zocnmax(glc_nzoc)) then
       zocn_class = glc_nzoc
       err_code = GLC_ZOCNCLASS_ERR_TOO_HIGH
    else
       err_code = GLC_ZOCNCLASS_ERR_NONE
       zocn_class = 0
       do zc = 1, glc_nzoc
          if (zlev >= zocnmax(zc - 1) .and. zlev < zocnmax(zc)) then
             zocn_class = zc
             exit
          end if
       end do

       SHR_ASSERT(zocn_class > 0, subname//' z-ocean class was not assigned')
    end if

  end subroutine glc_get_zocn_class

  !-----------------------------------------------------------------------
  function glc_get_zocnclass_bounds() result(zocnclass_bounds)
    !
    ! !DESCRIPTION:
    ! Get the boundaries of all z-ocean classes.
    !
    ! This returns an array of size glc_nzoc+1, since it contains both the lower and upper
    ! bounds of each z-ocean class.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: zocnclass_bounds(0:glc_nzoc)  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_get_zocnclass_bounds'
    !-----------------------------------------------------------------------

    zocnclass_bounds(:) = zocnmax(:)

  end function glc_get_zocnclass_bounds

  !-----------------------------------------------------------------------
  function glc_zocnclass_as_string(zocn_class) result(zc_string)
    !
    ! !DESCRIPTION:
    ! Returns a string corresponding to a given elevation class.
    !
    ! This string can be used as a suffix for fields in MCT attribute vectors.
    !
    ! ! NOTE(wjs, 2015-01-19) This function doesn't fully belong in this module, since it
    ! doesn't refer to the data stored in this module. However, I can't think of a more
    ! appropriate place for it.
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
