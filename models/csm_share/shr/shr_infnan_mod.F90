
! Flag representing compiler support of Fortran 2003's
! ieee_arithmetic intrinsic module.
#if defined CPRIBM || defined CPRPGI || defined CPRINTEL ||  defined CPRCRAY || defined CPRNAG
#define HAVE_IEEE_ARITHMETIC
#endif

module shr_infnan_mod
!---------------------------------------------------------------------
! Module to test for IEEE Inf and NaN values, which also provides a
! method of setting +/-Inf and signaling or quiet NaN.
!
! All functions are elemental, and thus work on arrays.
!---------------------------------------------------------------------
! To test for these values, just call the corresponding function, e.g:
!
!   var_is_nan = shr_infnan_isnan(x)
!
! You can also use it on arrays:
!
!   array_contains_nan = any(shr_infnan_isnan(my_array))
!
!---------------------------------------------------------------------
! To generate these values, assign one of the provided derived-type
! variables to a real:
!
!   use shr_infnan_mod, only: nan => shr_infnan_nan, &
!                             inf => shr_infnan_inf, &
!                             assignment(=)
!   real(r4) :: my_nan
!   real(r8) :: my_inf_array(2,2)
!   my_nan = nan
!   my_inf_array = inf
!
! Keep in mind that "shr_infnan_nan" and "shr_infnan_inf" cannot be
! passed to functions that expect real arguments. To pass a real
! NaN, you will have to use shr_infnan_nan to set a local real of
! the correct kind.
!---------------------------------------------------------------------

use shr_kind_mod, only: &
     r4 => SHR_KIND_R4, &
     r8 => SHR_KIND_R8

#ifdef HAVE_IEEE_ARITHMETIC

! If we have IEEE_ARITHMETIC, NaN tests are provided
! for us.
use, intrinsic :: ieee_arithmetic, only: &
     shr_infnan_isnan => ieee_is_nan

#else

! If we don't use Fortran 2003 IEEE module,
! use the old code to do all NaN/Inf testing.

use shr_test_infnan_mod, only: &
     shr_infnan_isnan => shr_test_infnan_isnan, &
     shr_infnan_isinf => shr_test_infnan_isinf, &
     shr_infnan_isposinf => shr_test_infnan_isposinf, &
     shr_infnan_isneginf => shr_test_infnan_isneginf

#endif

implicit none
private
save

! Test functions for NaN/Inf values.
public :: shr_infnan_isnan
public :: shr_infnan_isinf
public :: shr_infnan_isposinf
public :: shr_infnan_isneginf

#ifdef HAVE_IEEE_ARITHMETIC

! If we do have the IEEE intrinsics, define the
! inf test functions here.
interface shr_infnan_isinf
   module procedure shr_infnan_isinf_r4
   module procedure shr_infnan_isinf_r8
end interface

interface shr_infnan_isposinf
   module procedure shr_infnan_isposinf_r4
   module procedure shr_infnan_isposinf_r8
end interface

interface shr_infnan_isneginf
   module procedure shr_infnan_isneginf_r4
   module procedure shr_infnan_isneginf_r8
end interface

#endif

! Derived types for generation of NaN/Inf
! Even though there's no reason to "use"
! the types directly, Lahey breaks if you make
! shr_infnan_nan public, but shr_infnan_nan_type
! private.
public :: shr_infnan_nan_type
public :: shr_infnan_inf_type
public :: assignment(=)

! Type representing Not A Number.
type :: shr_infnan_nan_type
   logical :: quiet = .false.
end type shr_infnan_nan_type

! Type representing +/-Infinity.
type :: shr_infnan_inf_type
   logical :: positive = .true.
end type shr_infnan_inf_type

! Allow assigning reals to NaN or Inf.
interface assignment(=)
   module procedure set_r4_nan
   module procedure set_r4_inf
   module procedure set_r8_nan
   module procedure set_r8_inf
end interface

! Initialize objects of NaN/Inf type for other modules to use.

! Default NaN is signaling, but also provide snan and qnan to choose
! explicitly.
type(shr_infnan_nan_type), public, parameter :: shr_infnan_nan = &
     shr_infnan_nan_type(.false.)
type(shr_infnan_nan_type), public, parameter :: shr_infnan_snan = &
     shr_infnan_nan_type(.false.)
type(shr_infnan_nan_type), public, parameter :: shr_infnan_qnan = &
     shr_infnan_nan_type(.true.)

! Default Inf is positive, but provide posinf to go with neginf.
type(shr_infnan_inf_type), public, parameter :: shr_infnan_inf = &
     shr_infnan_inf_type(.true.)
type(shr_infnan_inf_type), public, parameter :: shr_infnan_posinf = &
     shr_infnan_inf_type(.true.)
type(shr_infnan_inf_type), public, parameter :: shr_infnan_neginf = &
     shr_infnan_inf_type(.false.)

contains

#ifdef HAVE_IEEE_ARITHMETIC

!---------------------------------------------------------------------
! TEST FUNCTIONS
!---------------------------------------------------------------------
! The "isinf" functions simply call "isposinf" and "isneginf".
!---------------------------------------------------------------------

elemental function shr_infnan_isinf_r4(x) result(isinf)
  real(r4), intent(in) :: x
  logical :: isinf

  isinf = shr_infnan_isposinf(x) .or. shr_infnan_isneginf(x)

end function shr_infnan_isinf_r4

elemental function shr_infnan_isinf_r8(x) result(isinf)
  real(r8), intent(in) :: x
  logical :: isinf

  isinf = shr_infnan_isposinf(x) .or. shr_infnan_isneginf(x)

end function shr_infnan_isinf_r8

!---------------------------------------------------------------------
! The "isposinf" and "isneginf" functions get the IEEE class of a
! real, and test to see if the class is equal to ieee_positive_inf
! or ieee_negative_inf.
!---------------------------------------------------------------------

elemental function shr_infnan_isposinf_r4(x) result(isposinf)
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_class, &
       ieee_positive_inf, &
       operator(==)
       
  real(r4), intent(in) :: x
  logical :: isposinf

  isposinf = (ieee_positive_inf == ieee_class(x))

end function shr_infnan_isposinf_r4

elemental function shr_infnan_isposinf_r8(x) result(isposinf)
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_class, &
       ieee_positive_inf, &
       operator(==)
       
  real(r8), intent(in) :: x
  logical :: isposinf

  isposinf = (ieee_positive_inf == ieee_class(x))

end function shr_infnan_isposinf_r8

elemental function shr_infnan_isneginf_r4(x) result(isneginf)
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_class, &
       ieee_negative_inf, &
       operator(==)
       
  real(r4), intent(in) :: x
  logical :: isneginf

  isneginf = (ieee_negative_inf == ieee_class(x))

end function shr_infnan_isneginf_r4

elemental function shr_infnan_isneginf_r8(x) result(isneginf)
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_class, &
       ieee_negative_inf, &
       operator(==)
       
  real(r8), intent(in) :: x
  logical :: isneginf

  isneginf = (ieee_negative_inf == ieee_class(x))

end function shr_infnan_isneginf_r8

#endif

!---------------------------------------------------------------------
! GENERATION FUNCTIONS
!---------------------------------------------------------------------
! Two approaches for generation of NaN and Inf values:
!   1. With Fortran 2003, use the ieee_value intrinsic to get a value
!      from the corresponding class. These are:
!       - ieee_signaling_nan
!       - ieee_quiet_nan
!       - ieee_positive_inf
!       - ieee_negative_inf
!   2. Without Fortran 2003, set the IEEE bit patterns directly.
!      Use BOZ literals to get an integer with the correct bit
!      pattern, then use "transfer" to transfer those bits into a
!      real.
!---------------------------------------------------------------------

elemental subroutine set_r4_nan(output, nan)
#ifdef HAVE_IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_signaling_nan, &
       ieee_quiet_nan, &
       ieee_value
#else
  integer, parameter :: ssnan_pat = Z'7FA00000'
  integer, parameter :: sqnan_pat = Z'7FC00000'
#endif
  real(r4), intent(out) :: output
  type(shr_infnan_nan_type), intent(in) :: nan

#ifdef HAVE_IEEE_ARITHMETIC
  if (nan%quiet) then
     output = ieee_value(output,ieee_quiet_nan)
  else
     output = ieee_value(output,ieee_signaling_nan)
  end if
#else
  if (nan%quiet) then
     output = transfer(sqnan_pat, output)
  else
     output = transfer(ssnan_pat, output)
  end if
#endif

end subroutine set_r4_nan

elemental subroutine set_r4_inf(output, inf)
#ifdef HAVE_IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_positive_inf, &
       ieee_negative_inf, &
       ieee_value
#else
  integer, parameter :: sposinf_pat = Z'7F800000'
  integer, parameter :: sneginf_pat = Z'FF800000'
#endif
  real(r4), intent(out) :: output
  type(shr_infnan_inf_type), intent(in) :: inf

#ifdef HAVE_IEEE_ARITHMETIC
  if (inf%positive) then
     output = ieee_value(output,ieee_positive_inf)
  else
     output = ieee_value(output,ieee_negative_inf)
  end if
#else
  if (inf%positive) then
     output = transfer(sposinf_pat, output)
  else
     output = transfer(sneginf_pat, output)
  end if
#endif

end subroutine set_r4_inf

elemental subroutine set_r8_nan(output, nan)
#ifdef HAVE_IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_signaling_nan, &
       ieee_quiet_nan, &
       ieee_value
#else
  use shr_kind_mod, only: &
       i8 => shr_kind_i8
  integer(i8), parameter :: dsnan_pat = Z'7FF4000000000000'
  integer(i8), parameter :: dqnan_pat = Z'7FF8000000000000'
#endif
  real(r8), intent(out) :: output
  type(shr_infnan_nan_type), intent(in) :: nan

#ifdef HAVE_IEEE_ARITHMETIC
  if (nan%quiet) then
     output = ieee_value(output,ieee_quiet_nan)
  else
     output = ieee_value(output,ieee_signaling_nan)
  end if
#else
  if (nan%quiet) then
     output = transfer(dqnan_pat, output)
  else
     output = transfer(dsnan_pat, output)
  end if
#endif

end subroutine set_r8_nan

elemental subroutine set_r8_inf(output, inf)
#ifdef HAVE_IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_positive_inf, &
       ieee_negative_inf, &
       ieee_value
#else
  use shr_kind_mod, only: &
       i8 => shr_kind_i8
  integer(i8), parameter :: dposinf_pat = Z'7FF0000000000000'
  integer(i8), parameter :: dneginf_pat = Z'FFF0000000000000'
#endif
  real(r8), intent(out) :: output
  type(shr_infnan_inf_type), intent(in) :: inf

#ifdef HAVE_IEEE_ARITHMETIC
  if (inf%positive) then
     output = ieee_value(output,ieee_positive_inf)
  else
     output = ieee_value(output,ieee_negative_inf)
  end if
#else
  if (inf%positive) then
     output = transfer(dposinf_pat, output)
  else
     output = transfer(dneginf_pat, output)
  end if
#endif

end subroutine set_r8_inf

end module shr_infnan_mod
