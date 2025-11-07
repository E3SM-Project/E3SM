module shr_strconvert_mod

! This module defines toString, a generic function for creating character type
! representations of data, as implemented for the most commonly used intrinsic
! types:
!
! - 4 and 8 byte integer
! - 4 and 8 byte real
! - logical
!
! No toString implementation is provided for character input, but this may be
! added if some use case arises.
!
! Currently, only scalar inputs are supported. The return type of this function
! is character with deferred (allocatable) length.
!
! The functions for integers and reals allow an optional format_string argument,
! which can be used to control the padding and precision of output as with any
! write statement. However, the implementations internally must use a
! preallocated buffer, so a format_string that significantly increases the size
! of the output may cause a run-time error or undefined behavior in the program.
!
! Other modules may want to provide extensions of toString for their own derived
! types. In this case there are two guidelines to observe:
!
! - It is preferable to have only one mandatory argument, which is the object to
!   produce a string from. There may be other formatting options, but the
!   implementation should do something sensible without these.
!
! - Since the main purpose of toString is to provide a human-readable
!   representation of a type, especially for documentation or debugging
!   purposes, refrain from printing large array components in their entirety
!   (instead consider printing only the shape, or statistics such as
!   min/mean/max for arrays of numbers).

use shr_kind_mod, only: &
     i4 => shr_kind_i4, &
     i8 => shr_kind_i8, &
     r4 => shr_kind_r4, &
     r8 => shr_kind_r8, &
     cs => shr_kind_cs

use shr_infnan_mod, only: &
     isnan => shr_infnan_isnan

implicit none
private

! Human-readable representation of data.
public :: toString

interface toString
   module procedure i4ToString
   module procedure i8ToString
   module procedure r4ToString
   module procedure r8ToString
   module procedure logicalToString
end interface toString

contains

pure character(len=cs) function i4ToString(input, format_string)
  integer(i4), intent(in) :: input
  character(len=*), intent(in), optional :: format_string

  if (present(format_string)) then
     write(i4ToString, format_string) input
  else
     ! For most compilers, these two statements are equivalent to a format of
     ! '(I0)', but that's not technically in the standard.
     write(i4ToString, '(I11)') input
     i4ToString = adjustl(i4ToString)
  end if

  i4ToString = trim(i4ToString)

end function i4ToString

pure character(len=cs) function i8ToString(input, format_string)
  integer(i8), intent(in) :: input
  character(len=*), intent(in), optional :: format_string

  if (present(format_string)) then
     write(i8ToString, format_string) input
  else
     ! For most compilers, these two statements are equivalent to a format of
     ! '(I0)', but that's not technically in the standard.
     write(i8ToString, '(I20)') input
     i8ToString = adjustl(i8ToString)
  end if

  i8ToString = trim(i8ToString)

end function i8ToString

pure character(len=cs) function r4ToString(input, format_string)
  real(r4), intent(in) :: input
  character(len=*), intent(in), optional :: format_string

  if (present(format_string)) then
     write(r4ToString, format_string) input
  else
     write(r4ToString, '(ES15.8 E2)') input
     r4ToString = adjustl(r4ToString)
     ! Deal with the fact that the "+" sign is optional by simply adding it if
     ! it is not present, so that the default format is standardized across
     ! compilers.
     ! Assumes that compilers do not treat the sign bit on NaN values specially.
     if (.not. isnan(input) .and. all(r4ToString(1:1) /= ["-", "+"])) then
        r4ToString = "+" // trim(r4ToString)
     end if
  end if

  r4ToString = trim(r4ToString)

end function r4ToString

pure character(len=cs) function r8ToString(input, format_string)
  real(r8), intent(in) :: input
  character(len=*), intent(in), optional :: format_string

  if (present(format_string)) then
     write(r8ToString, format_string) input
  else
     write(r8ToString, '(ES24.16 E3)') input
     r8ToString = adjustl(r8ToString)
     ! Deal with the fact that the "+" sign is optional by simply adding it if
     ! it is not present, so that the default format is standardized across
     ! compilers.
     ! Assumes that compilers do not treat the sign bit on NaN values specially.
     if (.not. isnan(input) .and. all(r8ToString(1:1) /= ["-", "+"])) then
        r8ToString = "+" // trim(r8ToString)
     end if
  end if

  r8ToString = trim(r8ToString)

end function r8ToString

pure character(len=cs) function logicalToString(input)
  logical, intent(in) :: input
  if (input) then
     logicalToString = "T"
  else
     logicalToString = "F"
  end if

end function logicalToString

end module shr_strconvert_mod
