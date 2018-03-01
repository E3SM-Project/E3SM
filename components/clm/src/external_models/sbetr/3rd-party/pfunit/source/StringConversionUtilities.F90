!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: StringConversionUtilities
!
!> @brief
!! A collection of utilities used throughout the framework.
!!
!! @author
!! Tom Clune, NASA/GSFC 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 05 Sep 2014 - Added options for working with whitespace including
!               ignore, trim, or keep.  Note: trimAll trims both
!               sides, while trimTrailingWhitespace is more like
!               Fortran's trim. MLR
!
! 07 Nov 2013 - Added the prologue for the compliance with Doxygen.
!
!-------------------------------------------------------------------------------
! This module converts integers/real's to strings of a specific format
! for unit test diagnostics.  Basically just a wrapper for Fortran
! formatting, but functional programming provides a better style in many
! situations.
! 
! Further control of field width could be added at a later time.
!

module StringConversionUtilities_mod

  use Params_mod, only : r32, r64
  use Params_mod, only : i32, i64

   implicit none
   private
   
   public :: toString
   public :: appendWithSpace
   public :: MAXLEN_STRING
   public :: nullTerminate
   public :: unlessScalar
   public :: WhitespaceOptions, IGNORE_ALL, TRIM_ALL, KEEP_ALL, IGNORE_DIFFERENCES
   public :: whitespacep, trimAll, trimTrailingWhitespace

   integer, parameter :: MAXLEN_STRING = 80
!   integer, parameter :: MAXLEN_STRING = 80*5

   interface toString
      module Procedure toString_real64Scalar
      module Procedure toString_realScalar
      module Procedure toString_complex64Scalar
      module Procedure toString_complexScalar
      module Procedure toString_integerScalar_i32
      module Procedure toString_integer1D_i32
      module Procedure toString_integerScalar_i64
      module Procedure toString_integer1D_i64
   end interface

   character(len=*), parameter :: r32fmtStr = 'SP,G14.7'
   character(len=*), parameter :: r64fmtStr = 'SP,G14.7'
   character(len=*), parameter :: r32fmt1 = '('//r32fmtStr//')'
   character(len=*), parameter :: r64fmt1 = '('//r64fmtStr//')'

   character(len=*), parameter :: c32fmt1 = '("z=(",'//r32fmt1//',",",'//r32fmt1//',")")'
   character(len=*), parameter :: c64fmt1 = '("z=(",'//r64fmt1//',",",'//r64fmt1//',")")'

!   enum, bind(c) :: WhitespaceOptions
   type WhitespaceOptions
      integer value
   end type WhitespaceOptions
   enum, bind(c)
      enumerator :: IGNORE_ALL_, TRIM_ALL_, KEEP_ALL_, IGNORE_DIFFERENCES_
   end enum
   type (WhitespaceOptions), parameter :: &
        & IGNORE_ALL=WhitespaceOptions(IGNORE_ALL_), &
        & TRIM_ALL  =WhitespaceOptions(TRIM_ALL_), &
        & KEEP_ALL  =WhitespaceOptions(KEEP_ALL_), &
        & IGNORE_DIFFERENCES =WhitespaceOptions(IGNORE_DIFFERENCES_)

contains

   character(len=MAXLEN_STRING) function toString_complex64Scalar(value) result(buffer)
      complex(kind=r64), intent(in) :: value

!      write(buffer,'(2(SP,G14.7))') value
      write(buffer,c64fmt1) value
      buffer = adjustL(buffer)

    end function toString_complex64Scalar

   character(len=MAXLEN_STRING) function toString_complexScalar(value) result(buffer)
      complex, intent(in) :: value

!      write(buffer,'(2(SP,G14.7))') value
      write(buffer,c32fmt1) value
      buffer = adjustL(buffer)

   end function toString_complexScalar

   character(len=MAXLEN_STRING) function toString_real64Scalar(value) result(buffer)
      real(kind=r64), intent(in) :: value

      write(buffer,'(SP,G14.7)') value
!      write(buffer,r64fmt1) value
      buffer = adjustL(buffer)

    end function toString_real64Scalar

   character(len=MAXLEN_STRING) function toString_realScalar(value) result(buffer)
      real(kind=r32), intent(in) :: value

      write(buffer,'(SP,G14.7)') value
!      print *,'r32fmt1: ',r32fmt1
!      print *,'       : ','(SP,G14.7)'
!      print *,'=?     : ','(SP,G14.7)'.EQ.r32fmt1
!      write(buffer,r32fmt1)
      buffer = adjustL(buffer)

   end function toString_realScalar

!-   character(len=MAXLEN_STRING) function toString_integerScalar(value) result(buffer)
!-      integer, intent(in) :: value
!-      character(len=20) :: fmt
!-
!-      fmt = '(I0)'
!-      write(buffer,trim(fmt)) value
!-      buffer = adjustL(buffer)
!-
!-   end function toString_integerScalar
!-
!-   function toString_integer1D(arrayShape) result(string)
!-      integer, intent(in) :: arrayShape(:)
!-      character(len=MAXLEN_STRING) :: string
!-
!-!      integer :: i
!-      
!-      select case (size(arrayShape)) ! rank
!-      case (0) ! scalar
!-         string = '0'
!-      case (1)
!-         write(string,'(i0)') arrayShape(1)
!-      case (2:)
!-         write(string,'(i0,14(",",i0:))') arrayShape(1:)
!-      end select
!-
!-      string = '[' // trim(string) // ']'
!-   end function toString_integer1D

   character(len=MAXLEN_STRING) function toString_integerScalar_i32(value) result(buffer)
      integer(kind=i32), intent(in) :: value
      character(len=20) :: fmt

      fmt = '(I0)'
      write(buffer,trim(fmt)) value
      buffer = adjustL(buffer)

    end function toString_integerScalar_i32

   function toString_integer1D_i32(arrayShape) result(string)
      integer(kind=i32), intent(in) :: arrayShape(:)
      character(len=MAXLEN_STRING) :: string

!      integer :: i
      
      select case (size(arrayShape)) ! rank
      case (0) ! scalar
         string = '0'
      case (1)
         write(string,'(i0)') arrayShape(1)
      case (2:)
         write(string,'(i0,14(",",i0:))') arrayShape(1:)
      end select

      string = '[' // trim(string) // ']'
    end function toString_integer1D_i32

   character(len=MAXLEN_STRING) function toString_integerScalar_i64(value) result(buffer)
      integer(kind=i64), intent(in) :: value
      character(len=20) :: fmt

      fmt = '(I0)'
      write(buffer,trim(fmt)) value
      buffer = adjustL(buffer)

    end function toString_integerScalar_i64

   function toString_integer1D_i64(arrayShape) result(string)
      integer(kind=i64), intent(in) :: arrayShape(:)
      character(len=MAXLEN_STRING) :: string

!      integer :: i
      
      select case (size(arrayShape)) ! rank
      case (0) ! scalar
         string = '0'
      case (1)
         write(string,'(i0)') arrayShape(1)
      case (2:)
         write(string,'(i0,14(",",i0:))') arrayShape(1:)
      end select

      string = '[' // trim(string) // ']'
    end function toString_integer1D_i64

   
   ! Joins two strings with a space separator unless first string is
   ! empty.
   function appendWithSpace(a, b) result(ab)
      character(len=*), intent(in) :: a
      character(len=*), intent(in) :: b
      character(len=len_trim(a)+1+len_trim(b)) :: ab

      if (len_trim(a) > 0) then
         ab = trim(a) // ' ' // trim(b)
      else
         ab = trim(b)
      end if

   end function appendWithSpace

   function nullTerminate(string) result(nullTerminatedString)
      use iso_c_binding
      character(len=*), intent(in) :: string
      character(len=:), allocatable :: nullTerminatedString

      nullTerminatedString = trim(string) // C_NULL_CHAR

   end function nullTerminate

   function unlessScalar(vShape,string) result(retString)
     integer, intent(in), dimension(:) :: vShape
     character(len=*), intent(in) :: string
     character(len=:), allocatable :: retString
     retString=""
     if(size(vShape).ne.0)then
        retString=string
     end if
   end function unlessScalar

   logical function whitespacep(c)
     character, intent(in) :: c
     integer, parameter :: iachar_spc = 32, iachar_tab = 9
     whitespacep = &
          & iachar(c) .eq. iachar_spc .or. &
          & iachar(c) .eq. iachar_tab
   end function whitespacep

   function trimAll(s) result(trimmed)
     character(len=*), intent(in) :: s
     character(len=:), allocatable :: trimmed
     integer :: i,lenS,leadingWhite,trailingWhite,lenTrimmed

     lenS = len(s)

     leadingWhite = 0
     do i = 1,lenS
        if (whitespacep(s(i:i))) then
           leadingWhite = leadingWhite + 1
        else
           exit
        end if
     end do

     trailingWhite = 0
     do i = lenS,leadingWhite+1,-1
        if (whitespacep(s(i:i))) then
           trailingWhite = trailingWhite + 1
        else
           exit
        end if
     end do
     lenTrimmed = lenS-leadingWhite-trailingWhite
     
     allocate(character(lenTrimmed) :: trimmed)
     do i = 1,lenTrimmed
        trimmed(i:i) = s(i+leadingWhite:i+leadingWhite)
     end do
   end function trimAll

   function trimTrailingWhitespace(s) result(trimmed)
     character(len=*), intent(in) :: s
     character(len=:), allocatable :: trimmed
     integer :: i,lenS
     integer :: trailingWhite,lenTrimmed
     integer :: leadingWhite

     lenS = len(s)

     leadingWhite = 0
     do i = 1,lenS
        if (whitespacep(s(i:i))) then
           leadingWhite = leadingWhite + 1
        else
           exit
        end if
     end do

     trailingWhite = 0
     do i = lenS,leadingWhite+1,-1
        if (whitespacep(s(i:i))) then
           trailingWhite = trailingWhite + 1
        else
           exit
        end if
     end do

     lenTrimmed = lenS-trailingWhite
     allocate(character(lenTrimmed) :: trimmed)
     do i = 1,lenTrimmed
        trimmed(i:i) = s(i:i)
     end do

   end function trimTrailingWhitespace



     



end module StringConversionUtilities_mod
