!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: SourceLocation
!
!> @brief
!! <BriefDescription>
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
! 07 Nov 2013 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------
! This module just provides a data type - not a class.
! Meant to be shared for easy access.

module SourceLocation_mod
   implicit none
   private

   public :: SourceLocation
   public :: UNKNOWN_SOURCE_LOCATION
   public :: UNKNOWN_FILE_NAME
   public :: UNKNOWN_LINE_NUMBER

   integer, parameter :: MAXLEN_FILE_NAME = 100
   character(len=MAXLEN_FILE_NAME), parameter :: UNKNOWN_FILE_NAME= '<unknown file>'
   integer, parameter :: UNKNOWN_LINE_NUMBER = -1

   type :: SourceLocation
      character(len=MAXLEN_FILE_NAME) :: fileName = UNKNOWN_FILE_NAME
      integer :: lineNumber = UNKNOWN_LINE_NUMBER
   contains
      procedure :: toString
   end type SourceLocation

   type (SourceLocation), parameter :: UNKNOWN_SOURCE_LOCATION = &
        & SourceLocation()

contains

   function toString(this) result(string)
      class (SourceLocation), intent(inout) :: this
      character(len=80) :: string
      
      if (this%fileName == UNKNOWN_FILE_NAME) then
         if (this%lineNumber == UNKNOWN_LINE_NUMBER) then
            string = '<unknown location>'
         else
            write(string,'(a,":",i0)') trim(UNKNOWN_FILE_NAME), this%lineNumber
         end if
      else
         if (this%lineNumber == UNKNOWN_LINE_NUMBER) then
            string = trim(this%fileName)
         else
            write(string,'(a,":",i0)') trim(this%fileName), this%lineNumber
         end if
      end if

      string = '[' // trim(string) // ']'

   end function toString

end module SourceLocation_mod
