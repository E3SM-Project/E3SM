!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: UnixPipeInterfaces
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
!-----------
! This module provides Bind(C) interfaces for
! standard unix pipe utilities:  popen(), fgets(), and pclose().
!
! These utilities work at a very low level, and should generally not
! be used directly from a user application.  E.g. strings must be null
! terminated.
! UnixProcess.F90 provides a (customized) higher-level interface that
! should be safer for routine use.
!-------------

module UnixPipeInterfaces_mod
   use, intrinsic :: ISO_C_BINDING
   private

   public :: popen
   public :: fgets
   public :: pclose
   public :: getline
   public :: getdelim
   public :: free

   ! error codes
   public :: CLOSE_FAILED

   integer(C_INT), parameter :: CLOSE_FAILED = -1
 
   interface

      function popen(command, mode) result(file) bind(C, name='popen')
         use, intrinsic :: iso_c_binding
         type (C_PTR) :: file
         character(kind=C_CHAR), dimension(*), intent(in) :: command
         character(kind=C_CHAR), dimension(*), intent(in) :: mode
      end function popen

      function fgets(str, size, stream) bind(C, name='fgets')
         use, intrinsic :: iso_c_binding
         type (C_PTR) :: fgets
         character(kind=C_CHAR), dimension(*), intent(inout) :: str
         integer(kind=C_INT), value, intent(in) :: size
         type (C_PTR), value :: stream
      end function fgets

      function pclose(stream) bind(C, name='pclose')
         use, intrinsic :: iso_c_binding
         integer(C_INT) :: pclose
         type (c_ptr), value :: stream
      end function pclose

      function getline(linep, linecapp, stream) bind(C, name='getline')
         use, intrinsic :: iso_c_binding
         integer (kind=C_SIZE_T) :: getline
         type (C_PTR) :: linep
         integer (kind=C_SIZE_T) :: linecapp
         type (C_PTR), value :: stream
      end function getline

      function getdelim(linep, linecapp, delimeter, stream) bind(C, name='getdelim')
         use, intrinsic :: iso_c_binding
         integer (kind=C_SIZE_T) :: getdelim
         type (C_PTR) :: linep
         integer (kind=C_SIZE_T) :: linecapp
         integer(kind=C_INT), value :: delimeter
         type (C_PTR), value :: stream
      end function getdelim

      subroutine free(ptr) bind(C, name='free')
         use, intrinsic :: iso_c_binding
         type (C_PTR), value :: ptr
      end subroutine free

   end interface

end module UnixPipeInterfaces_mod
