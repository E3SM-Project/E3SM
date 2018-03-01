!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: ParallelContext
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
! Default implementation is for a serial process
! Subclass for MPI operations
module ParallelContext_mod
   implicit none
   private

   public :: ParallelContext

   type, abstract :: ParallelContext
   contains
      procedure :: isActive
      procedure :: isRootProcess
      procedure(getNumProcesses), deferred :: getNumProcesses
      procedure(processRank), deferred :: processRank
      procedure(sum), deferred :: sum
      generic :: gather => gatherString
      generic :: gather => gatherInteger
      generic :: gather => gatherLogical
      procedure(gatherString), deferred :: gatherString
      procedure(gatherInteger), deferred :: gatherInteger
      procedure(gatherLogical), deferred :: gatherLogical
      procedure :: labelProcess
      procedure :: barrier
      procedure(allReduceLogical), deferred :: allReduce
   end type ParallelContext

   abstract interface
      integer function getNumProcesses(this)
         import ParallelContext
         class(ParallelContext), intent(in) :: this
      end function getNumProcesses

      integer function processRank(this)
         import ParallelContext
         class(ParallelContext), intent(in) :: this
      end function processRank

      integer function sum(this, value)
         import ParallelContext
         class (ParallelContext), intent(in) :: this
         integer, intent(in) :: value
       end function sum
       
      subroutine gatherString(this, values, list)
         import ParallelContext
         class (ParallelContext), intent(in) :: this
         character(len=*), intent(in) :: values(:)
         character(len=*), intent(out) :: list(:)
      end subroutine gatherString

      subroutine gatherInteger(this, values, list)
         import ParallelContext
         class (ParallelContext), intent(in) :: this
         integer, intent(in) :: values(:)
         integer, intent(out) :: list(:)
      end subroutine gatherInteger

      subroutine gatherLogical(this, values, list)
         import ParallelContext
         class (ParallelContext), intent(in) :: this
         logical, intent(in) :: values(:)
         logical, intent(out) :: list(:)
      end subroutine gatherLogical

      logical function allReduceLogical(this, q) result(anyQ)
         import ParallelContext
         class (ParallelContext), intent(in) :: this
         logical, intent(in) :: q
      end function allReduceLogical

   end interface

contains

   logical function isActive(this)
      class (ParallelContext),  intent(in) :: this
      isActive = .true.
   end function isActive

   subroutine barrier(this)
      class (ParallelContext), intent(in) :: this
   end subroutine barrier

   logical function isRootProcess(this)
      class (ParallelContext), intent(in) :: this
      isRootProcess = .true.
   end function isRootProcess

   subroutine labelProcess(this, message)
      class (ParallelContext), intent(in) :: this
      character(len=*), intent(inout) :: message
   end subroutine labelProcess

end module ParallelContext_mod
