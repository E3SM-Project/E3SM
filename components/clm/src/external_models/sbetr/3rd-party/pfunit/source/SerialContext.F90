!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: SerialContext
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
module SerialContext_mod
   use ParallelContext_mod
   implicit none
   private

   public :: SerialContext
   public :: newSerialContext
   public :: THE_SERIAL_CONTEXT

   type, extends(ParallelContext) :: SerialContext
   contains
      procedure :: getNumProcesses
      procedure :: processRank
      procedure :: sum
      procedure :: gatherString
      procedure :: gatherInteger
      procedure :: gatherLogical
      procedure :: allReduce
!TODO - NAG does not yet support FINAL keyword
!!$$      final :: clean
   end type SerialContext

   type (SerialContext), parameter :: THE_SERIAL_CONTEXT = SerialContext()

contains

   function newSerialContext() result(context)
      type (SerialContext) :: context
   end function newSerialContext

   integer function getNumProcesses(this)
      class (SerialContext),  intent(in) :: this

      getNumProcesses = 1

   end function getNumProcesses

   integer function processRank(this)
      class (SerialContext),  intent(in) :: this
      processRank = 0
   end function processRank

   integer function sum(this, value)
      class (SerialContext), intent(in) :: this
      integer, intent(in) :: value

      sum = value

   end function sum

   subroutine gatherString(this, values, list)
      class (SerialContext), intent(in) :: this
      character(len=*), intent(in) :: values(:)
      character(len=*), intent(out) :: list(:)

      list = values
   end subroutine gatherString

   subroutine gatherInteger(this, values, list)
      class (SerialContext), intent(in) :: this
      integer, intent(in) :: values(:)
      integer, intent(out) :: list(:)

      list = values

   end subroutine gatherInteger

   subroutine gatherLogical(this, values, list)
      class (SerialContext), intent(in) :: this
      logical, intent(in) :: values(:)
      logical, intent(out) :: list(:)

      list = values
   end subroutine gatherLogical

   logical function allReduce(this, q) result(anyQ)
      class (SerialContext), intent(in) :: this
      logical, intent(in) :: q
      anyQ = q
   end function allReduce

   subroutine clean(this)
      type (SerialContext), intent(inout) :: this
   end subroutine clean

end module SerialContext_mod
