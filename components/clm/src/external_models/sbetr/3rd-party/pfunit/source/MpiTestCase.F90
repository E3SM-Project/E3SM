!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: MpiTestCase
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

module MpiTestCase_mod
   use MpiContext_mod
   use TestCase_mod
   use AbstractTestParameter_mod
   use MpiTestParameter_mod
   use ParameterizedTestCase_mod, only: ParameterizedTestCase
!!$   use ParameterizedTestCase_mod, only: MAX_LEN_LABEL
   implicit none
   private

   public :: MpiTestCase

   type, abstract, extends(ParameterizedTestCase) :: MpiTestCase
      integer :: processRank
      type (MpiContext) :: context
      type (MpiContext) :: parentContext
   contains
      procedure :: countTestCases => countTestCases_mpi
      procedure :: run
      procedure :: runBare
      procedure :: getNumProcesses
      procedure :: getNumProcessesRequested
      procedure :: getProcessRank
      procedure :: getMpiCommunicator
      procedure :: getContext
   end type MpiTestCase

contains

   integer function countTestCases_mpi(this) result(countTestCases)
      class (MpiTestCase), intent(in) :: this
      countTestCases = 1
   end function countTestCases_mpi

   recursive subroutine run(this, tstResult, context)
      use TestResult_mod, only: TestResult
      use ParallelContext_mod
      use Exception_mod
      use SurrogateTestCase_mod
      class (MpiTestCase), intent(inout) :: this
      class (TestResult), intent(inout) :: tstResult
      class (ParallelContext), intent(in) :: context

      select type (context)
      type is (MpiContext)
         this%parentContext = context
      class default
         call throw('MPI test cannot run in a non-MPI context.')
         return
      end select

      call tstResult%run(this%getSurrogate(), context)

   end subroutine run

   recursive subroutine runBare(this)
      use Exception_mod
      use ParallelException_mod
      class (MpiTestCase), intent(inout) :: this

      logical :: discard

      ! create subcommunicator
      this%context = this%parentContext%makeSubcontext(this%getNumProcessesRequested())

      if (.not. anyExceptions(this%parentContext)) then
         if (this%context%isActive()) then
            call this%setUp()

            if (.not. anyExceptions(this%context)) then
               call this%runMethod()
               call this%tearDown()
            end if
         end if
      else
         ! only report context failure on root PE
         if (.not. this%parentContext%isRootProcess()) then
            discard = catch()
         end if
      end if

      call gather(this%parentContext)

   end subroutine runBare

   integer function getMpiCommunicator(this) result(mpiCommunicator)
      class (MpiTestCase), intent(in) :: this
      mpiCommunicator = this%context%getMpiCommunicator()
   end function getMpiCommunicator

   integer function getNumProcesses(this) result(numProcesses)
      class (MpiTestCase), intent(in) :: this
      numProcesses = this%context%getNumProcesses()
   end function getNumProcesses

   integer function getNumProcessesRequested(this) result(numProcessesRequested)
      use Exception_mod
      class (MpiTestCase), intent(in) :: this
      select type (p => this%testParameter)
      class is (MpiTestParameter)
         numProcessesRequested = p%getNumProcessesRequested()
      class default
         call throw('Incorrect type of test parameter in MpiTestCase::getNumProcessesRequested()')
      end select
   end function getNumProcessesRequested

   integer function getProcessRank(this) result(processRank)
      class (MpiTestCase), intent(in) :: this
      processRank = this%context%processRank()
   end function getProcessRank

   function getContext(this) result(context)
      class (MpiTestCase), intent(in) :: this
      class (MpiContext), allocatable :: context

      allocate(context, source=this%context)

   end function getContext

end module MpiTestCase_mod
