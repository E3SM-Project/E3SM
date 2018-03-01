!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: TestCase
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
! Serial TestCase 
module TestCase_mod
   use Exception_mod,         only : throw
   use Params_mod,            only : MAX_LENGTH_NAME
   use SurrogateTestCase_mod
   use TestResult_mod
   use Test_mod

   private

   public :: TestCase
   public :: TestCaseReference

   type, extends(SurrogateTestCase) :: ConcreteSurrogate
      private
      class (TestCase), pointer :: tCase => null()
   contains
      procedure :: runBare => runBare_surrogate
      procedure :: setName => setName_surrogate
      procedure :: getName => getName_surrogate
   end type ConcreteSurrogate
   
   type, abstract, extends(Test) :: TestCase
      private
      type (ConcreteSurrogate) :: surrogate
#ifdef DEFERRED_LENGTH_CHARACTER
      character(:), allocatable :: name
#else
      character(len=MAX_LENGTH_NAME) :: name
#endif
   contains
      procedure :: setSurrogate
      procedure :: baseName
      procedure :: getName 
      procedure :: setName
      procedure :: countTestCases
      procedure :: run
      procedure :: runBare
      procedure :: setUp
      procedure :: tearDown
      procedure :: getSurrogate
      procedure :: runMethod
   end type TestCase

   type TestCaseReference
      class (TestCase), allocatable :: test
   end type TestCaseReference

contains

   function baseName(this) result(name)
      class (TestCase), intent(in) :: this
      character(:), allocatable :: name
      name = this%name
   end function baseName

   function getName(this) result(name)
      class (TestCase), intent(in) :: this
      character(:), allocatable :: name
      name = this%baseName()
   end function getName

   subroutine setName(this, name)
      class (TestCase), intent(inout) :: this
      character(len=*),intent(in) :: name

#ifndef DEFERRED_LENGTH_CHARACTER
      integer :: nameLength

      nameLength = len_trim( name )
      if (nameLength > MAX_LENGTH_NAME) then
         call throw( 'TestCase.setName: Too long: ' // name )
         nameLength = MAX_LENGTH_NAME
      end if
      this%name = name(1:nameLength)
#else
      this%name = trim(name)
#endif

   end subroutine setName

   integer function countTestCases(this)
      class (TestCase), intent(in) :: this
      countTestCases = 1
   end function countTestCases

! Implement deferred method from class Test
   recursive subroutine run(this, tstResult, context)
      use SerialContext_mod
      use TestResult_mod
      use ParallelContext_mod
      class (TestCase), intent(inout) :: this
      class (TestResult), intent(inout) :: tstResult
      class (ParallelContext), intent(in) :: context

      ! Always run serial tests in a serial context.
      if (context%isRootProcess()) then
         call tstResult%run(this%getSurrogate(), THE_SERIAL_CONTEXT)
      end if

      call context%barrier()

   end subroutine run

   recursive subroutine runBare(this)
      use Exception_mod, only: noExceptions
      class (TestCase), intent(inout) :: this

      call this%setUp()
      if (noExceptions()) then
         call this%runMethod()
         call this%tearDown()
      end if

   end subroutine runBare

   recursive subroutine runBare_surrogate(this)
      class (ConcreteSurrogate), intent(inout) :: this
      class (TestCase), pointer :: p
      p => this%tCase
      call p%runBare()
   end subroutine runBare_surrogate

   function getName_surrogate(this) result(name)
      class (ConcreteSurrogate), intent(in) :: this
      character(:), allocatable :: name
      name = this%tCase%getName()
   end function getName_surrogate

   subroutine setName_surrogate(this, name)
      class (ConcreteSurrogate), intent(inout) :: this
      character(len=*),intent(in) :: name
      call this%tCase%setName(trim(name))
   end subroutine setName_surrogate

   subroutine setUp(this)
      class (TestCase), intent(inOut) :: this
   end subroutine setUp

   subroutine tearDown(this)
      class (TestCase), intent(inOut) :: this
   end subroutine tearDown

   function getSurrogate(this) result(surrogate)
      class (TestCase), target, intent(inout) :: this
      class (SurrogateTestCase), pointer :: surrogate
      call this%setSurrogate()
      surrogate => this%surrogate
   end function getSurrogate

   subroutine setSurrogate(this)
      class (TestCase), target :: this
      this%surrogate%tCase => this
   end subroutine setSurrogate

   recursive subroutine runMethod(this)
      use Exception_mod, only: throw
      class (TestCase), intent(inout) :: this
      call throw('TestCase::runMethod() must be overridden.')
   end subroutine runMethod

end module TestCase_mod
