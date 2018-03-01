!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: MockRepository
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
module MockRepository_mod
   use Expectation_mod, only : Expectation, newExpectation
   use Expectation_mod, only : Subject, newSubject, newSubjectNameOnly
   use Expectation_mod, only : Predicate, newPredicate
   use Expectation_mod, only : wasCalled, wasNotCalled, wasCalledOnce

   implicit none
   private

   public :: MockRepository
   public :: newMockRepository
   public :: MockRepositoryPointer
   public :: MAX_LEN_METHOD_NAME
   public :: MAX_LEN_CALL_REGISTRATION

   integer, parameter :: MAX_LEN_METHOD_NAME = 32
   integer, parameter :: MAX_LEN_CALL_REGISTRATION = 32
   type MockRepository
      ! mlr todo Allocatable strings are available...
      ! mlr todo test under nag
      character(len=MAX_LEN_METHOD_NAME) :: method = ' '
      type(Expectation) :: Expectations(2)
      integer :: lastExpectation = 0
      !   character(len=MAX_LEN_CALL_REGISTRATION) :: callRegistry(2) ! make dynamic...
      type(Expectation) :: callRegistry(2)
      integer :: lastRegistration = 0

   contains
      procedure :: verifyMocking

      procedure :: expectCall
      procedure :: hasCalled

      generic   :: addExpectationThat => addExpectationThat_sub_
      generic   :: addExpectationThat => addExpectationThat_subNameOnly_

      procedure :: addExpectationThat_sub_
      procedure :: addExpectationThat_subNameOnly_

      generic   :: registerMockCallBy => registerMockCallBy_subName_

      procedure :: registerMockCallBy_subName_

      procedure :: verify

      ! final?
      procedure :: delete

   end type MockRepository

   interface addExpectationThat_
      module procedure addExpectationThat_sub_
      module procedure addExpectationThat_subNameOnly_
   end interface addExpectationThat_

   interface registerMockCallBy_
      module procedure registerMockCallBy_subName_
   end interface registerMockCallBy_

   interface
      subroutine subVoid
      end subroutine subVoid
   end interface

!   interface 
!      module procedure addExpectationThat_sub_
!   end interface

   class (MockRepository), pointer :: MockRepositoryPointer => null()

contains

!! Begin older code

   function newMockRepository() result(repository)
      type (MockRepository), pointer :: repository
!      type (MockRepository), allocatable, target :: repository
      if ( associated(MockRepositoryPointer) ) then
         print *,'MockRepository::ERROR::RepositoryAlreadyAllocated'
      end if
      allocate(repository)
      MockRepositoryPointer => repository
      repository%lastExpectation  = 0
      repository%lastRegistration = 0
   end function newMockRepository

   subroutine delete (this)
     class (MockRepository), intent(inout) :: this
     
     nullify(MockRepositoryPointer)

   end subroutine delete

   subroutine verifyMocking(this, object)
      use Exception_mod
      class (MockRepository), intent(inout) :: this
      class (*) :: object
      
      if (trim(this%method) /= '') then
         call throw('Expected method not called: method1() on object of class MockSUT.')
      end if

      ! Only need to verify once. Finish it off...
      call this%delete()

   end subroutine verifyMocking

   subroutine expectCall(this, obj, method)
      class (MockRepository), intent(inout) :: this
      class(*), intent(in) :: obj
      character(len=*), intent(in) :: method

      this%method = method
   end subroutine expectCall

   subroutine hasCalled(this, obj, method)
      class (MockRepository), intent(inout) :: this
      class(*), intent(in) :: obj
      character(len=*), intent(in) :: method

      if (trim(method) == trim(this%method)) then
         this%method=''
      end if
   end subroutine hasCalled

!! End older code

   subroutine addExpectationThat_sub_(this,sub,pred)
     class (MockRepository), intent(inout) :: this
     procedure(subVoid), pointer, intent(in) :: sub
!     procedure(subVoid), pointer, intent(in) :: subptr
     type(Predicate), intent(in) :: pred
     type(Expectation) exp

     exp = newExpectation( &
          & newSubject( 'dummy-sub-name', sub), &
          & pred )

     ! <exp ok> & <enough space in list>
     this%lastExpectation = this%lastExpectation + 1
     this%Expectations(this%lastExpectation) = exp

   end subroutine addExpectationThat_sub_

   subroutine addExpectationThat_subNameOnly_(this,subName,pred)
     class (MockRepository), intent(inout) :: this
     character(len=*), intent(in) :: subName
!     procedure(subVoid), pointer, intent(in) :: sub
!     procedure(subVoid), pointer, intent(in) :: subptr
     type(Predicate), intent(in) :: pred
     type(Expectation) exp

     exp = newExpectation( &
          & newSubjectNameOnly( subName ), &
          & pred )

     ! <exp ok> & <enough space in list>
     this%lastExpectation = this%lastExpectation + 1
     this%Expectations(this%lastExpectation) = exp

   end subroutine addExpectationThat_subNameOnly_

   subroutine registerMockCallBy_subName_(this,subName)
     class (MockRepository), intent(inout) :: this
     character(len=*), intent(in) :: subName
!     print *,'reg200: ',subName,this%lastRegistration
     ! Can we includ the calling sub here? For a better comparison with our Exp. list?
     ! <if space in registry>
     this%lastRegistration = this%lastRegistration + 1
     this%callRegistry(this%lastRegistration) &
          & = newExpectation( & ! mlr todo Expectation --> foundAction -- "Result"
          &                  newSubjectNameOnly(subName), &
          &                  wasCalled )
   end subroutine registerMockCallBy_subName_

   subroutine verify(this)
      use Exception_mod
      class (MockRepository), intent(inout), target :: this
      integer iExp, iReg
      class (Expectation), pointer :: exp, reg
      logical ok
      integer nCalls

      ! Go through expectation logic. Note:  Maybe use original list of strings approach.
      ! Need to work out more complex logic.  Ess. need logic analyzer.
      ! Eventually, expectations should probably be trees.

      ! Maybe rework the following into "expected vs. found" for greater alignment
      ! with existing usage in PFUNIT.  Also consider existing capabilities in PFUNIT.

! mlr todo -- verify should not be hardwired for the things it has to handle 
! todo -- refactor expectations & foundActionResult -- recall interpreter implementation
!


      do iExp=1,this%lastExpectation  ! 'with' syntax?
         exp => this%Expectations(iExp)
         ok = .false.

!         if(exp%pred%name .eq. 'wasCalled')then
!            ok = .false.
!            do iReg=1,this%lastRegistration
!               reg = this%callRegistry(iReg)
!!               print *,'verify1000: ', &
!!                    & trim(exp%subj%name)//'='// &
!!                    & trim(reg%subj%name)//', '// &
!!                    & trim(exp%pred%name)//'='// &
!!                    & trim(reg%pred%name)//'.'
!               if(exp%subj%name .eq. reg%subj%name) then
!                  if(exp%pred%name .eq. reg%pred%name) then
!   !???               if(exp%pred .eq. reg%pred) then
!                     ok=.true.
!                  end if
!               end if
!            end do
!         end if
!
!         if(exp%pred%name .eq. 'wasNotCalled')then
!            ok = .true.
!            do iReg=1,this%lastRegistration
!               reg = this%callRegistry(iReg)
!               if(exp%subj%name .eq. reg%subj%name) then
!                  if(trim(reg%pred%name).eq.'wasCalled')then
!                     ok = .false.
!                  end if
!               end if
!            end do
!         end if

         if( &
              & (exp%pred%name .eq. 'wasCalled') .or. &
              & (exp%pred%name .eq. 'wasCalledOnce') .or. &
              & (exp%pred%name .eq. 'wasNotCalled') &
              & )then
            ok = .true.
            nCalls = 0
            do iReg=1,this%lastRegistration
               reg => this%callRegistry(iReg)
               if(exp%subj%name .eq. reg%subj%name) then
                  if(trim(reg%pred%name).eq.'wasCalled')then
                     nCalls=nCalls+1
                  end if
               end if
            end do
            if(exp%pred%name .eq. 'wasCalled')then
               ok = nCalls.ge.1
            else if(exp%pred%name .eq. 'wasCalledOnce')then
               ok = nCalls.eq.1
            else if(exp%pred%name .eq. 'wasNotCalled')then
               ok = nCalls.eq.0
            end if

         end if

         if(.not.ok)then
            call throw('Expectations not met:  "'// &
                 & trim(exp%subj%name)//'" "'//trim(exp%pred%name)//'" does not hold.')
         end if
      end do

      ! call throw('exception%verify: Not implemented')
      
    end subroutine verify

   
end module MockRepository_mod
