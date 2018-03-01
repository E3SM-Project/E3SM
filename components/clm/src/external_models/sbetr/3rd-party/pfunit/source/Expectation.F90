

! Note: maybe have multiple expectation types for subroutines, classes, etc.
! 

module Expectation_mod
  use StringConversionUtilities_mod, only : MAXLEN_STRING
  implicit none
  private

  public :: Expectation, newExpectation
  public :: Predicate, newPredicate
  public :: Subject, newSubject, newSubjectNameOnly
  public :: wasCalled, wasNotCalled, wasCalledOnce

  type :: Subject
     ! mlr todo allocatable strings
     character(len=MAXLEN_STRING) :: name
     procedure(subVoid), pointer, nopass :: ptr
  end type Subject

  interface 
     subroutine subVoid
     end subroutine subVoid
  end interface


  type :: Predicate
     character(len=MAXLEN_STRING) :: name
  end type Predicate

! TDD
  type(Predicate), parameter :: wasCalled     = Predicate('wasCalled')
  type(Predicate), parameter :: wasNotCalled  = Predicate('wasNotCalled')
  type(Predicate), parameter :: wasCalledOnce = Predicate('wasCalledOnce')
! todo:  
!    checking expectation sub called with right value (important for sci.)
!    syntax for distinguishing arguments -- (position/keys)
!    combined expectations -- one on method, one on argument
!    -- or combined in the text...
! todo expectation augment
!    - vary numbers & kinds of arguments 
! todo:  automatic generation -- for proposal
! todo:  a trivial example of interleaved method calls
! 
! todo question: !    how to require mock functions to return certain values

  type :: Expectation
     type(Subject) :: subj
     type(Predicate) :: pred
  end type Expectation

contains

  type(Predicate) function newPredicate(name) result(pred_)
    character(*) :: name
    pred_%name = name
  end function newPredicate

  type(Subject) function newSubject(name,sub) result(subj_)
    character(*) :: name
    procedure(subVoid), pointer :: sub
    subj_%name = name
    subj_%ptr => sub
    ! maybe include a reference too
  end function newSubject

  type(Subject) function newSubjectNameOnly(name) result(subj_)
    character(*) :: name
    procedure(subVoid), pointer :: sub
    subj_%name = name
    ! subj_%ptr => sub ! Maybe nullify...
    nullify(subj_%ptr)
    ! maybe include a reference too
  end function newSubjectNameOnly

!  type(Subject) function newSubject(name) result(subj_)

  type(Expectation) function newExpectation(subj, pred) result(exp_)
    type(Subject), intent(in) :: subj
    type(Predicate), intent(in) :: pred
    exp_%subj = subj
    exp_%pred = pred
  end function newExpectation

end module Expectation_mod
