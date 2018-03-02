!#include "reflection.h"

!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_AssertComplex_mod
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune,  NASA/GSFC
!!
!! @date
!! 20 Mar 2015
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 20 Mar 2015 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------

!
! Test_AssertComplex_mod
!    Were tests for AssertReal using TestAssertReal.F90 as a model.
!
!    2013-0425 MLR: First version. 
!
!
! Note: `note name' is a mark to remind me of references to code being
!    tested.  I.e. non-boilerplate.

module Test_AssertComplex_mod ! note name
!  use Exception_mod, only: getNumExceptions, anyExceptions
  use TestSuite_mod
  use Params_mod, only : r32, i64, i32
  use StringConversionUtilities_mod, only: toString, appendWithSpace
  use AssertBasic_mod
  use Assert_mod
  use AssertArraysSupport_mod, only: differenceReport, valuesReport
!   AssertReal_mod, only: assertEqual, differenceReport, valuesReport
!   AssertComplex_mod, only: assertEqual
!   AssertComplex_mod, only: assertNotEqual
!   AssertComplex_mod, only: assertRelativelyEqual
  use ThrowFundamentalTypes_mod, only: locationFormat
! , differenceReport, valuesReport
  use SourceLocation_mod

  implicit none
  private

  public :: suite

  complex(kind=r32), parameter :: good = (42.0,24.0)
  complex(kind=r32), parameter :: bad  = (-666,-999)

contains

  function suite()
    use TestSuite_mod, only: TestSuite, newTestSuite 
    use TestMethod_mod, only: newTestMethod

    type (TestSuite) :: suite

    suite = newTestSuite('AssertComplexSuite') 

!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

    call suite%addTest( &
           &   newTestMethod('testEquals_C_complexScalar', &
           &                  testEquals_C_complexScalar))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_0D1D', &
           &                  testEquals_C_0D1D))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_1D_nonConformable1', &
           &                  testEquals_C_1D_nonConformable1))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_2D_SingleElementDifferent', &
           &                  testEquals_C_2D_SingleElementDifferent))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiD_SingleElementDifferent', &
           &                  testEquals_C_MultiD_SingleElementDifferent))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiD_SingleElementDifferent1', &
           &                  testEquals_C_MultiD_SingleElementDifferent1))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiD_SingleElementDifferent2', &
           &                  testEquals_C_MultiD_SingleElementDifferent2))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiD_SingleElementDifferent3', &
           &                  testEquals_C_MultiD_SingleElementDifferent3))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiD_SingleElementDifferent4', &
           &                  testEquals_C_MultiD_SingleElementDifferent4))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiD_SingleElementDifferent5', &
           &                  testEquals_C_MultiD_SingleElementDifferent5))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff1', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff1))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff2', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff2))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff3', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff3))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff4', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff4))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff5', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff5))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff6', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff6))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff7', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff7))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDMultiPrec_SingleEltDiff8', &
           &                  testEquals_C_MultiDMultiPrec_SingleEltDiff8))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDWithTolerance', &
           &                  testEquals_C_MultiDWithTolerance))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDWithTolerance1', &
           &                  testEquals_C_MultiDWithTolerance1))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDWithTolerance64', &
           &                  testEquals_C_MultiDWithTolerance64))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDWithTolerance64_1', &
           &                  testEquals_C_MultiDWithTolerance64_1))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDWithTolerance64_2', &
           &                  testEquals_C_MultiDWithTolerance64_2))
    call suite%addTest( &
           &   newTestMethod('testEquals_C_MultiDSourceLocation', &
           &                  testEquals_C_MultiDSourceLocation))
    call suite%addTest( &
           &   newTestMethod('testEquals_4DPComplex_DifferenceReport', &
           &                  testEquals_4DPComplex_DifferenceReport))
    call suite%addTest( &
           &   newTestMethod('testEquals_ComplexMultiD_SingleElementNE1', &
           &                  testEquals_ComplexMultiD_SingleElementNE1))
    call suite%addTest( &
           &   newTestMethod('testEquals_ComplexMultiD_SingleElementRE1', &
           &                  testEquals_ComplexMultiD_SingleElementRE1))
    call suite%addTest( &
           &   newTestMethod('testEquals_ComplexMultiD_SingleEltVarious1', &
           &                  testEquals_ComplexMultiD_SingleEltVarious1))

  end function suite


  subroutine testEquals_C_complexScalar()
    use Params_mod

    complex(kind=r32) :: expected
    complex(kind=r32) :: found

    character(len=:), allocatable :: msg

    expected = good
    found = bad

    allocate(msg,source='testEquals_C_complexScalar')

    call assertEqual(expected, found, message=msg)

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & '.' ) )

!mlr-         & ';  first difference at element  [1].') &
!mlr-         & )
    
    deallocate(msg)

  end subroutine testEquals_C_complexScalar


  ! Same rank, different shape.
  subroutine testEquals_C_0D1D()
    use Params_mod
!    use Assert_mod, only: assertEqual

    integer(kind=i64) :: expected
    integer, parameter :: good = 42
    complex(kind=r32), parameter :: z_good = good
    complex(kind=r32), dimension(1) :: found

    character(len=:), allocatable :: msg

    expected = good
    ! The location [1] below is the 1 here.
    found = bad

    ! The following should throw an exception...
    
    allocate(msg,source='testEquals_C_0D1D')

    call assertEqual(expected, found, message=msg)

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - z_good), 0.)) // &
         & ';  first difference at element [1].') &
         & )
    
    deallocate(msg)

  end subroutine testEquals_C_0D1D

  ! Same rank, different shape.
  subroutine testEquals_C_1D_nonConformable1()
    use Params_mod
!    use Assert_mod, only: assertEqual


    integer(kind=i32), dimension(2) :: expected
    complex(kind=r32), dimension(1) :: found

    character(len=:), allocatable :: msg

    ! The following should throw an exception...

    expected = good; found = good
    
    allocate(msg,source='testEquals_C_2D_nonConformable1')

    call assertEqual(expected, found, message=msg)

    call assertCatch( &
         & appendWithSpace(msg, &
         & 'nonconforming arrays - expected shape: ' // &
         & trim(toString(shape(expected))) // ' but found shape: ' // &
         & trim(toString(shape(found)))) &
         & )
    
    deallocate(msg)

  end subroutine testEquals_C_1D_nonConformable1

  subroutine testEquals_C_2D_SingleElementDifferent()
    use Params_mod
!    use Assert_mod, only: assertEqual


    complex, dimension(2,2) :: expected, found

    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2

    character(len=:), allocatable :: msg

    i1 = 1; i2 = 2; expected=good; found=good; found(i1,i2) = bad

    !dbg1 print *,'1000'

    allocate(msg,source='testEquals_C_2D_SingleElementDifferent')

    ! The following should throw an exception...
!    call assertEqual(expected,found, message=msg)
    call assertEqual(expected,found,msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good,bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_2D_SingleElementDifferent

  subroutine testEquals_C_MultiD_SingleElementDifferent()
    use Params_mod
!    use Assert_mod, only: assertEqual

    real(kind=r32) :: expected
    complex(kind=r32), dimension(:,:), allocatable :: found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2

    character(len=:), allocatable :: msg

    !dbg1 print *,'2000'

    n1 = 1; n2 = 2; allocate(found(n1,n2))
    expected=good; found=expected
    i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg, source='testEquals_C_MultiD_SingleElementDifferent:Rank0')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(real(good), bad)) // &
         & '; ' // trim(differenceReport(abs(bad - real(good)), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiD_SingleElementDifferent

  subroutine testEquals_C_MultiD_SingleElementDifferent1
    use Params_mod
!    use Assert_mod, only: assertEqual

! Don't do ths.
!    real(kind=r32), dimension(:,:), allocatable :: found
!    complex(kind=r32), dimension(:,:), allocatable :: expected

!  And if I just switch them...
    real(kind=r32), dimension(:,:), allocatable :: expected
    complex(kind=r32), dimension(:,:), allocatable :: found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2

    character(len=:), allocatable :: msg

    !dbg1 print *,'3000'

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    ! Note interesting interplay if we use found = expected...
    expected = real(good); found = real(good); i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank2')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(real(good),bad)) // &
         & '; ' // trim(differenceReport(abs(bad - real(good)), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiD_SingleElementDifferent1

  subroutine testEquals_C_MultiD_SingleElementDifferent2
    use Params_mod
!    use Assert_mod, only: assertEqual

! Don't do this...
!    complex(kind=r32), dimension(:,:,:), allocatable :: expected
!    real(kind=r32), dimension(:,:,:), allocatable :: found

! Try a simple switch...
    complex(kind=r32), dimension(:,:,:), allocatable :: found
    real(kind=r32), dimension(:,:,:), allocatable :: expected


    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3

    character(len=:), allocatable :: msg

    !dbg1 print *,'4000'

    n1 = 2; n2 = 3; n3 = 1; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good; found = real(good);
    i1 = 1; i2 = 1; i3 = 1; found(i1,i2,i3) = good

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank3')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3] )) [i1, i2, i3]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(real(good), good)) // &
         & '; ' // trim(differenceReport(abs(real(good) - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiD_SingleElementDifferent2

  subroutine testEquals_C_MultiD_SingleElementDifferent3
    use Params_mod
!    use Assert_mod, only: assertEqual

    complex(kind=r32), dimension(:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    character(len=:), allocatable :: msg

    !dbg1 print *,'5000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good; found = good;
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank4')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4] )) [i1, i2, i3, i4]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiD_SingleElementDifferent3

  subroutine testEquals_C_MultiD_SingleElementDifferent4
    use Params_mod
!    use Assert_mod, only: assertEqual

    complex(kind=r32), dimension(:,:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    character(len=:), allocatable :: msg

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5),found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank5')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiD_SingleElementDifferent4

  subroutine testEquals_C_MultiD_SingleElementDifferent5
    use Params_mod
!    use Assert_mod, only: assertEqual

    complex(kind=r32), dimension(:,:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    character(len=:), allocatable :: msg

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5))
    n1 = 2; n2 = 3; n3 = 2; n4 = 3; n5 = 2;
    allocate(found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    !dbg2 print *,'10000'

    allocate(msg,source=&
         & 'testEquals_C_MultiD_SingleElementDifferent:nonConformable')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code. Not needed for nonconf.
    ! write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
          & 'nonconforming arrays - expected shape: ' // &
          & trim(toString(shape(expected))) // ' but found shape: ' // &
          & trim(toString(shape(found)))) &
          & )

    deallocate(msg)

  end subroutine testEquals_C_MultiD_SingleElementDifferent5

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2

    character(len=:), allocatable :: msg

    !dbg3 print *,'11000'

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good
    i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank2')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )


    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff1()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3

    character(len=:), allocatable :: msg

    !dbg3 print *,'12000'

    n1 = 2; n2 = 3; n3 = 1; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; found(i1,i2,i3) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank3')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3] )) [i1, i2, i3]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff1

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff2()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    character(len=:), allocatable :: msg

    !dbg3 print *,'13000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good; found = good;
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank4')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4] )) [i1, i2, i3, i4]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )


    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff2

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff3()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    character(len=:), allocatable :: msg

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5),found(n1,n2,n3,n4,n5))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank5')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff3

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff4()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    character(len=:), allocatable :: msg


    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5))
    n1 = 2; n2 = 3; n3 = 2; n4 = 3; n5 = 2;
    allocate(found(n1,n2,n3,n4,n5))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source=&
         & 'testEquals_C_MultiD_SingleElementDifferent:Rank5:NonConformable')

    ! The following should throw an exception...
    call assertEqual(expected, found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code. Not needed for nonconf.
    ! write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
          & 'nonconforming arrays - expected shape: ' // &
          & trim(toString(shape(expected))) // ' but found shape: ' // &
          & trim(toString(shape(found)))) &
          & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff4

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff5()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2

    character(len=:), allocatable :: msg


    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good
    i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank2')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff5

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff6()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3

    character(len=:), allocatable :: msg

    n1 = 2; n2 = 3; n3 = 1; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; found(i1,i2,i3) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank3')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3] )) [i1, i2, i3]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff6

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff7()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    character(len=:), allocatable :: msg

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank4')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4] )) [i1, i2, i3, i4]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff7

  subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff8()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    character(len=:), allocatable :: msg

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5),found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source='testEquals_C_MultiD_SingleElementDifferent:Rank5')

    ! The following should throw an exception...
    call assertEqual(expected,found,message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(expected,found)
    deallocate(msg)

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5))
    n1 = 2; n2 = 3; n3 = 2; n4 = 3; n5 = 2;
    allocate(found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad


    allocate(msg, source = &
         & 'testEquals_C_MultiD_SingleElementDifferent:Rank5:NonConformable')

    ! The following should throw an exception...
    call assertEqual(expected,found, msg)

    ! "locationInArray" is not used in the original AssertEqual code. Not needed for nonconf.
    ! write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
          & 'nonconforming arrays - expected shape: ' // &
          & trim(toString(shape(expected))) // ' but found shape: ' // &
          & trim(toString(shape(found)))) &
          & )


    deallocate(msg)

  end subroutine testEquals_C_MultiDMultiPrec_SingleEltDiff8

  subroutine testEquals_ScalarWithTolerance()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r32) :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray

    real(kind=r32) :: tolerance32
    complex(kind=r32) :: bad32

    expected = good; found = good;

    tolerance32 = 0.01
    bad32 = good + tolerance32*2.0

    found = bad32

    allocate(msg,source= &
         & 'testEquals_ScalarWithTolerance')

    ! The following should throw an exception...
    ! call assertEqual(expected,found,tolerance = tolerance32, message=msg)
    call assertEqual(expected,found,tolerance32, message=msg)

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad32)) // &
         & '; ' // trim(differenceReport(abs(bad32 - good), tolerance32)) // &
         & '.' ))

!mlr-         & ';  first difference at element  [
!mlr-         & )

    deallocate(msg)

  end subroutine testEquals_ScalarWithTolerance

  subroutine testEquals_C_MultiDWithTolerance()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r32), dimension(:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    real(kind=r32) :: tolerance32, bad32

    character(len=:), allocatable :: msg

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good;

    tolerance32 = 0.01
    bad32 = good + tolerance32*2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad32

    allocate(msg,source=  &
         & 'testEquals_C_MultiDSingleEltTol32-Throw:Rank2,Tolerance32')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance =tolerance32, message = msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, cmplx(bad32))) // &
         & '; ' // trim(differenceReport(abs(bad32 - good), tolerance32)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDWithTolerance

  subroutine testEquals_C_MultiDWithTolerance1()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r32), dimension(:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    complex(kind=r32)    :: bad32
    real(kind=r32)       :: tolerance32

    character(len=:), allocatable :: msg

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good;

    tolerance32 = 0.01
    bad32 = good + tolerance32/2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad32

    allocate(msg,source= &
         & 'testEquals_C_MultiDSingleEltTol32-NoThrow:Rank2,Tolerance32')

    ! The following should not throw an exception...
    call assertEqual(expected,found,tolerance =tolerance32, message = msg) 

    call assertCatch( "" )

    deallocate(msg)

  end subroutine testEquals_C_MultiDWithTolerance1

  subroutine testEquals_C_MultiDWithTolerance64()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:), allocatable :: expected, found
    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    real(kind=r64)    :: tolerance64
    complex(kind=r64) :: good64, bad64

    character(len=:), allocatable :: msg

    good64 = good

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good;

    tolerance64 = 0.01
    bad64 = good64 + tolerance64*2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad64

    allocate(msg,source=&
         & 'testEquals_C_MultiDSingleEltTol64-Throw:Rank2,Tolerance64')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance=tolerance64, message =msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

! Fix the need for the real below.  Note we're just reporting at this stage, not calculating.
    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good64, bad64)) // &
         & '; ' // trim(differenceReport(abs(bad64 - good64), tolerance64)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

end subroutine testEquals_C_MultiDWithTolerance64

  subroutine testEquals_C_MultiDWithTolerance64_1()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    real(kind=r64)    :: tolerance64, good64, bad64

    character(len=:), allocatable :: msg

    good64 = good

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good64; found = good64;

    tolerance64 = 0.01
    bad64 = good64 + tolerance64/2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad64

    allocate(msg,source= &
         & 'testEquals_C_MultiDSingleEltTol64-NoThrow:Rank2,Tolerance64')

    ! The following should not throw an exception...
    call assertEqual(expected,found,tolerance=tolerance64, message=msg)

    call assertCatch( "" )

    deallocate(msg)

  end subroutine testEquals_C_MultiDWithTolerance64_1


  subroutine testEquals_C_MultiDWithTolerance64_2()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    complex(kind=r64), dimension(:,:,:), allocatable :: expected, found

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3
    real(kind=r64)    :: tolerance64
    complex(kind=r64) :: good64, bad64

    character(len=:), allocatable :: msg

    good64 = good

    n1 = 1; n2 = 2; n3 = 2; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good64; found = good64

    tolerance64 = 0.01
    bad64 = good64 + tolerance64*2.0

    i1 = 1; i2 = 2; i3 = 1
    found(i1,i2,i3) = bad64

    allocate(msg,source=&
         & 'testEquals_C_MultiDSingleEltTol64-Throw:Rank3,Tolerance64')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance=tolerance64,message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3] )) [i1, i2, i3]

    call assertCatch( &
         & appendWithSpace(msg,&
         & trim(valuesReport(good64, bad64)) // &
         & '; ' // trim(differenceReport(abs(bad64 - good64), tolerance64)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDWithTolerance64_2

  subroutine testEquals_C_MultiDSourceLocation()
    use Params_mod
    implicit none

    complex(kind=r64), dimension(:,:), allocatable :: expected, found
    complex(kind=r64) :: good64, bad64
    real(kind=r64) :: tolerance64

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    type (SourceLocation) :: location

    character(len=:), allocatable :: msg

    good64 = good
    
    n1 = 2; n2 = 3; allocate(expected(n1,n2),found(n1,n2))
    expected = good64; found = good64

    tolerance64 = 0.01
    bad64 = good64 + tolerance64*2.0

    i1 = 2; i2 = 3; found(i1,i2) = bad64

    location = SourceLocation(lineNumber=999,fileName='AFileName')


    allocate(msg,source= &
         & 'testEquals_C_MultiDSourceLocation')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance = tolerance64, message = msg, &
         & location=location)

    ! location = SourceLocation(lineNumber=998,fileName='AFileName2')

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

! Note use of real...  Consider overloading the reporting functions...
    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good64, bad64)) // &
         & '; ' // trim(differenceReport(abs(bad64 - good64), tolerance64)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.'), &
         & location=location &
         & )

    deallocate(msg)

  end subroutine testEquals_C_MultiDSourceLocation

  subroutine testEquals_4DPComplex_DifferenceReport()
    use Params_mod
    implicit none
    
    complex(kind=r64), dimension(4) :: expected, found
    real(kind=r64), parameter :: tolerance64 = 0.2220446E-13

    character(len=:), allocatable :: msg

    expected = [ (-2.18668712341731,-5.06286061261393), &
         &       (-0.434093364126748,-1.76502068116290), &
         &       (0.553931617203788,0.705041772163443), &
         &       (-0.833863952864297,-1.68080268623140) &
         &     ]

    found = [ (-2.18668712341731,-5.06286061261393), &
         &    (-0.434093364126748,-1.76502068116290), &
         &    (5.991912653851966E-002,-0.529989454499727), &
         &    (-0.833863952864297,-1.68080268623140) &
         &  ]

    found(1) = found(1) + (1.0E-15,0.0)


    allocate(msg,source= 'testEquals_4DPComplex_DifferenceReport')


    call assertEqual(expected,found,tolerance = tolerance64, &
         & message = msg )

    call assertCatch( &
         & appendWithSpace(msg,&
         & trim(valuesReport(expected(3),found(3))) // &
         & '; ' // trim(differenceReport(abs(found(3)-expected(3)),tolerance64)) // &
         & ';  first difference at element [3].') &
         & )

!    print *,'1000'//' DifferenceReportTesting: ' // &
!         & trim(differenceReport(abs(found(3)-expected(3)),tolerance64))
!    print *,'1001'//' DifferenceReportTesting: ' // &
!         & trim(differenceReport((found(3)-expected(3)),tolerance64))

    deallocate(msg)

  end subroutine testEquals_4DPComplex_DifferenceReport

  subroutine testEquals_ComplexMultiD_SingleElementNE1
    use Params_mod
!    use Assert_mod, only: assertEqual

    complex, parameter :: good = 1

    complex(kind=r32), dimension(:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    !dbg1 print *,'5000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good + 1 ; found = good;
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = expected(i1,i2,i3,i4)

    allocate(msg,source='testEquals_MultiD_SingleElementNE1')

    ! The following should not throw an exception...
    call assertNotEqual(expected,found, message=msg)

!    "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4] )) [i1, i2, i3, i4]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good+1, good+1, & 
         &      ePrefix='NOT: expected', &
         &      fPrefix='but found:' )) // &
         & '; '//trim(differenceReport(0.,0.))// &
         & ';  first equality at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_ComplexMultiD_SingleElementNE1

  subroutine testEquals_ComplexMultiD_SingleElementRE1
    use Params_mod
!    use Assert_mod, only: assertEqual

    complex, parameter :: good = 1

    complex(kind=r32), dimension(:,:,:,:), allocatable :: expected, found

    real(kind=r32) :: tolerance32,numerator,denominator

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    !dbg1 print *,'5000'

    tolerance32 = 0.001

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good ; found = good + tolerance32*0.5;
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = expected(i1,i2,i3,i4) + 1.0

    numerator   = abs(found(i1,i2,i3,i4) - expected(i1,i2,i3,i4) )
    denominator = abs(expected(i1,i2,i3,i4)) ! zero? 

    allocate(msg,source='testEquals_MultiD_SingleElementRE1')

    ! The following should not throw an exception...
    call assertRelativelyEqual(expected,found, tolerance32, message=msg)

!    "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4] )) [i1, i2, i3, i4]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(expected(i1,i2,i3,i4), found(i1,i2,i3,i4), & 
         &      ePrefix='RELEQ: expected', &
         &      fPrefix='to be near:' )) // &
         & '; '//trim(differenceReport(numerator/denominator,tolerance32))// &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_ComplexMultiD_SingleElementRE1


  subroutine testEquals_ComplexMultiD_SingleEltVarious1
    use Params_mod
!    use Assert_mod, only: assertEqual

    complex, parameter :: good = 1

    complex(kind=r32), dimension(:,:,:,:), allocatable :: expected, found
    real(kind=r32) :: tolerance32

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    !dbg1 print *,'5000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    i1 = 1; i2 = 2; i3 = 1; i4 = 2

    allocate(msg,source='testEquals_ComplexMultiD_SingleEltVarious1')

    ! The following should not throw an exception...

    expected = good; found = good;
    call assertEqual(expected, found, message=msg)

    expected = good+1; found = good;
    call assertNotEqual(expected, found, message=msg)

!    expected = good+1; found = good;
!    call assertGreaterThan(expected, found, message=msg)
!
!    expected = good+1; found = good; found(i1,i2,i3,i4) = good+1
!    call assertGreaterThanOrEqual(expected, found, message=msg)
!
!    expected = good; found = good+1;
!    call assertLessThan(expected, found, message=msg)
!
!    expected = good; found = good+1; found(i1,i2,i3,i4) = good
!    call assertLessThanOrEqual(expected, found, message=msg)

    tolerance32 = 0.001
    expected = good; found = good + tolerance32*0.5;
    call assertRelativelyEqual(expected, found, tolerance32, message=msg )

!    "locationInArray" is not used in the original AssertEqual code.
!    write(locationInArray,locationFormat( [i1,i2,i3,i4] )) [i1, i2, i3, i4]

!    call assertCatch( &
!         & appendWithSpace(msg, &
!         & trim(valuesReport(good+1, good, & 
!         &      ePrefix='expected', &
!         &      fPrefix='to be greater than:' )) // &
!         & ';  first difference at element ' // trim(locationInArray) // '.') &
!         & )

    deallocate(msg)

  end subroutine testEquals_ComplexMultiD_SingleEltVarious1


  ! Check to see that the test result is as expected...
  subroutine assertCatch(string,location)
    use Params_mod
    use Exception_mod, only: getNumExceptions, Exception, catchNext
    use Assert_mod, only: assertEqual
    character(len=*), intent(in) :: string
    type (SourceLocation), optional, intent(in) :: location
    type (Exception) :: anException

    if (getNumExceptions() > 0) then
       anException = catchNext()

       !, 'exceptions do not match')
       call assertEqual(string,anException%getMessage()) ! ,message='Exception message test')
       if(present(location))then
          call assertEqual( &
               & location%lineNumber,anException%getLineNumber(), &
               & message='Source line number test')
          call assertEqual(location%fileName,anException%getFileName(), &
               & message='Source file name test')
       end if
    else
       !, 'missing exception')
       call assertEqual(string, ' ')
    end if
  end subroutine assertCatch

  
end module Test_AssertComplex_mod

  
