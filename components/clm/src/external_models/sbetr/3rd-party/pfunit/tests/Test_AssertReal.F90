!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_AssertReal_mod
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
! Test_AssertReal_mod
!    Tests for AssertReal using TestAssertReal.F90 as a model.
!
!    2013-0304 MLR: First version. 
!
!
! Note: `note name' is a mark to remind me of references to code being
!    tested.  I.e. non-boilerplate.


module Test_AssertReal_mod ! note name
  use TestSuite_mod
  use Params_mod, only : r32, i64
  use StringConversionUtilities_mod, only: toString, appendWithSpace
  use AssertBasic_mod
  use Assert_mod, only: assertEqual
  use Assert_mod, only: assertGreaterThan
  use Assert_mod, only: assertGreaterThanOrEqual
  use Assert_mod, only: assertLessThan
  use Assert_mod, only: assertLessThanOrEqual
  use Assert_mod, only: assertRelativelyEqual
  use AssertArraysSupport_mod, only: differenceReport, valuesReport
  use ThrowFundamentalTypes_mod, only: locationFormat
  use SourceLocation_mod

  implicit none
  private

  public :: suite

  real(kind=r32), parameter :: good = 1
  real(kind=r32), parameter :: bad  = -999

contains

  function suite()
    use TestSuite_mod, only: TestSuite, newTestSuite 
    use TestMethod_mod, only: newTestMethod

    type (TestSuite) :: suite

    suite = newTestSuite('AssertRealSuite') 

!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

    call suite%addTest( &
           &   newTestMethod('testEquals_0D1D', &
           &                  testEquals_0D1D))
    call suite%addTest( &
           &   newTestMethod('testEquals_1D_nonConformable1', &
           &                  testEquals_1D_nonConformable1))
    call suite%addTest( &
           &   newTestMethod('testEquals_2D_SingleElementDifferent', &
           &                  testEquals_2D_SingleElementDifferent))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementDifferent', &
           &                  testEquals_MultiD_SingleElementDifferent))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementDifferent1', &
           &                  testEquals_MultiD_SingleElementDifferent1))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementDifferent2', &
           &                  testEquals_MultiD_SingleElementDifferent2))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementDifferent3', &
           &                  testEquals_MultiD_SingleElementDifferent3))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementDifferent4', &
           &                  testEquals_MultiD_SingleElementDifferent4))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementDifferent5', &
           &                  testEquals_MultiD_SingleElementDifferent5))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff1', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff1))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff2', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff2))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff3', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff3))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff4', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff4))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff5', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff5))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff6', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff6))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff7', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff7))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDMultiPrec_SingleEltDiff8', &
           &                  testEquals_MultiDMultiPrec_SingleEltDiff8))
    call suite%addTest( &
           &   newTestMethod('testEquals_ScalarWithTolerance', &
           &                  testEquals_ScalarWithTolerance))
    call suite%addTest( &
           &   newTestMethod('testEquals_ScalarWithToleranceNoMsg', &
           &                  testEquals_ScalarWithToleranceNoMsg))
    call suite%addTest( &
           &   newTestMethod('testEquals_VectorWithToleranceNoMsg', &
           &                  testEquals_VectorWithToleranceNoMsg))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDWithTolerance', &
           &                  testEquals_MultiDWithTolerance))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDWithTolerance1', &
           &                  testEquals_MultiDWithTolerance1))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDWithTolerance64', &
           &                  testEquals_MultiDWithTolerance64))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDWithTolerance64_1', &
           &                  testEquals_MultiDWithTolerance64_1))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDWithTolerance64_2', &
           &                  testEquals_MultiDWithTolerance64_2))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiDSourceLocation', &
           &                  testEquals_MultiDSourceLocation))
    call suite%addTest( &
           &   newTestMethod('testEquals_ScalarAndLocation', &
           &                  testEquals_ScalarAndLocation))
    call suite%addTest( &
           &   newTestMethod('testEquals_ScalarInfinity_equal', &
           &                  testEquals_ScalarInfinity_equal))
    call suite%addTest( &
           &   newTestMethod('testEquals_ScalarInfinity_unequal_A', &
           &                  testEquals_ScalarInfinity_unequal_A))
    call suite%addTest( &
           &   newTestMethod('testEquals_ScalarInfinity_unequal_B', &
           &                  testEquals_ScalarInfinity_unequal_B))
    call suite%addTest( &
           &   newTestMethod('testEquals_ScalarInfinity_unequal_C', &
           &                  testEquals_ScalarInfinity_unequal_C))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementGT1', &
           &                  testEquals_MultiD_SingleElementGT1))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleElementGT2', &
           &                  testEquals_MultiD_SingleElementGT2))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleEltVarious1', &
           &                  testEquals_MultiD_SingleEltVarious1))
    call suite%addTest( &
           &   newTestMethod('testEquals_MultiD_SingleEltVarious2', &
           &                  testEquals_MultiD_SingleEltVarious2))

  end function suite

  ! Same rank, different shape.
  subroutine testEquals_0D1D()
    use Params_mod
!    use Assert_mod, only: assertEqual

!    integer :: expected
    integer(kind=i64) :: expected
    real(kind=r32), dimension(1) :: found

    character(len=:), allocatable :: msg

    expected = good
    ! The location [1] below is the 1 here.
    found = bad

    ! The following should throw an exception...

    allocate(msg, source='testEquals_0D1D')
    
    call assertEqual(expected, found, msg)

    ! Using expected below instead of good to get the type correct. Only affects printout.
    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(expected, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element [1].') &
         & )

    deallocate(msg)
    
  end subroutine testEquals_0D1D

  ! Same rank, different shape.
  subroutine testEquals_1D_nonConformable1()
    use Params_mod
!    use Assert_mod, only: assertEqual


    integer, dimension(2) :: expected
    real(kind=r32), dimension(1) :: found

    character(len=:), allocatable :: msg

    allocate(msg, source='testEquals_2D_nonConformable1')

    ! The following should throw an exception...

    expected = good; found = good
    
    call assertEqual(expected, found, msg)

    call assertCatch( &
         & appendWithSpace(msg, & 
         & 'nonconforming arrays - expected shape: ' // &
         & trim(toString(shape(expected))) // ' but found shape: ' // &
         & trim(toString(shape(found)))) &
         & )

    deallocate(msg)
    
  end subroutine testEquals_1D_nonConformable1

  subroutine testEquals_2D_SingleElementDifferent()
    use Params_mod
!    use Assert_mod, only: assertEqual


    real, dimension(2,2) :: expected, found

    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2

    character(len=:), allocatable :: msg

    i1 = 1; i2 = 2; expected=good; found=good; found(i1,i2) = bad

    !dbg1 print *,'1000'

    allocate(msg, source='testEquals_2D_SingleElementDifferent')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good,bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

  end subroutine testEquals_2D_SingleElementDifferent

  subroutine testEquals_MultiD_SingleElementDifferent()
    use Params_mod
!    use Assert_mod, only: assertEqual

    real(kind=r32) :: expected
    real(kind=r32), dimension(:,:), allocatable :: found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2

    !dbg1 print *,'2000'

    n1 = 1; n2 = 2; allocate(found(n1,n2))
    expected=good; found=good
    i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg, source='testEquals_MultiD_SingleElementDifferent:Rank0')

    ! The following should throw an exception...
    call assertEqual(expected,found,message=msg)
    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_MultiD_SingleElementDifferent

  subroutine testEquals_MultiD_SingleElementDifferent1
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r32), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2

    !dbg1 print *,'3000'

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good; i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg, source='testEquals_MultiD_SingleElementDifferent:Rank2')

    ! The following should throw an exception...
    call assertEqual(expected,found, message=msg )

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad)) // &
         & '; ' // trim(differenceReport(abs(bad - good), 0.)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_MultiD_SingleElementDifferent1

  subroutine testEquals_MultiD_SingleElementDifferent2
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r32), dimension(:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3

    !dbg1 print *,'4000'

    n1 = 2; n2 = 3; n3 = 1; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good; found = good;
    i1 = 1; i2 = 2; i3 = 1; found(i1,i2,i3) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank3')

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

  end subroutine testEquals_MultiD_SingleElementDifferent2

  subroutine testEquals_MultiD_SingleElementDifferent3
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r32), dimension(:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    !dbg1 print *,'5000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good; found = good;
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank4')

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

  end subroutine testEquals_MultiD_SingleElementDifferent3

  subroutine testEquals_MultiD_SingleElementDifferent4
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r32), dimension(:,:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5),found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank5')

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

  end subroutine testEquals_MultiD_SingleElementDifferent4

  subroutine testEquals_MultiD_SingleElementDifferent5
    use Params_mod
!    use Assert_mod, only: assertEqual

    real(kind=r32), dimension(:,:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5))
    n1 = 2; n2 = 3; n3 = 2; n4 = 3; n5 = 2;
    allocate(found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    !dbg2 print *,'10000'

    allocate(msg,source=&
         & 'testEquals_MultiD_SingleElementDifferent:nonConformable')

    ! The following should throw an exception...
    call assertEqual(expected,found,message=msg)

    ! "locationInArray" is not used in the original AssertEqual code. Not needed for nonconf.
    ! write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
          & 'nonconforming arrays - expected shape: ' // &
          & trim(toString(shape(expected))) // ' but found shape: ' // &
          & trim(toString(shape(found)))) &
          & )

    deallocate(msg)

  end subroutine testEquals_MultiD_SingleElementDifferent5

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r64), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2

    !dbg3 print *,'11000'

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good
    i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank2')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff1()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r64), dimension(:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3

    !dbg3 print *,'12000'

    n1 = 2; n2 = 3; n3 = 1; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; found(i1,i2,i3) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank3')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff1

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff2()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r64), dimension(:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    !dbg3 print *,'13000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good; found = good;
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank4')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff2

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff3()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real, parameter :: good = 1
    real, parameter :: bad  = -999

    real(kind=r64), dimension(:,:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5),found(n1,n2,n3,n4,n5))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank5')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff3

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff4()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5


    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5))
    n1 = 2; n2 = 3; n3 = 2; n4 = 3; n5 = 2;
    allocate(found(n1,n2,n3,n4,n5))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source=&
         & 'testEquals_MultiD_SingleElementDifferent:Rank5:NonConformable')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff4

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff5()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2


    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good
    i1 = 1; i2 = 2; found(i1,i2) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank2')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff5

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff6()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3

    n1 = 2; n2 = 3; n3 = 1; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; found(i1,i2,i3) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank3')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff6

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff7()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good; found = good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank4')

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

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff7

  subroutine testEquals_MultiDMultiPrec_SingleEltDiff8()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4, i5
    integer :: n1, n2, n3, n4, n5

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5),found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementDifferent:Rank5')

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

    deallocate(expected,found)
    n1 = 2; n2 = 3; n3 = 2; n4 = 2; n5 = 2;
    allocate(expected(n1,n2,n3,n4,n5))
    n1 = 2; n2 = 3; n3 = 2; n4 = 3; n5 = 2;
    allocate(found(n1,n2,n3,n4,n5))
    expected=good; found=good
    i1 = 1; i2 = 2; i3 = 1; i4 = 2; i5 = 1
    found(i1,i2,i3,i4,i5) = bad

    deallocate(msg)
    allocate(msg, source= &
         & 'testEquals_MultiD_SingleElementDifferent:Rank5:NonConformable')

    ! The following should throw an exception...
    call assertEqual(expected,found,message=msg)

    ! "locationInArray" is not used in the original AssertEqual code. Not needed for nonconf.
    ! write(locationInArray,locationFormat( [i1,i2,i3,i4,i5] )) [i1, i2, i3, i4, i5]

    call assertCatch( &
         & appendWithSpace(msg, &
         & 'nonconforming arrays - expected shape: ' // &
         & trim(toString(shape(expected))) // ' but found shape: ' // &
         & trim(toString(shape(found)))) &
         & )

    deallocate(msg)

  end subroutine testEquals_MultiDMultiPrec_SingleEltDiff8

  subroutine testEquals_ScalarWithTolerance()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r32) :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray

    real(kind=r32) :: tolerance32, bad32

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

!mlr-         & ';  first difference at element  [1].') &
!mlr-         & )

    deallocate(msg)

  end subroutine testEquals_ScalarWithTolerance

  subroutine testEquals_ScalarWithToleranceNoMsg()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

!    real(kind=r32) :: expected, found, tolerance32
    real :: expected, found, tolerance32

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray

    real(kind=r32) :: bad32

    expected = good; found = good;

    tolerance32 = 0.01
    bad32 = good + tolerance32*2.0

    found = bad32

    allocate(msg,source='')

    ! The following should throw an exception...
    ! call assertEqual(expected,found,tolerance = tolerance32, message=msg)
    ! call assertEqual(expected,found,tolerance32, message=msg)
    call assertEqual(expected,found,tolerance32)

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad32)) // &
         & '; ' // trim(differenceReport(abs(bad32 - good), tolerance32)) // &
         & '.' ))

!mlr-         & ';  first difference at element  [1].') &
!mlr-         & )

    deallocate(msg)

  end subroutine testEquals_ScalarWithToleranceNoMsg

  subroutine testEquals_VectorWithToleranceNoMsg()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r32), dimension(3) :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray

    real(kind=r32) :: tolerance32, bad32

    expected = good; found = good;

    tolerance32 = 0.01
    bad32 = good + tolerance32*2.0

    found(2) = bad32

    allocate(msg,source='')

    ! The following should throw an exception...
    ! call assertEqual(expected,found,tolerance = tolerance32, message=msg)
    ! call assertEqual(expected,found,tolerance32, message=msg)
    call assertEqual(expected,found,tolerance32)

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad32)) // &
         & '; ' // trim(differenceReport(abs(bad32 - good), tolerance32)) // &
         & ';  first difference at element [2].') &
         & )

    deallocate(msg)

  end subroutine testEquals_VectorWithToleranceNoMsg

  subroutine testEquals_MultiDWithTolerance()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r32), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    real(kind=r32) :: tolerance32, bad32

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good;

    tolerance32 = 0.01
    bad32 = good + tolerance32*2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad32

    allocate(msg,source= &
         & 'testEquals_MultiDSingleEltTol32-Throw:Rank2,Tolerance32')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance = tolerance32, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good, bad32)) // &
         & '; ' // trim(differenceReport(abs(bad32 - good), tolerance32)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_MultiDWithTolerance


  subroutine testEquals_MultiDWithTolerance1()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r32), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    real(kind=r32)    :: tolerance32, bad32

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good;

    tolerance32 = 0.01
    bad32 = good + tolerance32/2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad32

    allocate(msg,source= &
         & 'testEquals_MultiDSingleEltTol32-NoThrow:Rank2,Tolerance32')

    ! The following should not throw an exception...
    call assertEqual(expected,found,tolerance = tolerance32, message=msg)

    deallocate(msg)

  end subroutine testEquals_MultiDWithTolerance1

  subroutine testEquals_MultiDWithTolerance64()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg
    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    real(kind=r64)    :: tolerance64, good64, bad64

    good64 = good

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good; found = good;

    tolerance64 = 0.01
    bad64 = good64 + tolerance64*2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad64

    allocate(msg,source= &
         & 'testEquals_MultiDSingleEltTol64-Throw:Rank2,Tolerance64')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance = tolerance64, message=msg)


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

end subroutine testEquals_MultiDWithTolerance64

  subroutine testEquals_MultiDWithTolerance64_1()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    real(kind=r64)    :: tolerance64, good64, bad64

    good64 = good

    n1 = 1; n2 = 2; allocate(expected(n1,n2),found(n1,n2))
    expected = good64; found = good64;

    tolerance64 = 0.01
    bad64 = good64 + tolerance64/2.0

    i1 = 1; i2 = 2; 
    found(i1,i2) = bad64

    allocate(msg,source= &
         & 'testEquals_MultiDSingleEltTol64-NoThrow:Rank2,Tolerance64')

    ! The following should not throw an exception...
    call assertEqual(expected,found,tolerance = tolerance64, message=msg)

    deallocate(msg)

  end subroutine testEquals_MultiDWithTolerance64_1


  subroutine testEquals_MultiDWithTolerance64_2()
    use Params_mod
!    use Assert_mod, only: assertEqual
    implicit none

    real(kind=r64), dimension(:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3
    integer :: n1, n2, n3
    real(kind=r64)    :: tolerance64, good64, bad64

    good64 = good

    n1 = 1; n2 = 2; n3 = 2; allocate(expected(n1,n2,n3),found(n1,n2,n3))
    expected = good64; found = good64

    tolerance64 = 0.01
    bad64 = good64 + tolerance64*2.0

    i1 = 1; i2 = 2; i3 = 1
    found(i1,i2,i3) = bad64

    allocate(msg,source= &
         & 'testEquals_MultiDSingleEltTol64-Throw:Rank3,Tolerance64')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance=tolerance64, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3] )) [i1, i2, i3]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good64,bad64)) // &
         & '; ' // trim(differenceReport(abs(bad64 - good64), tolerance64)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_MultiDWithTolerance64_2

  subroutine testEquals_MultiDSourceLocation()
    use Params_mod
    implicit none

    real(kind=r64), dimension(:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg
    real(kind=r64) :: tolerance64, good64, bad64

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2
    integer :: n1, n2
    type (SourceLocation) :: location

    good64 = good
    
    n1 = 2; n2 = 3; allocate(expected(n1,n2),found(n1,n2))
    expected = good64; found = good64

    tolerance64 = 0.01
    bad64 = good64 + tolerance64*2.0

    i1 = 2; i2 = 3; found(i1,i2) = bad64

    location = SourceLocation(lineNumber=999,fileName='AFileName')

    allocate(msg,source='testEquals_MultiDSourceLocation')

    ! The following should throw an exception...
    call assertEqual(expected,found,tolerance = tolerance64, message=msg, &
         & location=location)

    ! location = SourceLocation(lineNumber=998,fileName='AFileName2')

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2] )) [i1, i2]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good64, bad64)) // &
         & '; ' // trim(differenceReport(abs(bad64 - good64), tolerance64)) // &
         & ';  first difference at element ' // trim(locationInArray) // '.'), &
         & location=location &
         & )

    deallocate(msg)

  end subroutine testEquals_MultiDSourceLocation

  subroutine testEquals_ScalarAndLocation()
    use Params_mod
    implicit none

    real(kind=r64) :: expected, found

    character(len=:), allocatable :: msg
    real(kind=r64) :: tolerance64

    !mlr maybe move this to a larger scope...
!    integer, parameter :: MAXLEN_SHAPE = 80
!    character(len=MAXLEN_SHAPE) :: locationInArray
!    type (SourceLocation) :: location

    tolerance64 = 0.0
    expected = 2.0
    found = 3.0

!    location = SourceLocation(lineNumber=999,fileName='AFileName')

    allocate(msg,source='')

    call assertEqual(expected,found,msg)

    ! location = SourceLocation(lineNumber=998,fileName='AFileName2')
    ! "locationInArray" is not used in the original AssertEqual code.

! Note use of real...  Consider overloading the reporting functions...
    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(expected, found)) // &
         & '; ' // trim(differenceReport(abs(expected - found), tolerance64)) // &
         &  '.' ) )
    deallocate(msg)

  end subroutine testEquals_ScalarAndLocation

  subroutine testEquals_ScalarInfinity_equal()
    use MakeInfinity_mod, only:  makeInf_64, makeInf_32
    
    call assertEqual(makeInf_32(), makeInf_32(), 'equal inf 32')
    call assertEqual(makeInf_64(), makeInf_64(), 'equal inf 64')

    call assertEqual(makeInf_32(), [makeInf_32(), makeInf_32()], 'equal inf array')

  end subroutine testEquals_ScalarInfinity_equal

  subroutine testEquals_ScalarInfinity_unequal_A()
    use MakeInfinity_mod, only: makeInf_64
    
    call assertEqual(1.d0, makeInf_64(), 'unequal')
    call assertCatch( &
         & appendWithSpace('unequal', &
         & trim(valuesReport(1., makeInf_64())) // &
         & '; ' // trim(differenceReport(makeInf_64(), 0.)) // &
         &  '.' ) )

  end subroutine testEquals_ScalarInfinity_unequal_A

  subroutine testEquals_ScalarInfinity_unequal_B()
    use MakeInfinity_mod, only: makeInf_64
    
    call assertEqual(makeInf_64(), 1.0d0, 'unequal')
    call assertCatch( &
         & appendWithSpace('unequal', &
         & trim(valuesReport(makeInf_64(), 1.0d0)) // &
         & '; ' // trim(differenceReport(makeInf_64(), 0.)) // &
         &  '.' ) )

  end subroutine testEquals_ScalarInfinity_unequal_B

  subroutine testEquals_ScalarInfinity_unequal_C()
    use MakeInfinity_mod, only: makeInf_64
    
    call assertEqual(1.d0, [makeInf_64(), 1.d0], 'unequal')
    call assertCatch( &
         & appendWithSpace('unequal', &
         & trim(valuesReport(1., makeInf_64())) // &
         & '; ' // trim(differenceReport(makeInf_64(), 0.)) // &
         &  ';  first difference at element [1].' ) )


    call assertEqual(1.d0, [1.d0, makeInf_64()], 'unequal')
    call assertCatch( &
         & appendWithSpace('unequal', &
         & trim(valuesReport(1., makeInf_64())) // &
         & '; ' // trim(differenceReport(makeInf_64(), 0.)) // &
         &  ';  first difference at element [2].' ) )

  end subroutine testEquals_ScalarInfinity_unequal_C

  subroutine testEquals_MultiD_SingleElementGT1
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1
    real, parameter :: bad  = 999

    real(kind=r32), dimension(:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    !dbg1 print *,'5000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good+1; found = good;
    i1 = 1; i2 = 2; i3 = 1; i4 = 2
    found(i1,i2,i3,i4) = bad

    allocate(msg,source='testEquals_MultiD_SingleElementGT1')

    ! The following should throw an exception...
    call assertGreaterThan(expected,found, message=msg)

    ! "locationInArray" is not used in the original AssertEqual code.
    write(locationInArray,locationFormat( [i1,i2,i3,i4] )) [i1, i2, i3, i4]

    call assertCatch( &
         & appendWithSpace(msg, &
         & trim(valuesReport(good+1, bad, & 
         &      ePrefix='expected', &
         &      fPrefix='to be greater than:' )) // &
         & ';  first difference at element ' // trim(locationInArray) // '.') &
         & )

    deallocate(msg)

  end subroutine testEquals_MultiD_SingleElementGT1

  subroutine testEquals_MultiD_SingleElementGT2
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1

    real(kind=r32), dimension(:,:,:,:), allocatable :: expected, found

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray
!    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4

    !dbg1 print *,'5000'

    n1 = 2; n2 = 3; n3 = 2; n4 = 2; 
    allocate(expected(n1,n2,n3,n4),found(n1,n2,n3,n4))
    expected = good+1; found = good;
!    i1 = 1; i2 = 2; i3 = 1; i4 = 2

    allocate(msg,source='testEquals_MultiD_SingleElementGT2')

    ! The following should not throw an exception...
    call assertGreaterThan(expected,found, message=msg)

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

  end subroutine testEquals_MultiD_SingleElementGT2

  subroutine testEquals_MultiD_SingleEltVarious1
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1

    real(kind=r32), dimension(:,:,:,:), allocatable :: expected, found
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

    allocate(msg,source='testEquals_MultiD_SingleEltVarious1')

    ! The following should not throw an exception...

    expected = good; found = good;
    call assertEqual(expected, found, message=msg)

    expected = good+1; found = good;
    call assertGreaterThan(expected, found, message=msg)

    expected = good+1; found = good; found(i1,i2,i3,i4) = good+1
    call assertGreaterThanOrEqual(expected, found, message=msg)

    expected = good; found = good+1;
    call assertLessThan(expected, found, message=msg)

    expected = good; found = good+1; found(i1,i2,i3,i4) = good
    call assertLessThanOrEqual(expected, found, message=msg)

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

  end subroutine testEquals_MultiD_SingleEltVarious1

  subroutine testEquals_MultiD_SingleEltVarious2
    use Params_mod
!    use Assert_mod, only: assertEqual

    real, parameter :: good = 1

    real(kind=r32) :: expected, found
    real(kind=r32) :: tolerance32

    character(len=:), allocatable :: msg

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80
    character(len=MAXLEN_SHAPE) :: locationInArray

    !dbg1 print *,'5000'

    allocate(msg,source='testEquals_MultiD_SingleEltVarious2')

    ! The following should not throw an exception...

    expected = good; found = good;
    call assertEqual(expected, found, message=msg)

    expected = good+1; found = good;
    call assertGreaterThan(expected, found, message=msg)

    expected = good+1; found = good; found = good+1
    call assertGreaterThanOrEqual(expected, found, message=msg)

    expected = good; found = good+1;
    call assertLessThan(expected, found, message=msg)

    expected = good; found = good+1; found = good
    call assertLessThanOrEqual(expected, found, message=msg)

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

  end subroutine testEquals_MultiD_SingleEltVarious2



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

  
end module Test_AssertReal_mod

  
