!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_AssertBasic_mod
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
module Test_AssertBasic_mod
   use Exception_mod, only: NULL_MESSAGE
   use AssertBasic_mod
   use TestSuite_mod, only: TestSuite, newTestSuite
   use Exception_mod, only: catch
   use Exception_mod, only: getNumExceptions
   implicit none
   private

   public :: suite

contains

   function suite()
      use TestSuite_mod, only: TestSuite, newTestSuite
      use TestMethod_mod, only: newTestMethod
      use Test_mod

      type (TestSuite) :: suite

      suite = newTestSuite('AssertIntegerTests')
!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

      call suite%addTest( &
           &   newTestMethod('testAssertTrueF', &
           &                  testAssertTrueF))
      call suite%addTest( &
           &   newTestMethod('testAssertTrueF1', &
           &                  testAssertTrueF1))
      call suite%addTest( &
           &   newTestMethod('testAssertTrueF2', &
           &                  testAssertTrueF2))
      call suite%addTest( &
           &   newTestMethod('testAssertTrueT', &
           &                  testAssertTrueT))
      call suite%addTest( &
           &   newTestMethod('testAssertTrueT1', &
           &                  testAssertTrueT1))
      call suite%addTest( &
           &   newTestMethod('testAssertFalseT', &
           &                  testAssertFalseT))
      call suite%addTest( &
           &   newTestMethod('testAssertFalseT1', &
           &                  testAssertFalseT1))
      call suite%addTest( &
           &   newTestMethod('testAssertFalseF', &
           &                  testAssertFalseF))
      call suite%addTest( &
           &   newTestMethod('testAssertFalseF1', &
           &                  testAssertFalseF1))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringSame', &
           &                  testAssertEqualStringSame))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringDifferent', &
           &                  testAssertEqualStringDifferent))

      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnAllWhite1', &
           &                  testAssertEqualStrIgnAllWhite1))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnAllWhite2', &
           &                  testAssertEqualStrIgnAllWhite2))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnAllWhite3', &
           &                  testAssertEqualStrIgnAllWhite3))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnAllWhite4', &
           &                  testAssertEqualStrIgnAllWhite4))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnAllWhite5', &
           &                  testAssertEqualStrIgnAllWhite5))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnAllWhite6', &
           &                  testAssertEqualStrIgnAllWhite6))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnAllWhite7', &
           &                  testAssertEqualStrIgnAllWhite7))

      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif1', &
           &                  testAssertEqualStrIgnWhiDif1))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif2', &
           &                  testAssertEqualStrIgnWhiDif2))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif3', &
           &                  testAssertEqualStrIgnWhiDif3))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif4', &
           &                  testAssertEqualStrIgnWhiDif4))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif5', &
           &                  testAssertEqualStrIgnWhiDif5))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif6', &
           &                  testAssertEqualStrIgnWhiDif6))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif7', &
           &                  testAssertEqualStrIgnWhiDif7))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif8', &
           &                  testAssertEqualStrIgnWhiDif8))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStrIgnWhiDif9', &
           &                  testAssertEqualStrIgnWhiDif9))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringTrimWhitespace1', &
           &                  testAssertEqualStringTrimWhitespace1))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringTrimWhitespace2', &
           &                  testAssertEqualStringTrimWhitespace2))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringTrimWhitespace3', &
           &                  testAssertEqualStringTrimWhitespace3))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringTrimWhitespace4', &
           &                  testAssertEqualStringTrimWhitespace4))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringKeepWhitespace1', &
           &                  testAssertEqualStringKeepWhitespace1))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringKeepWhitespace2', &
           &                  testAssertEqualStringKeepWhitespace2))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringKeepWhitespace3', &
           &                  testAssertEqualStringKeepWhitespace3))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualStringKeepWhitespace4', &
           &                  testAssertEqualStringKeepWhitespace4))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualNonzeroBlanks1', &
           &                  testAssertEqualNonzeroBlanks1 ))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualNonzeroBlanks2', &
           &                  testAssertEqualNonzeroBlanks2 ))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualNonzeroBlanks3', &
           &                  testAssertEqualNonzeroBlanks3 ))
      call suite%addTest( &
           &   newTestMethod('testAssertEqualNonzeroBlanks4', &
           &                  testAssertEqualNonzeroBlanks4 ))
      call suite%addTest( &
           &   newTestMethod('testAssertAny', &
           &                  testAssertAny))
      call suite%addTest( &
           &   newTestMethod('testAssertAnyFail', &
           &                  testAssertAnyFail))
      call suite%addTest( &
           &   newTestMethod('testAssertAll', &
           &                  testAssertAll))
      call suite%addTest( &
           &   newTestMethod('testAssertAllFail', &
           &                  testAssertAllFail))
      call suite%addTest( &
           &   newTestMethod('testAssertNone', &
           &                  testAssertNone))
      call suite%addTest( &
           &   newTestMethod('testAssertNoneFail', &
           &                  testAssertNoneFail))
      call suite%addTest( &
           &   newTestMethod('testAssertNotAll', &
           &                  testAssertNotAll))
      call suite%addTest( &
           &   newTestMethod('testAssertNotAllFail', &
           &                  testAssertNotAllFail))

      call suite%addTest( &
           &   newTestMethod('testAssertIsNaN', &
           &                  testAssertIsNaN))
      call suite%addTest( &
           &   newTestMethod('testAssertIsFinite', &
           &                  testAssertIsFinite))

      call suite%addTest( &
           &   newTestMethod('testAssertFail', &
           &                  testAssertFail))
      call suite%addTest( &
           &   newTestMethod('testAssertExceptionRaised', &
           &                  testAssertExceptionRaised))
   end function suite

   subroutine testAssertTrueF()
      call assertTrue(.false.)
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertTrueF

   subroutine testAssertTrueF1()
      call assertTrue([.false.].eqv.[.true.])
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertTrueF1

   subroutine testAssertTrueF2()
      call assertTrue([.true.,.false.].eqv.[.true.,.true.])
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertTrueF2

   subroutine testAssertTrueT()
      call assertTrue(.true.)
   end subroutine testAssertTrueT

   subroutine testAssertTrueT1()
      call assertTrue([.true.,.true.].eqv.[.true.,.true.])
   end subroutine testAssertTrueT1

   subroutine testAssertFalseF()
      call assertFalse(.false.)
   end subroutine testAssertFalseF

   subroutine testAssertFalseF1()
      call assertFalse([.false.,.false.].eqv.[.true.,.true.])
   end subroutine testAssertFalseF1

   subroutine testAssertFalseT()
      call assertFalse(.true.)
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertFalseT

   subroutine testAssertFalseT1()
      call assertFalse([.true.,.true.,.true.])
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertFalseT1

   subroutine testAssertEqualStringSame()
      call assertEqual(expected="string A", found="string A")
   end subroutine testAssertEqualStringSame

   subroutine testAssertEqualStringDifferent()
      call assertEqual(expected="string A", found="string B")
      call assertTrue(catch('String assertion failed:' // new_line('A') // &
           & '    expected: <"string A">' // new_line('A') // &
           & '   but found: <"string B">' // new_line('A') // &
           & '  first diff:   -------^'),'Unexpected Equal String')
   end subroutine testAssertEqualStringDifferent

   subroutine testAssertEqualStrIgnAllWhite1()
     call assertEqual(expected="stringA", found="string A", &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualStrIgnAllWhite1

   subroutine testAssertEqualStrIgnAllWhite2()
     call assertEqual(expected="string A", found="stringA", &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualStrIgnAllWhite2

   subroutine testAssertEqualStrIgnAllWhite3()
     call assertEqual(expected="stringA ", found="string A", &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualStrIgnAllWhite3

   subroutine testAssertEqualStrIgnAllWhite4()
     call assertEqual(expected=" string A", found="string A", &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualStrIgnAllWhite4

   subroutine testAssertEqualStrIgnAllWhite5()
     call assertEqual(expected=" string A ", found="stringA", &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualStrIgnAllWhite5

   subroutine testAssertEqualStrIgnAllWhite6()
     character tab; tab = char(9)
     call assertEqual(expected=tab//"string A ", found="stringA", &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualStrIgnAllWhite6

   subroutine testAssertEqualStrIgnAllWhite7()
     character tab; tab = char(9)
     call assertEqual(expected=tab//"string A ", found=" stringA"//tab, &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualStrIgnAllWhite7

   subroutine testAssertEqualStrIgnWhiDif1()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected=tab, found=spc, &
          & whitespace=IGNORE_DIFFERENCES)
   end subroutine

   subroutine testAssertEqualStrIgnWhiDif2()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected=spc, found=tab, &
          & whitespace=IGNORE_DIFFERENCES)
   end subroutine

   subroutine testAssertEqualStrIgnWhiDif3()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected=tab//"string A", found="string A", &
          & whitespace=IGNORE_DIFFERENCES)
   end subroutine

   subroutine testAssertEqualStrIgnWhiDif4()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected=tab//"string A", found=spc//"string A", &
          & whitespace=IGNORE_DIFFERENCES)
   end subroutine

   subroutine testAssertEqualStrIgnWhiDif5()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected=tab//"string A", found=spc//"string A"//spc, &
          & whitespace=IGNORE_DIFFERENCES)
   end subroutine

   subroutine testAssertEqualStrIgnWhiDif6()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual( &
          & expected =spc//"string"//tab//spc//"6", &
          & found    ="string"//spc//spc//spc//spc//"6"//spc, &
          & whitespace=IGNORE_DIFFERENCES)
   end subroutine

   subroutine testAssertEqualStrIgnWhiDif7()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected=tab//"string 7 ", found="string7", &
          & whitespace=IGNORE_DIFFERENCES)

     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <"'//tab//'string 7">' // new_line('A') // &
          & '   but found: <"string7">'  // new_line('A') // &
          & '  first diff:   ------^'))

   end subroutine

   subroutine testAssertEqualStrIgnWhiDif8()
     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected="string8A", found="string8B", &
          & whitespace=IGNORE_DIFFERENCES)

     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <"string8A">' // new_line('A') // &
          & '   but found: <"string8B">'  // new_line('A') // &
          & '  first diff:   -------^'))

   end subroutine

   subroutine testAssertEqualStrIgnWhiDif9()
!     character, parameter :: tab = char(9), spc = char(32)
     call assertEqual(expected="", found=" ", &
          & whitespace=IGNORE_DIFFERENCES)

!     call assertTrue(catch( &
!          & 'String assertion failed:' // new_line('A') // &
!          & '    expected: <"string8A">' // new_line('A') // &
!          & '   but found: <"string8B">'  // new_line('A') // &
!          & '  first diff:   -------^'))

   end subroutine

   subroutine testAssertEqualStringTrimWhitespace1()
     character tab; tab = char(9)
     call assertEqual(expected=tab//"string A ", found="string A", &
          & whitespace=TRIM_ALL )
   end subroutine testAssertEqualStringTrimWhitespace1

   subroutine testAssertEqualStringTrimWhitespace2()
     character tab; tab = char(9)
     ! Should fail !
     call assertEqual(expected=tab//"string A ", found="stringA", &
          & whitespace=TRIM_ALL )
     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <"'//tab//'string A">' // new_line('A') // &
          & '   but found: <"stringA">'  // new_line('A') // &
          & '  first diff:   ------^'))
   end subroutine testAssertEqualStringTrimWhitespace2
   
   subroutine testAssertEqualStringTrimWhitespace3()
     character tab; tab = char(9)
     ! Should fail !
     call assertEqual(expected="", found=" ", &
          & whitespace=TRIM_ALL )
   end subroutine testAssertEqualStringTrimWhitespace3

   subroutine testAssertEqualStringTrimWhitespace4()
     call assertEqual( &
          & expected = "i= 1 f= F s=word x=  1.23", &
          & found    = "i= 1  f= F s=word x=  1.23", &
!          & expected = "i= 1 f= F s=word x=  1.23", &
!          & found    = “i= 1  f= F s=word x=  1.23", &
          & whitespace=TRIM_ALL )
!          & whitespace=IGNORE_DIFFERENCES )
     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <"i= 1 f= F s=word x=  1.23">' // new_line('A') // &
          & '   but found: <"i= 1  f= F s=word x=  1.23">' // new_line('A') // &
          & '  first diff:   -----^'))
   end subroutine

   subroutine testAssertEqualStringKeepWhitespace1()
     character tab; tab = char(9)
     call assertEqual(expected=tab//"string A ", found=tab//"string A ", &
          & whitespace=KEEP_ALL )
   end subroutine testAssertEqualStringKeepWhitespace1

   subroutine testAssertEqualStringKeepWhitespace2()
     character tab; tab = char(9)
     ! A strict interpretation of keep:  TAB != SPC !
     call assertEqual(expected=tab//"string A ", found=" string A ", &
          & whitespace=KEEP_ALL )
     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <"	string A ">' // new_line('A') // &
          & '   but found: <" string A ">' // new_line('A') // &
          & '  first diff:   ^'))
   end subroutine testAssertEqualStringKeepWhitespace2

   subroutine testAssertEqualStringKeepWhitespace3()
     character tab; tab = char(9)
     ! Should fail !
     call assertEqual(expected="", found=" ", &
          & whitespace=KEEP_ALL )
     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <"">' // new_line('A') // &
          & '   but found: <" ">'  // new_line('A') // &
          & '  first diff:   ^'))
   end subroutine testAssertEqualStringKeepWhitespace3

   subroutine testAssertEqualStringKeepWhitespace4()
     character tab; tab = char(9)
     ! Should fail !
     call assertEqual(expected=" ", found="", &
          & whitespace=KEEP_ALL )
     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <" ">' // new_line('A') // &
          & '   but found: <"">'  // new_line('A') // &
          & '  first diff:   ^'))
   end subroutine testAssertEqualStringKeepWhitespace4

   subroutine testAssertEqualNonzeroBlanks1
     call assertEqual(expected=" ", found="      ", &
          & whitespace=KEEP_ALL )
     call assertTrue(catch( &
          & 'String assertion failed:' // new_line('A') // &
          & '    expected: <" ">' // new_line('A') // &
          & '   but found: <"      ">'  // new_line('A') // &
          & '  first diff:   -^'))
   end subroutine testAssertEqualNonzeroBlanks1

   subroutine testAssertEqualNonzeroBlanks2
     call assertEqual(expected=" ", found="      ", &
          & whitespace=TRIM_ALL )
   end subroutine testAssertEqualNonzeroBlanks2

   subroutine testAssertEqualNonzeroBlanks3
     call assertEqual(expected=" ", found="      ", &
          & whitespace=IGNORE_DIFFERENCES )
   end subroutine testAssertEqualNonzeroBlanks3

   subroutine testAssertEqualNonzeroBlanks4
     call assertEqual(expected=" ", found="      ", &
          & whitespace=IGNORE_ALL )
   end subroutine testAssertEqualNonzeroBlanks4
   
   ! Fail only if all .false.
   subroutine testAssertAny()
      call assertAny([.true.])
      call assertAny([.true., .true.])
      call assertAny([.true.,.false.])
      call assertAny([.false.,.true.])
   end subroutine testAssertAny

   subroutine testAssertAnyFail()
      call assertAny([.false.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertAny([.false.,.false.])
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertAnyFail

   ! Fail if any .false.
   subroutine testAssertAll()
      call assertAll([.true.])
      call assertAll([.true., .true.])
   end subroutine testAssertAll

   subroutine testAssertAllFail()
      call assertAll([.false.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertAll([.false.,.false.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertAll([.true.,.false.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertAll([.false.,.true.])
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertAllFail

   ! Fail if any .true.
   subroutine testAssertNone()
      call assertNone([.false.])
      call assertNone([.false., .false.])
   end subroutine testAssertNone

   subroutine testAssertNoneFail()
      call assertNone([.true.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertNone([.false.,.true.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertNone([.true.,.false.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertNone([.true.,.true.])
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertNoneFail


   ! Fail if any .true.
   subroutine testAssertNotAll()
      call assertNotAll([.false.])
      call assertNotAll([.false., .true.])
      call assertNotAll([.true., .false.])
      call assertNotAll([.false., .false.])
   end subroutine testAssertNotAll

   subroutine testAssertNotAllFail()
      call assertNotAll([.true.])
      call assertTrue(catch(NULL_MESSAGE))
      call assertNotAll([.true.,.true.])
      call assertTrue(catch(NULL_MESSAGE))
   end subroutine testAssertNotAllFail

   subroutine testAssertIsNaN()
      use MakeNaN_mod, only: makeNaN_32, makeNaN_64

      call assertIsNaN(1.e0, 'not NaN')
      call assertExceptionRaised('not NaN')
      call assertIsNaN(1.d0, 'not NaN')
      call assertExceptionRaised('not NaN')

      call assertIsNaN(makeNaN_32())
      call assertIsNaN(makeNaN_64())
 
   end subroutine testAssertIsNaN

   subroutine testAssertIsFinite()
      use MakeInfinity_mod, only: makeInf_32, makeInf_64

      call assertIsFinite(1.e0, 'finite')
      call assertIsFinite(1.d0, 'finite')

      call assertIsFinite(makeInf_32(), 'not finite')
      call assertExceptionRaised('not finite')
      call assertIsFinite(makeInf_64(), 'not finite')
      call assertExceptionRaised('not finite')

   end subroutine testAssertIsFinite

   subroutine testAssertExceptionRaised()
      use Exception_mod, only: throw
      use SourceLocation_mod

      character(len=*), parameter :: message = 'a message'

      call throw(message)
      call assertExceptionRaised(message)

      call throw(message)
      call assertExceptionRaised(message,SourceLocation('here',5))

   end subroutine testAssertExceptionRaised

   subroutine testAssertFail()
      use SourceLocation_mod

      character(len=*), parameter :: message = 'a message'

      call assertFail(message)
      call assertExceptionRaised(message)

      call assertFail(message)
      call assertExceptionRaised(message,SourceLocation('here',5))
   end subroutine testAssertFail

end module Test_AssertBasic_mod
