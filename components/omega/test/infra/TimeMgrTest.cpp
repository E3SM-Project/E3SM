//===-- Test driver for OMEGA Time Manager------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Time Manager
///
/// This driver tests the OMEGA time manager module that tracks simulation time
/// during integrations and manages alarms to trigger events at precise moments.
/// This unit test driver consists of six functions that test each one of the
/// classes that comprise the time manager, confirming the functionality of
/// all defined constructors, accessors, methods, and operators.
//
//===-----------------------------------------------------------------------===/

#include "TimeMgr.h"
#include "DataTypes.h"
#include "Logging.h"

//------------------------------------------------------------------------------
// TimeFrac test

int testTimeFrac(void) {

   LOG_INFO("TimeMgrTest: TimeFrac tests ------------------------------------");

   // Initialize error codes
   OMEGA::I4 Err1{0};
   OMEGA::I4 Err2{0};
   OMEGA::I4 ErrAll{0};

   // Initialize some reference values for the fractional
   // representation of 2 1/3 seconds.
   OMEGA::I8 WRef{2};
   OMEGA::I8 NRef{1};
   OMEGA::I8 DRef{3};
   OMEGA::R8 RRef{2.3333333333333333};

   OMEGA::I8 WTst{2};
   OMEGA::I8 NTst{1};
   OMEGA::I8 DTst{3};
   OMEGA::R8 RTst{2.3333333333333333};

   // Test default constructor to create a reference fraction
   // Also implicitly tests one form of the get routine.

   OMEGA::TimeFrac RefTF;

   Err1 = RefTF.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == 0 && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeFrac: default constructor and get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: default constructor or get: FAIL");
   }

   // Test set/get by each component to set reference values

   Err1 = RefTF.setWhole(WRef);
   if (Err1 != 0) {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: setWhole: FAIL");
   }
   Err1 = RefTF.setNumer(NRef);
   if (Err1 != 0) {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: setNumer: FAIL");
   }
   Err1 = RefTF.setDenom(DRef);
   if (Err1 != 0) {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: setDenom: FAIL");
   }

   WTst = RefTF.getWhole();
   NTst = RefTF.getNumer();
   DTst = RefTF.getDenom();

   if (WTst == WRef) {
      LOG_INFO("TimeMgrTest/TimeFrac: setWhole/getWhole: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: setWhole/getWhole: FAIL");
   }

   if (NTst == NRef) {
      LOG_INFO("TimeMgrTest/TimeFrac: setNumer/getNumer: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: setNumer/getNumer: FAIL");
   }

   if (DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeFrac: setDenom/getDenom: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: setDenom/getDenom: FAIL");
   }

   // Test component constructor

   OMEGA::TimeFrac Tst1TF(WRef, NRef, DRef);

   Err1 = Tst1TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == WRef && NTst == NRef && DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeFrac: component constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: component constructor: FAIL");
   }

   // Can now test equivalence operator

   if (Tst1TF == RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(==): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(==): FAIL");
   }

   // Test unified set call

   Err1 = Tst1TF.set(0, 0, 1);
   WTst = Tst1TF.getWhole();
   NTst = Tst1TF.getNumer();
   DTst = Tst1TF.getDenom();

   if (Err1 == 0 && WTst == 0 && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeFrac: unified set: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: unified set: FAIL");
   }

   // Test non-equivalence

   if (Tst1TF != RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(!=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(!=): FAIL");
   }

   // Test < operator (and < part of <= operator)

   if (Tst1TF < RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(<): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(<): FAIL");
   }

   if (Tst1TF <= RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(<=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(<=): FAIL");
   }

   // Test > operator (and > part of >= operator)

   if (RefTF > Tst1TF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(>): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(>): FAIL");
   }

   if (RefTF >= Tst1TF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(>=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(>=): FAIL");
   }

   // Test assignment operator and = part of above comparisons

   Tst1TF = RefTF;

   if (Tst1TF == RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: assignment operator: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: assignment operator: FAIL");
   }

   if (Tst1TF <= RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(<=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(<=): FAIL");
   }

   if (RefTF >= Tst1TF) {
      LOG_INFO("TimeMgrTest/TimeFrac: operator(>=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: operator(>=): FAIL");
   }

   // Test copy constuctor

   OMEGA::TimeFrac Tst2TF(RefTF);

   if (Tst2TF == RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: copy constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: copy constructor: FAIL");
   }

   // Test addition

   Tst1TF.set(2, 2, 3);
   Tst2TF.set(1, 1, 5);

   OMEGA::TimeFrac Tst3TF;
   Tst3TF = Tst1TF + Tst2TF;

   Err1 = Tst3TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == 3 && NTst == 13 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: addition: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: addition: FAIL");
   }

   // Test increment

   Tst1TF.set(2, 2, 3);
   Tst2TF.set(1, 1, 5);

   Tst3TF = Tst1TF;
   Tst3TF += Tst2TF;

   Err1 = Tst3TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == 3 && NTst == 13 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: increment: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: increment: FAIL");
   }

   // Test subtraction

   Tst1TF.set(2, 2, 3);
   Tst2TF.set(1, 1, 5);

   Tst3TF = Tst1TF - Tst2TF;

   Err1 = Tst3TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == 1 && NTst == 7 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: subtraction: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: subtraction: FAIL");
   }

   // Test decrement

   Tst1TF.set(2, 2, 3);
   Tst2TF.set(1, 1, 5);

   Tst3TF = Tst1TF;
   Tst3TF -= Tst2TF;

   Err1 = Tst3TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == 1 && NTst == 7 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: decrement: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: decrement: FAIL");
   }

   // Test multiply by int functions

   Tst1TF.set(2, 2, 3);
   Tst3TF         = Tst1TF;
   OMEGA::I4 ITst = 5;

   Tst2TF = Tst1TF * ITst;
   Tst3TF *= ITst;

   Err1 = Tst2TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 13 && NTst == 1 && DTst == 3) {
      LOG_INFO("TimeMgrTest/TimeFrac: int multiply: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: int multiply: FAIL");
   }

   Err1 = Tst3TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 13 && NTst == 1 && DTst == 3) {
      LOG_INFO("TimeMgrTest/TimeFrac: int multiply in place: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: int multiply in place: FAIL");
   }

   // Test multiply by real functions

   Tst1TF.set(2, 2, 3);
   Tst3TF = Tst1TF;
   RTst   = 7.55;

   Tst2TF = Tst1TF * RTst;
   Tst3TF *= RTst;

   Err1 = Tst2TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 20 && NTst == 2 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: real multiply: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: real multiply: FAIL");
   }

   Err1 = Tst3TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 20 && NTst == 2 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: real multiply in place: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: real multiply in place: FAIL");
   }

   // Test divide by int functions

   Tst1TF.set(2, 2, 3);
   Tst3TF = Tst1TF;
   WTst   = 5;

   Tst2TF = Tst1TF / WTst;
   Tst3TF /= WTst;

   Err1 = Tst2TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 0 && NTst == 8 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: divide: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: divide: FAIL");
   }

   Err1 = Tst3TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 0 && NTst == 8 && DTst == 15) {
      LOG_INFO("TimeMgrTest/TimeFrac: divide in place: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: divide in place: FAIL");
   }

   // Test divide fractions

   Tst1TF.set(2, 2, 3);
   Tst2TF.set(3, 1, 5);

   RTst = Tst1TF / Tst2TF;

   if (fabs(RTst - 0.8333333333333333) < 1.e-15) {
      LOG_INFO("TimeMgrTest/TimeFrac: divide fractions: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: divide fractions: FAIL");
   }

   // Test modulo functions

   Tst1TF.set(3, 4, 15);
   Tst2TF.set(5, 7, 10);

   Tst3TF = Tst1TF % Tst2TF;
   Tst1TF %= Tst2TF;

   Err1 = Tst3TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 0 && NTst == 98 && DTst == 171) {
      LOG_INFO("TimeMgrTest/TimeFrac: modulo: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: modulo: FAIL");
   }

   Err1 = Tst1TF.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 0 && NTst == 98 && DTst == 171) {
      LOG_INFO("TimeMgrTest/TimeFrac: modulo in place: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: modulo in place: FAIL");
   }

   // Test get/set for integer hour, minute, second interfaces

   OMEGA::I4 HRef{13};
   OMEGA::I4 MRef{87};
   OMEGA::I4 SRef{36};
   OMEGA::I4 HTst{0};
   OMEGA::I4 MTst{0};
   OMEGA::I4 STst{0};

   Err1 = Tst1TF.setHMS(HRef, MRef, SRef);
   Err2 = Tst1TF.getHMS(HTst, MTst, STst);

   if (Err1 == 0 && Err2 == 0 && HTst == 14 && MTst == 27 && STst == 36) {
      LOG_INFO("TimeMgrTest/TimeFrac: getHMS/setHMS: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: getHMS/setHMS: FAIL");
   }

   // Test real hour minute second interfaces
   // First test real seconds.

   RRef = 7.8;
   RTst = 0.0;

   Err1 = Tst1TF.setSeconds(RRef);
   RTst = Tst1TF.getSeconds();
   Err2 = Tst1TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && Err2 == 0 && fabs(RTst - RRef) < 1.e-15 && WTst == 7 &&
       NTst == 4 && DTst == 5) {
      LOG_INFO("TimeMgrTest/TimeFrac: get/set real seconds: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: get/set real seconds: FAIL");
   }

   // Test the related constructor from real seconds.

   OMEGA::TimeFrac Tst4TF(RRef);

   if (Tst4TF == Tst1TF) {
      LOG_INFO("TimeMgrTest/TimeFrac: real seconds constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: real seconds constructor: FAIL");
   }

   // Test real hours.

   RRef = 3.55;
   RTst = 0.0;

   Err1 = Tst1TF.setHours(RRef);
   RTst = Tst1TF.getHours();
   Err2 = Tst1TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && Err2 == 0 && fabs(RTst - RRef) < 1.e-15 && WTst == 12780 &&
       NTst == 0 & DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeFrac: get/set real hours: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: get/set real hours: FAIL");
   }

   // Test real minutes.

   RRef = 5.0875;
   RTst = 0.0;

   Err1 = Tst1TF.setMinutes(RRef);
   RTst = Tst1TF.getMinutes();
   Err2 = Tst1TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && Err2 == 0 && fabs(RTst - RRef) < 1.e-15 && WTst == 305 &&
       NTst == 1 & DTst == 4) {
      LOG_INFO("TimeMgrTest/TimeFrac: get/set real minutes: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: get/set real minutes: FAIL");
   }

   // Test simplify function (not exhaustive test)

   Tst2TF.set(2, 5, 3);
   Err2 = Tst2TF.simplify();
   Err1 = Tst2TF.get(WTst, NTst, DTst);

   if (Err1 == 0 && DTst == 3 && WTst == 3 && NTst == 2) {
      LOG_INFO("TimeMgrTest/TimeFrac: simplify: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: simplify: FAIL");
   }

   // Test convert function

   Tst2TF.set(2, 5, 3);
   DTst = 15;
   Err1 = Tst2TF.convert(DTst);
   Err2 = Tst2TF.get(WTst, NTst, DTst);

   // note that convert leaves an improper fraction - no simplify
   if (Err1 == 0 && Err2 == 0 && DTst == 15 && WTst == 0 && NTst == 55) {
      LOG_INFO("TimeMgrTest/TimeFrac: convert: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: convert: FAIL");
   }

   return ErrAll;

} // end testTimeFrac

//------------------------------------------------------------------------------
// The test driver.

int main(int argc, char *argv[]) {

   OMEGA::I4 Err{0};
   OMEGA::I4 TotErr{0};

   Err = testTimeFrac();
   TotErr += Err;

   if (TotErr == 0) {
      LOG_INFO("TimeMgrTest: Successful completion");
   } else {
      LOG_INFO("TimeMgrTest: Failed");
   }

} // end of main
//===-----------------------------------------------------------------------===/

