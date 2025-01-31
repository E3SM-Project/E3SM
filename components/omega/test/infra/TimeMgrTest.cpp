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
#include "MachEnv.h"
#include "mpi.h"

using namespace OMEGA;

//------------------------------------------------------------------------------
// TimeFrac test

int testTimeFrac(void) {

   LOG_INFO("TimeMgrTest: TimeFrac tests ------------------------------------");

   // Initialize error codes
   I4 Err1{0};
   I4 Err2{0};
   I4 ErrAll{0};

   // Initialize some reference values for the fractional
   // representation of 2 1/3 seconds.
   I8 WRef{2};
   I8 NRef{1};
   I8 DRef{3};
   R8 RRef{2.3333333333333333};

   I8 WTst{2};
   I8 NTst{1};
   I8 DTst{3};
   R8 RTst{2.3333333333333333};

   // Test default constructor to create a reference fraction
   // Also implicitly tests one form of the get routine.

   TimeFrac RefTF;

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

   TimeFrac Tst1TF(WRef, NRef, DRef);

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

   TimeFrac Tst2TF(RefTF);

   if (Tst2TF == RefTF) {
      LOG_INFO("TimeMgrTest/TimeFrac: copy constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeFrac: copy constructor: FAIL");
   }

   // Test addition

   Tst1TF.set(2, 2, 3);
   Tst2TF.set(1, 1, 5);

   TimeFrac Tst3TF;
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
   Tst3TF  = Tst1TF;
   I4 ITst = 5;

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

   I4 HRef{13};
   I4 MRef{87};
   I4 SRef{36};
   I4 HTst{0};
   I4 MTst{0};
   I4 STst{0};

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

   TimeFrac Tst4TF(RRef);

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
// Calendar test

int testCalendar(void) {

   LOG_INFO("TimeMgrTest: Calendar tests ------------------------------------");

   // Initialize error codes
   I4 Err1{0};
   I4 Err2{0};
   I4 ErrAll{0};

   // Test custom calendar
   // Also tests the get routine.

   CalendarKind Kind0 = CalendarCustom;
   I4 MonthsPerYear0  = 12;
   std::vector<I4> DaysPerMonth0(MonthsPerYear0, 0);
   DaysPerMonth0[0]   = 10;
   DaysPerMonth0[1]   = 10;
   DaysPerMonth0[2]   = 10;
   DaysPerMonth0[3]   = 10;
   DaysPerMonth0[4]   = 10;
   DaysPerMonth0[5]   = 10;
   DaysPerMonth0[6]   = 10;
   DaysPerMonth0[7]   = 10;
   DaysPerMonth0[8]   = 10;
   DaysPerMonth0[9]   = 10;
   DaysPerMonth0[10]  = 10;
   DaysPerMonth0[11]  = 14;
   I4 SecondsPerDay0  = 100;
   I4 SecondsPerYear0 = 12400;
   I4 DaysPerYear0    = 124;

   CalendarKind Kind1 = CalendarCustom;
   I4 MonthsPerYear1  = 12;
   std::vector<I4> DaysPerMonth1(MonthsPerYear1, 0);
   DaysPerMonth1[0]   = 99;
   DaysPerMonth1[1]   = 99;
   DaysPerMonth1[2]   = 99;
   DaysPerMonth1[3]   = 99;
   DaysPerMonth1[4]   = 99;
   DaysPerMonth1[5]   = 99;
   DaysPerMonth1[6]   = 99;
   DaysPerMonth1[7]   = 99;
   DaysPerMonth1[8]   = 99;
   DaysPerMonth1[9]   = 99;
   DaysPerMonth1[10]  = 99;
   DaysPerMonth1[11]  = 99;
   I4 SecondsPerDay1  = 999;
   I4 SecondsPerYear1 = 999;
   I4 DaysPerYear1    = 999;

   Calendar::init(DaysPerMonth0, SecondsPerDay0, SecondsPerYear0, DaysPerYear0);

   Calendar *CalCustom = Calendar::get();
   Kind1               = Calendar::getKind();
   DaysPerMonth1       = Calendar::getDaysPerMonth();
   MonthsPerYear1      = Calendar::getMonthsPerYear();
   SecondsPerDay1      = Calendar::getSecondsPerDay();
   SecondsPerYear1     = Calendar::getSecondsPerYear();
   DaysPerYear1        = Calendar::getDaysPerYear();

   if (Err1 != 0 || Kind1 != Kind0 || MonthsPerYear1 != MonthsPerYear0 ||
       SecondsPerDay1 != SecondsPerDay0 || SecondsPerYear1 != SecondsPerYear0 ||
       DaysPerYear1 != DaysPerYear0) {
      Err2 = 1;
   } else {
      Err2 = 0;
   }

   for (int I = 0; I < MonthsPerYear0; I++) {
      if (DaysPerMonth0[I] != DaysPerMonth1[I])
         Err2 = 1;
   }

   if (Err2 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: custom constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: custom constructor: FAIL");
   }

   // Test a date in custom calendar (1957-10-4, 00:01:24.25)
   I8 TmpSeconds = (1957 * (I8)12400) + (9 * 10 + 4 - 1) * 100 + 60 + 24;
   TimeFrac TstTime(TmpSeconds, 1, 4);
   I8 TstYear    = 1957;
   I8 TstMonth   = 10;
   I8 TstDay     = 4;
   I8 TstHour    = 0;
   I8 TstMinute  = 1;
   I8 TstSecondW = 24;
   I8 TstSecondN = 1;
   I8 TstSecondD = 4;

   TimeFrac ChkTime(0, 0, 1);
   I8 ChkYear    = 0;
   I8 ChkMonth   = 0;
   I8 ChkDay     = 0;
   I8 ChkHour    = 0;
   I8 ChkMinute  = 0;
   I8 ChkSecondW = 0;
   I8 ChkSecondN = 0;
   I8 ChkSecondD = 1;

   ChkTime =
       Calendar::getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert custom date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert custom date to "
                "elapsed time: FAIL");
   }

   Err1 = Calendar::getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
                                ChkMinute, ChkSecondW, ChkSecondN, ChkSecondD);
   if (Err1 == 0 && ChkYear == TstYear && ChkMonth == TstMonth &&
       ChkDay == TstDay && ChkHour == TstHour && ChkMinute == TstMinute &&
       ChkSecondW == TstSecondW && ChkSecondN == TstSecondN &&
       ChkSecondD == TstSecondD) {
      LOG_INFO("TimeMgrTest/Calendar: convert elapsed time to "
               "Custom date: PASS");

   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert elapsed time to "
                "Custom date: FAIL");
   }

   // Test validate

   Err1 = CalCustom->validate();

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: validate: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: validate: FAIL");
   }

   //------------------------
   // Test Gregorian calendar

   Calendar::reset(); // reset calendar
   Kind0             = CalendarGregorian;
   DaysPerMonth0[0]  = 31;
   DaysPerMonth0[1]  = 28;
   DaysPerMonth0[2]  = 31;
   DaysPerMonth0[3]  = 30;
   DaysPerMonth0[4]  = 31;
   DaysPerMonth0[5]  = 30;
   DaysPerMonth0[6]  = 31;
   DaysPerMonth0[7]  = 31;
   DaysPerMonth0[8]  = 30;
   DaysPerMonth0[9]  = 31;
   DaysPerMonth0[10] = 30;
   DaysPerMonth0[11] = 31;
   MonthsPerYear0    = 12;
   SecondsPerDay0    = 86400;
   SecondsPerYear0   = 31536000;
   DaysPerYear0      = 365;

   Kind1             = CalendarNoCalendar;
   DaysPerMonth1[0]  = 99;
   DaysPerMonth1[1]  = 99;
   DaysPerMonth1[2]  = 99;
   DaysPerMonth1[3]  = 99;
   DaysPerMonth1[4]  = 99;
   DaysPerMonth1[5]  = 99;
   DaysPerMonth1[6]  = 99;
   DaysPerMonth1[7]  = 99;
   DaysPerMonth1[8]  = 99;
   DaysPerMonth1[9]  = 99;
   DaysPerMonth1[10] = 99;
   DaysPerMonth1[11] = 99;
   MonthsPerYear1    = 999;
   SecondsPerDay1    = 999;
   SecondsPerYear1   = 999;
   DaysPerYear1      = 999;

   Calendar::init("Gregorian");
   Calendar *CalGreg = Calendar::get();

   Kind1           = Calendar::getKind();
   DaysPerMonth1   = Calendar::getDaysPerMonth();
   MonthsPerYear1  = Calendar::getMonthsPerYear();
   SecondsPerDay1  = Calendar::getSecondsPerDay();
   SecondsPerYear1 = Calendar::getSecondsPerYear();
   DaysPerYear1    = Calendar::getDaysPerYear();

   if (Kind1 != Kind0 || MonthsPerYear1 != MonthsPerYear0 ||
       SecondsPerDay1 != SecondsPerDay0 || SecondsPerYear1 != SecondsPerYear0 ||
       DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Gregorian constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Gregorian constructor: FAIL");
   }

   // Test leap year with a non-leap year
   TstYear = 1981;
   if (!Calendar::isLeapYear(TstYear)) {
      if (Err1 == 0) {
         LOG_INFO("TimeMgrTest/Calendar: non-leap year Gregorian: PASS");
      } else {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Calendar: non-leap year Gregorian: FAIL");
      }
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: non-leap year Gregorian: FAIL");
   }
   // Test leap year with a leap year
   TstYear = 1984;
   if (Calendar::isLeapYear(TstYear)) {
      if (Err1 == 0) {
         LOG_INFO("TimeMgrTest/Calendar: 1984 leap year Gregorian: PASS");
      } else {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Calendar: 1984 leap year Gregorian: FAIL");
      }
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: 1984 leap year Gregorian: FAIL");
   }
   // Test special Gregorian leap year exceptions
   TstYear = 1900;
   if (!Calendar::isLeapYear(TstYear)) {
      if (Err1 == 0) {
         LOG_INFO("TimeMgrTest/Calendar: leap year exception "
                  "100 Gregorian: PASS");
      } else {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Calendar: leap year exception "
                   "100 Gregorian: FAIL");
      }
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: leap year exception "
                "100 Gregorian: FAIL");
   }
   TstYear = 2000;
   if (Calendar::isLeapYear(TstYear)) {
      if (Err1 == 0) {
         LOG_INFO("TimeMgrTest/Calendar: leap year exception "
                  "400 Gregorian: PASS");
      } else {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Calendar: leap year exception "
                   "400 Gregorian: FAIL");
      }
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: leap year exception "
                "400 Gregorian: FAIL");
   }

   // Test calendar date to/from elapsed time conversion
   // This is not an exhaustive test - just testing a single time
   // and for internal consistency

   // For Gregorian, check that Oct 4.81 1957 converts to the
   // Julian Day of 2436116.31 * SECONDS_PER_DAY for elapsed time
   Err1       = TstTime.set(210480449184, 1, 4);
   TstYear    = 1957;
   TstMonth   = 10;
   TstDay     = 4;
   TstHour    = 19;
   TstMinute  = 26;
   TstSecondW = 24;
   TstSecondN = 1;
   TstSecondD = 4;

   Err1       = ChkTime.set(0, 0, 1);
   ChkYear    = 0;
   ChkMonth   = 0;
   ChkDay     = 0;
   ChkHour    = 0;
   ChkMinute  = 0;
   ChkSecondW = 0;
   ChkSecondN = 0;
   ChkSecondD = 1;

   ChkTime =
       Calendar::getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);

   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert Gregorian date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert Gregorian date to "
                "elapsed time: FAIL");
   }

   Err1 = Calendar::getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
                                ChkMinute, ChkSecondW, ChkSecondN, ChkSecondD);
   if (Err1 == 0 && ChkYear == TstYear && ChkMonth == TstMonth &&
       ChkDay == TstDay && ChkHour == TstHour && ChkMinute == TstMinute &&
       ChkSecondW == TstSecondW && ChkSecondN == TstSecondN &&
       ChkSecondD == TstSecondD) {
      LOG_INFO("TimeMgrTest/Calendar: convert elapsed time to "
               "Gregorian date: PASS");

   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert elapsed time to "
                "Gregorian date: FAIL");
   }

   // Test calendar date increment function
   // Check normal year increment/decrement

   TstYear  = 1983;
   TstMonth = 6;
   TstDay   = 15;
   ChkYear  = 1985;
   ChkMonth = 6;
   ChkDay   = 15;

   Err1 =
       CalGreg->incrementDate(2, TimeUnits::Years, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by year: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date by year: FAIL");
   }

   // Check normal month increment/decrement
   TstYear  = 1984;
   TstMonth = 6;
   TstDay   = 15;
   ChkYear  = 1984;
   ChkMonth = 8;
   ChkDay   = 15;

   Err1 =
       CalGreg->incrementDate(2, TimeUnits::Months, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by month: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date "
                "by month: FAIL");
   }

   // Check year rollover for longer month intervals

   TstYear  = 1984;
   TstMonth = 10;
   TstDay   = 15;
   ChkYear  = 1986;
   ChkMonth = 4;
   ChkDay   = 15;

   Err1 =
       CalGreg->incrementDate(18, TimeUnits::Months, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by "
               "18 months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date by "
                "18 months: FAIL");
   }

   TstYear  = 1984;
   TstMonth = 10;
   TstDay   = 15;
   ChkYear  = 1983;
   ChkMonth = 4;
   ChkDay   = 15;

   Err1 = CalGreg->incrementDate(-18, TimeUnits::Months, TstYear, TstMonth,
                                 TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement Gregorian date by "
               "18 months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement Gregorian date by "
                "18 months: FAIL");
   }

   // Check error case when day exceeds max day of new month

   TstYear  = 1984;
   TstMonth = 8;
   TstDay   = 31;

   Err1 =
       CalGreg->incrementDate(1, TimeUnits::Months, TstYear, TstMonth, TstDay);
   if (Err1 != 0) {
      LOG_INFO("TimeMgrTest/Calendar: increment catch bad day range: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment catch bad day range: FAIL");
   }

   // Test normal daily increments/decrements including a leap day

   TstYear  = 1984;
   TstMonth = 2;
   TstDay   = 25;
   ChkYear  = 1984;
   ChkMonth = 3;
   ChkDay   = 6;

   Err1 =
       CalGreg->incrementDate(10, TimeUnits::Days, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by "
               "10 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date by "
                "10 days: FAIL");
   }

   TstYear  = 1984;
   TstMonth = 3;
   TstDay   = 6;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 =
       CalGreg->incrementDate(-10, TimeUnits::Days, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement Gregorian date by "
               "10 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement Gregorian date by "
                "10 days: FAIL");
   }

   // Test longer daily intervals

   TstYear  = 1984;
   TstMonth = 2;
   TstDay   = 25;
   ChkYear  = 1985;
   ChkMonth = 3;
   ChkDay   = 31;

   Err1 =
       CalGreg->incrementDate(400, TimeUnits::Days, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by "
               "400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date by "
                "400 days: FAIL");
   }

   TstYear  = 1985;
   TstMonth = 3;
   TstDay   = 31;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 =
       CalGreg->incrementDate(-400, TimeUnits::Days, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement Gregorian date by "
               "400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement Gregorian date by "
                "400 days: FAIL");
   }

   //------------------------
   // Test No Leap calendar

   Calendar::reset(); // reset calendar
   Calendar::init("No Leap");
   Calendar *CalNoLeap = Calendar::get();

   Kind0 = CalendarNoLeap;
   // Other quantities identical to Gregorian above

   Kind1             = CalendarNoCalendar;
   DaysPerMonth1[0]  = 99;
   DaysPerMonth1[1]  = 99;
   DaysPerMonth1[2]  = 99;
   DaysPerMonth1[3]  = 99;
   DaysPerMonth1[4]  = 99;
   DaysPerMonth1[5]  = 99;
   DaysPerMonth1[6]  = 99;
   DaysPerMonth1[7]  = 99;
   DaysPerMonth1[8]  = 99;
   DaysPerMonth1[9]  = 99;
   DaysPerMonth1[10] = 99;
   DaysPerMonth1[11] = 99;
   MonthsPerYear1    = 999;
   SecondsPerDay1    = 999;
   SecondsPerYear1   = 999;
   DaysPerYear1      = 999;

   Kind1           = Calendar::getKind();
   DaysPerMonth1   = Calendar::getDaysPerMonth();
   MonthsPerYear1  = Calendar::getMonthsPerYear();
   SecondsPerDay1  = Calendar::getSecondsPerDay();
   SecondsPerYear1 = Calendar::getSecondsPerYear();
   DaysPerYear1    = Calendar::getDaysPerYear();

   if (Kind1 != Kind0 || MonthsPerYear1 != MonthsPerYear0 ||
       SecondsPerDay1 != SecondsPerDay0 || SecondsPerYear1 != SecondsPerYear0 ||
       DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: No leap constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: No leap constructor: FAIL");
   }

   // Verify the calendar does not have a leap year
   TstYear = 1984;
   if (!Calendar::isLeapYear(TstYear)) {
      if (Err1 == 0) {
         LOG_INFO("TimeMgrTest/Calendar: 1984 leap year NoLeap: PASS");
      } else {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Calendar: 1984 leap year NoLeap: FAIL");
      }
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: 1984 leap year NoLeap: FAIL");
   }

   // Test calendar date to/from elapsed time conversion
   // This is not an exhaustive test - just testing a single time
   // and for internal consistency

   TmpSeconds =
       (1957 * (I8)86400 * 365) +
       (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 4 - 1) * (I8)86400 +
       19 * 3600 + 26 * 60 + 24;
   Err1       = TstTime.set(TmpSeconds, 1, 4);
   TstYear    = 1957;
   TstMonth   = 10;
   TstDay     = 4;
   TstHour    = 19;
   TstMinute  = 26;
   TstSecondW = 24;
   TstSecondN = 1;
   TstSecondD = 4;

   Err1       = ChkTime.set(0, 0, 1);
   ChkYear    = 0;
   ChkMonth   = 0;
   ChkDay     = 0;
   ChkHour    = 0;
   ChkMinute  = 0;
   ChkSecondW = 0;
   ChkSecondN = 0;
   ChkSecondD = 1;

   ChkTime =
       Calendar::getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert NoLeap date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert NoLeap date to "
                "elapsed time: FAIL");
   }

   Err1 = Calendar::getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
                                ChkMinute, ChkSecondW, ChkSecondN, ChkSecondD);
   if (Err1 == 0 && ChkYear == TstYear && ChkMonth == TstMonth &&
       ChkDay == TstDay && ChkHour == TstHour && ChkMinute == TstMinute &&
       ChkSecondW == TstSecondW && ChkSecondN == TstSecondN &&
       ChkSecondD == TstSecondD) {
      LOG_INFO("TimeMgrTest/Calendar: convert elapsed time to "
               "NoLeap date: PASS");

   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert elapsed time to "
                "NoLeap date: FAIL");
   }

   // Test calendar date increment function
   // Check normal year increment/decrement
   TstYear  = 1983;
   TstMonth = 6;
   TstDay   = 15;
   ChkYear  = 1985;
   ChkMonth = 6;
   ChkDay   = 15;

   Err1 =
       CalNoLeap->incrementDate(2, TimeUnits::Years, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by year: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by year: FAIL");
   }

   TstYear = 1983;
   ChkYear = 1981;
   Err1    = CalNoLeap->incrementDate(-2, TimeUnits::Years, TstYear, TstMonth,
                                      TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement NoLeap date by year: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement NoLeap date by year: FAIL");
   }

   // Check normal month increment/decrement
   TstYear  = 1984;
   TstMonth = 6;
   TstDay   = 15;
   ChkYear  = 1984;
   ChkMonth = 8;
   ChkDay   = 15;

   TstMonth = 6;
   Err1     = CalNoLeap->incrementDate(2, TimeUnits::Months, TstYear, TstMonth,
                                       TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by month: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by month: FAIL");
   }

   TstMonth = 6;
   ChkMonth = 4;
   Err1     = CalNoLeap->incrementDate(-2, TimeUnits::Months, TstYear, TstMonth,
                                       TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement NoLeap date by month: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement NoLeap date by month: FAIL");
   }

   // Test normal daily increments/decrements

   TstYear  = 1984;
   TstMonth = 2;
   TstDay   = 25;
   ChkYear  = 1984;
   ChkMonth = 3;
   ChkDay   = 7;

   Err1 =
       CalNoLeap->incrementDate(10, TimeUnits::Days, TstYear, TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by 10 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by 10 days: FAIL");
   }

   TstYear  = 1984;
   TstMonth = 3;
   TstDay   = 7;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 = CalNoLeap->incrementDate(-10, TimeUnits::Days, TstYear, TstMonth,
                                   TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement NoLeap date by 10 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement NoLeap date by 10 days: FAIL");
   }

   // Test longer daily intervals

   TstYear  = 1984;
   TstMonth = 2;
   TstDay   = 25;
   ChkYear  = 1985;
   ChkMonth = 4;
   ChkDay   = 1;

   Err1 = CalNoLeap->incrementDate(400, TimeUnits::Days, TstYear, TstMonth,
                                   TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by 400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by "
                "400 days: FAIL");
   }

   TstYear  = 1985;
   TstMonth = 4;
   TstDay   = 1;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 = CalNoLeap->incrementDate(-400, TimeUnits::Days, TstYear, TstMonth,
                                   TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement NoLeap date by 400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement NoLeap date by "
                "400 days: FAIL");
   }

   //------------------------
   // Straight Julian calendar the same
   Calendar::reset(); // reset calendar
   Calendar::init("Julian");
   Calendar *CalJulian = Calendar::get();

   Kind1             = CalendarNoCalendar;
   DaysPerMonth1[0]  = 99;
   DaysPerMonth1[1]  = 99;
   DaysPerMonth1[2]  = 99;
   DaysPerMonth1[3]  = 99;
   DaysPerMonth1[4]  = 99;
   DaysPerMonth1[5]  = 99;
   DaysPerMonth1[6]  = 99;
   DaysPerMonth1[7]  = 99;
   DaysPerMonth1[8]  = 99;
   DaysPerMonth1[9]  = 99;
   DaysPerMonth1[10] = 99;
   DaysPerMonth1[11] = 99;
   MonthsPerYear1    = 999;
   SecondsPerDay1    = 999;
   SecondsPerYear1   = 999;
   DaysPerYear1      = 999;

   Kind0 = CalendarJulian;

   Kind1           = Calendar::getKind();
   DaysPerMonth1   = Calendar::getDaysPerMonth();
   MonthsPerYear1  = Calendar::getMonthsPerYear();
   SecondsPerDay1  = Calendar::getSecondsPerDay();
   SecondsPerYear1 = Calendar::getSecondsPerYear();
   DaysPerYear1    = Calendar::getDaysPerYear();

   if (Kind1 != Kind0 || MonthsPerYear1 != MonthsPerYear0 ||
       SecondsPerDay1 != SecondsPerDay0 || SecondsPerYear1 != SecondsPerYear0 ||
       DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Julian constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Julian constructor: FAIL");
   }

   // Leap year calendar normal leap years
   TstYear = 1984;
   if (Calendar::isLeapYear(TstYear)) {
      if (Err1 == 0) {
         LOG_INFO("TimeMgrTest/Calendar: 1984 leap year Julian: PASS");
      } else {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Calendar: 1984 leap year Julian: FAIL");
      }
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: 1984 leap year Julian: FAIL");
   }

   // Test calendar date to/from elapsed time conversion
   // For Julian calendar, don't have a good reference date
   // so just check internal consistency
   TstYear    = 1957;
   TstMonth   = 10;
   TstDay     = 4;
   TstHour    = 19;
   TstMinute  = 26;
   TstSecondW = 24;
   TstSecondN = 1;
   TstSecondD = 4;

   Err1       = ChkTime.set(0, 0, 1);
   ChkYear    = 0;
   ChkMonth   = 0;
   ChkDay     = 0;
   ChkHour    = 0;
   ChkMinute  = 0;
   ChkSecondW = 0;
   ChkSecondN = 0;
   ChkSecondD = 1;

   ChkTime =
       Calendar::getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);

   Err1 = Calendar::getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
                                ChkMinute, ChkSecondW, ChkSecondN, ChkSecondD);
   if (Err1 == 0 && ChkYear == TstYear && ChkMonth == TstMonth &&
       ChkDay == TstDay && ChkHour == TstHour && ChkMinute == TstMinute &&
       ChkSecondW == TstSecondW && ChkSecondN == TstSecondN &&
       ChkSecondD == TstSecondD) {
      LOG_INFO("TimeMgrTest/Calendar: convert elapsed time to "
               "Julian date: PASS");

   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert elapsed time to "
                "Julian date: FAIL");
   }

   //------------------------
   // Test calendar construction for 360 day
   Calendar::reset(); // reset calendar
   Calendar::init("360 Day");
   Calendar *Cal360Day = Calendar::get();

   Kind0             = Calendar360Day;
   DaysPerMonth0[0]  = 30;
   DaysPerMonth0[1]  = 30;
   DaysPerMonth0[2]  = 30;
   DaysPerMonth0[3]  = 30;
   DaysPerMonth0[4]  = 30;
   DaysPerMonth0[5]  = 30;
   DaysPerMonth0[6]  = 30;
   DaysPerMonth0[7]  = 30;
   DaysPerMonth0[8]  = 30;
   DaysPerMonth0[9]  = 30;
   DaysPerMonth0[10] = 30;
   DaysPerMonth0[11] = 30;
   MonthsPerYear0    = 12;
   SecondsPerDay0    = 86400;
   SecondsPerYear0   = 31104000;
   DaysPerYear0      = 360;

   Kind1             = CalendarNoCalendar;
   DaysPerMonth1[0]  = 99;
   DaysPerMonth1[1]  = 99;
   DaysPerMonth1[2]  = 99;
   DaysPerMonth1[3]  = 99;
   DaysPerMonth1[4]  = 99;
   DaysPerMonth1[5]  = 99;
   DaysPerMonth1[6]  = 99;
   DaysPerMonth1[7]  = 99;
   DaysPerMonth1[8]  = 99;
   DaysPerMonth1[9]  = 99;
   DaysPerMonth1[10] = 99;
   DaysPerMonth1[11] = 99;
   MonthsPerYear1    = 999;
   SecondsPerDay1    = 999;
   SecondsPerYear1   = 999;
   DaysPerYear1      = 999;

   Kind1           = Calendar::getKind();
   DaysPerMonth1   = Calendar::getDaysPerMonth();
   MonthsPerYear1  = Calendar::getMonthsPerYear();
   SecondsPerDay1  = Calendar::getSecondsPerDay();
   SecondsPerYear1 = Calendar::getSecondsPerYear();
   DaysPerYear1    = Calendar::getDaysPerYear();

   if (Kind1 != Kind0 || MonthsPerYear1 != MonthsPerYear0 ||
       SecondsPerDay1 != SecondsPerDay0 || SecondsPerYear1 != SecondsPerYear0 ||
       DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: 360 Day constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: 360 Day constructor: FAIL");
   }

   // Test calendar date to/from elapsed time conversion
   // This is not an exhaustive test - just testing a single time
   // and for internal consistency

   TmpSeconds = (1957 * (I8)86400 * 360) + (9 * 30 + 4 - 1) * (I8)86400 +
                19 * 3600 + 26 * 60 + 24;
   Err1       = TstTime.set(TmpSeconds, 1, 4);
   TstYear    = 1957;
   TstMonth   = 10;
   TstDay     = 4;
   TstHour    = 19;
   TstMinute  = 26;
   TstSecondW = 24;
   TstSecondN = 1;
   TstSecondD = 4;

   Err1       = ChkTime.set(0, 0, 1);
   ChkYear    = 0;
   ChkMonth   = 0;
   ChkDay     = 0;
   ChkHour    = 0;
   ChkMinute  = 0;
   ChkSecondW = 0;
   ChkSecondN = 0;
   ChkSecondD = 1;

   ChkTime =
       Calendar::getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert 360-day date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert 360-day date to "
                "elapsed time: FAIL");
   }

   Err1 = Calendar::getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
                                ChkMinute, ChkSecondW, ChkSecondN, ChkSecondD);
   if (Err1 == 0 && ChkYear == TstYear && ChkMonth == TstMonth &&
       ChkDay == TstDay && ChkHour == TstHour && ChkMinute == TstMinute &&
       ChkSecondW == TstSecondW && ChkSecondN == TstSecondN &&
       ChkSecondD == TstSecondD) {
      LOG_INFO("TimeMgrTest/Calendar: convert elapsed time to "
               "360-day date: PASS");

   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert elapsed time to "
                "360-day date: FAIL");
   }

   //------------------------
   // Test calendar construction for Julian day
   Calendar::reset(); // reset calendar
   Calendar::init("Julian Day");
   Calendar *CalJulianDay = Calendar::get();

   Kind0             = CalendarJulianDay;
   DaysPerMonth0[0]  = 0;
   DaysPerMonth0[1]  = 0;
   DaysPerMonth0[2]  = 0;
   DaysPerMonth0[3]  = 0;
   DaysPerMonth0[4]  = 0;
   DaysPerMonth0[5]  = 0;
   DaysPerMonth0[6]  = 0;
   DaysPerMonth0[7]  = 0;
   DaysPerMonth0[8]  = 0;
   DaysPerMonth0[9]  = 0;
   DaysPerMonth0[10] = 0;
   DaysPerMonth0[11] = 0;
   MonthsPerYear0    = 12;
   SecondsPerDay0    = 86400;
   SecondsPerYear0   = 0;
   DaysPerYear0      = 0;

   Kind1             = CalendarNoCalendar;
   DaysPerMonth1[0]  = 99;
   DaysPerMonth1[1]  = 99;
   DaysPerMonth1[2]  = 99;
   DaysPerMonth1[3]  = 99;
   DaysPerMonth1[4]  = 99;
   DaysPerMonth1[5]  = 99;
   DaysPerMonth1[6]  = 99;
   DaysPerMonth1[7]  = 99;
   DaysPerMonth1[8]  = 99;
   DaysPerMonth1[9]  = 99;
   DaysPerMonth1[10] = 99;
   DaysPerMonth1[11] = 99;
   MonthsPerYear1    = 999;
   SecondsPerDay1    = 999;
   SecondsPerYear1   = 999;
   DaysPerYear1      = 999;

   Kind1           = Calendar::getKind();
   DaysPerMonth1   = Calendar::getDaysPerMonth();
   MonthsPerYear1  = Calendar::getMonthsPerYear();
   SecondsPerDay1  = Calendar::getSecondsPerDay();
   SecondsPerYear1 = Calendar::getSecondsPerYear();
   DaysPerYear1    = Calendar::getDaysPerYear();

   if (Kind1 != Kind0 || MonthsPerYear1 != MonthsPerYear0 ||
       SecondsPerDay1 != SecondsPerDay0 || SecondsPerYear1 != SecondsPerYear0 ||
       DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Julian Day constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Julian Day constructor: FAIL");
   }

   // Calendar with no leap year
   TstYear = 1981;
   if (!Calendar::isLeapYear(TstYear)) {
      if (Err1 == 0) {
         LOG_INFO("TimeMgrTest/Calendar: non-leap year JulianDay: PASS");
      } else {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Calendar: non-leap year JulianDay: FAIL");
      }
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: non-leap year JulianDay: FAIL");
   }

   // Test calendar date to/from elapsed time conversion
   // This is not an exhaustive test - just testing a single time
   // and for internal consistency
   // Julian Day - test Julian Day of 2436116.31 (plus 1/4 second)
   TmpSeconds = (2436116 - 1) * (I8)86400 + 7 * 3600 + 26 * 60 + 24;
   Err1       = TstTime.set(TmpSeconds, 1, 4);
   TstYear    = 0;
   TstMonth   = 0;
   TstDay     = 2436116;
   TstHour    = 7;
   TstMinute  = 26;
   TstSecondW = 24;
   TstSecondN = 1;
   TstSecondD = 4;

   Err1       = ChkTime.set(0, 0, 1);
   ChkYear    = 0;
   ChkMonth   = 0;
   ChkDay     = 0;
   ChkHour    = 0;
   ChkMinute  = 0;
   ChkSecondW = 0;
   ChkSecondN = 0;
   ChkSecondD = 1;

   ChkTime =
       Calendar::getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert Julian day to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert Julian day to "
                "elapsed time: FAIL");
   }

   Err1 = Calendar::getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
                                ChkMinute, ChkSecondW, ChkSecondN, ChkSecondD);
   if (Err1 == 0 && ChkYear == TstYear && ChkMonth == TstMonth &&
       ChkDay == TstDay && ChkHour == TstHour && ChkMinute == TstMinute &&
       ChkSecondW == TstSecondW && ChkSecondN == TstSecondN &&
       ChkSecondD == TstSecondD) {
      LOG_INFO("TimeMgrTest/Calendar: convert elapsed time to "
               "Julian day: PASS");

   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert elapsed time to "
                "Julian day: FAIL");
   }

   //------------------------
   // Modified Julian day identical
   Calendar::reset(); // reset calendar
   Calendar::init("Modified Julian Day");
   Calendar *CalModJulianDay = Calendar::get();

   Kind0 = CalendarModJulianDay;

   Kind1             = CalendarNoCalendar;
   DaysPerMonth1[0]  = 99;
   DaysPerMonth1[1]  = 99;
   DaysPerMonth1[2]  = 99;
   DaysPerMonth1[3]  = 99;
   DaysPerMonth1[4]  = 99;
   DaysPerMonth1[5]  = 99;
   DaysPerMonth1[6]  = 99;
   DaysPerMonth1[7]  = 99;
   DaysPerMonth1[8]  = 99;
   DaysPerMonth1[9]  = 99;
   DaysPerMonth1[10] = 99;
   DaysPerMonth1[11] = 99;
   MonthsPerYear1    = 999;
   SecondsPerDay1    = 999;
   SecondsPerYear1   = 999;
   DaysPerYear1      = 999;

   Kind1           = Calendar::getKind();
   DaysPerMonth1   = Calendar::getDaysPerMonth();
   MonthsPerYear1  = Calendar::getMonthsPerYear();
   SecondsPerDay1  = Calendar::getSecondsPerDay();
   SecondsPerYear1 = Calendar::getSecondsPerYear();
   DaysPerYear1    = Calendar::getDaysPerYear();

   if (Kind1 != Kind0 || MonthsPerYear1 != MonthsPerYear0 ||
       SecondsPerDay1 != SecondsPerDay0 || SecondsPerYear1 != SecondsPerYear0 ||
       DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Modified Julian Day constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Modified Julian Day constructor: FAIL");
   }

   // Test calendar date to/from elapsed time conversion
   // This is not an exhaustive test - just testing a single time
   // and for internal consistency
   // Mod Julian Day - test Mod Julian Day of 2436116.31 (plus 1/4 second)
   TmpSeconds = (2436116 - 1) * (I8)86400 + 7 * 3600 + 26 * 60 + 24;
   Err1       = TstTime.set(TmpSeconds, 1, 4);
   TstYear    = 0;
   TstMonth   = 0;
   TstDay     = 2436116;
   TstHour    = 7;
   TstMinute  = 26;
   TstSecondW = 24;
   TstSecondN = 1;
   TstSecondD = 4;

   Err1       = ChkTime.set(0, 0, 1);
   ChkYear    = 0;
   ChkMonth   = 0;
   ChkDay     = 0;
   ChkHour    = 0;
   ChkMinute  = 0;
   ChkSecondW = 0;
   ChkSecondN = 0;
   ChkSecondD = 1;

   ChkTime =
       Calendar::getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert Mod Julian day to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert Mod Julian day to "
                "elapsed time: FAIL");
   }

   Err1 = Calendar::getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
                                ChkMinute, ChkSecondW, ChkSecondN, ChkSecondD);
   if (Err1 == 0 && ChkYear == TstYear && ChkMonth == TstMonth &&
       ChkDay == TstDay && ChkHour == TstHour && ChkMinute == TstMinute &&
       ChkSecondW == TstSecondW && ChkSecondN == TstSecondN &&
       ChkSecondD == TstSecondD) {
      LOG_INFO("TimeMgrTest/Calendar: convert elapsed time to "
               "Mod Julian day: PASS");

   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert elapsed time to "
                "Mod Julian day: FAIL");
   }

   // All done calendar tests - switch to Gregorian calendar before returning
   Calendar::reset();
   Calendar::init("Gregorian");
   return ErrAll;

} // end testCalendar

//------------------------------------------------------------------------------
// TimeInterval test

int testTimeInterval(void) {

   LOG_INFO("TimeMgrTest: TimeInterval tests --------------------------------");

   // Initialize error codes
   I4 Err1{0};
   I4 Err2{0};
   I4 ErrAll{0};

   // Initialize some reference values for the fractional
   // representation of a TimeInterval 5 3/8 seconds.
   I8 WRef{5};
   I8 NRef{3};
   I8 DRef{8};
   R8 RRef{5.375};
   I8 IRef{WRef};

   I8 WTst{5};
   I8 NTst{3};
   I8 DTst{8};
   R8 RTst{5.375};
   I8 ITst{WTst};

   // Test default constructor to create a reference fraction
   // Also implicitly tests one form of the get routine.

   TimeInterval TiRef;

   // Test get function for a fractional second representation
   Err1 = TiRef.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == 0 && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: default constructor and get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: default constructor and get: FAIL");
   }

   // Test constructor from fractional seconds

   TimeInterval TiTstFS(WRef, NRef, DRef);

   Err1 = TiTstFS.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == WRef && NTst == NRef && DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: fractional second constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: fractional second "
                "constructor: FAIL");
   }

   // Can now test assignment and equivalence operator

   TiRef = TiTstFS;
   if (TiTstFS == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(==): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(==): FAIL");
   }

   // Test time interval constructor from real seconds
   TimeInterval TiTstRSec(RRef, TimeUnits::Seconds);

   if (TiTstRSec == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real seconds constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real seconds constructor: FAIL");
   }

   // Test time interval constructor from real minutes
   TimeInterval TiTstRMin((RRef / SECONDS_PER_MINUTE), TimeUnits::Minutes);
   if (TiTstRMin == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real minutes constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real minutes constructor: FAIL");
   }

   // Test time interval constructor from real hours
   TimeInterval TiTstRHour((RRef / SECONDS_PER_HOUR), TimeUnits::Hours);
   if (TiTstRHour == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real hours constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real hours constructor: FAIL");
   }

   // Test time interval constructor from real days
   TimeInterval TiTstRDay(RRef, TimeUnits::Days);

   Err1 = TiTstRDay.get(RTst, TimeUnits::Days);
   if (Err1 == 0 && RTst == 5.0) {
      LOG_INFO("TimeMgrTest/TimeInterval: real days constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real days constructor: FAIL");
   }

   // Test time interval constructor from real months
   TimeInterval TiTstRMonth(RRef, TimeUnits::Months);

   Err1 = TiTstRMonth.get(RTst, TimeUnits::Months);
   if (Err1 == 0 && RTst == 5.0) {
      LOG_INFO("TimeMgrTest/TimeInterval: real months constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real months constructor: FAIL");
   }

   // Test time interval constructor from real years
   TimeInterval TiTstRYear(RRef, TimeUnits::Years);

   Err1 = TiTstRYear.get(RTst, TimeUnits::Years);
   if (Err1 == 0 && RTst == 5.0) {
      LOG_INFO("TimeMgrTest/TimeInterval: real years constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real years constructor: FAIL");
   }

   // Test time interval constructor from integer seconds
   TimeInterval TiTstISec(IRef, TimeUnits::Seconds);

   Err1 = TiTstISec.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == WRef && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer seconds constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer seconds constructor: FAIL");
   }

   // Test time interval constructor from integer minutes
   TimeInterval TiTstIMin(IRef, TimeUnits::Minutes);

   Err1 = TiTstIMin.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == WRef * SECONDS_PER_MINUTE && NTst == 0 &&
       DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer minutes constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer minutes constructor: FAIL");
   }

   // Test non-equivalence comparison operator for TimeInterval

   if (TiTstIMin != TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(!=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(!=): FAIL");
   }

   // Test time interval constructor from integer hours
   TimeInterval TiTstIHour(IRef, TimeUnits::Hours);

   Err1 = TiTstIHour.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == WRef * SECONDS_PER_HOUR && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer hours constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer hours constructor: FAIL");
   }

   // Test calendar-based time interval constructor in years
   TimeInterval TiTstIYear(IRef, TimeUnits::Years);

   Err1 = TiTstIYear.get(ITst, TimeUnits::Years);

   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: year constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: year constructor: FAIL");
   }

   // Test calendar-based time interval constructor in months
   TimeInterval TiTstIMonth(IRef, TimeUnits::Months);

   Err1 = TiTstIMonth.get(ITst, TimeUnits::Months);

   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: month constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: month constructor: FAIL");
   }

   // Test calendar-based time interval constructor in days
   TimeInterval TiTstIDay(IRef, TimeUnits::Days);

   Err1 = TiTstIDay.get(ITst, TimeUnits::Days);

   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: day constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: day constructor: FAIL");
   }

   // Test assignment operator for time interval
   TiTstFS = TiRef;

   if (TiTstFS == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: assignment operator: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: assignment operator: FAIL");
   }

   // Test accessor functions in pairs
   // Test fractional second put/get

   Err1 = TiTstFS.set(WRef, NRef, DRef);
   Err2 = TiTstFS.get(WTst, NTst, DTst);

   if (Err1 == 0 && Err2 == 0 && WTst == WRef && NTst == NRef && DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: fractional second put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: fractional second put/get: FAIL");
   }

   // Test real seconds put/get

   Err1 = TiTstRSec.set(RRef, TimeUnits::Seconds);
   Err2 = TiTstRSec.get(RTst, TimeUnits::Seconds);

   if (Err1 == 0 && Err2 == 0 && RTst == RRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real seconds put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real seconds put/get: FAIL");
   }

   // Test integer seconds put/get

   Err1 = TiTstISec.set(IRef, TimeUnits::Seconds);
   Err2 = TiTstISec.get(ITst, TimeUnits::Seconds);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer seconds put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer seconds put/get: FAIL");
   }

   // Test real hours put/get

   Err1 = TiTstRSec.set(RRef, TimeUnits::Hours);
   Err2 = TiTstRSec.get(RTst, TimeUnits::Hours);

   if (Err1 == 0 && Err2 == 0 && RTst == RRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real hours put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real hours put/get: FAIL");
   }

   // Test integer hours put/get

   Err1 = TiTstISec.set(IRef, TimeUnits::Hours);
   Err2 = TiTstISec.get(ITst, TimeUnits::Hours);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer hours put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer hours put/get: FAIL");
   }

   // Test real minutes put/get

   Err1 = TiTstRSec.set(RRef, TimeUnits::Minutes);
   Err2 = TiTstRSec.get(RTst, TimeUnits::Minutes);

   if (Err1 == 0 && Err2 == 0 && RTst == RRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real minutes put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real minutes put/get: FAIL");
   }

   // Test integer minutes put/get

   Err1 = TiTstISec.set(IRef, TimeUnits::Minutes);
   Err2 = TiTstISec.get(ITst, TimeUnits::Minutes);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer minutes put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer minutes put/get: FAIL");
   }

   // Test integer years put/get

   Err1 = TiTstISec.set(IRef + 1, TimeUnits::Years);
   Err2 = TiTstISec.get(ITst, TimeUnits::Years);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef + 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer years put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer years put/get: FAIL");
   }

   // Test integer months put/get

   Err1 = TiTstISec.set(IRef + 1, TimeUnits::Months);
   Err2 = TiTstISec.get(ITst, TimeUnits::Months);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef + 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer months put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer months put/get: FAIL");
   }

   // Test integer days put/get

   Err1 = TiTstISec.set(IRef + 1, TimeUnits::Days);
   Err2 = TiTstISec.get(ITst, TimeUnits::Days);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef + 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer days put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer days put/get: FAIL");
   }

   // Test < operator (and < part of <= operator) for non-calendars
   // Test both success and failure modes

   Err1 = TiTstFS.set(WRef - 1, NRef - 1, DRef);

   if (TiTstFS < TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<): FAIL");
   }

   if (TiTstFS <= TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=): FAIL");
   }

   if (!(TiRef < TiTstFS)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<): FAIL");
   }

   if (!(TiRef <= TiTstFS)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=): FAIL");
   }

   // Test < operator (and < part of <= operator) for calendar intervals

   Err1 = TiTstIYear.set(IRef, TimeUnits::Years);
   Err1 = TiTstIMonth.set(IRef, TimeUnits::Months);
   Err1 = TiTstIDay.set(IRef, TimeUnits::Days);

   TimeInterval TiTstIYear2(IRef - 1, TimeUnits::Years);
   TimeInterval TiTstIMonth2(IRef - 1, TimeUnits::Months);
   TimeInterval TiTstIDay2(IRef - 1, TimeUnits::Days);

   if (TiTstIYear2 < TiTstIYear) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<) years: FAIL");
   }

   if (TiTstIYear2 <= TiTstIYear) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) years: FAIL");
   }

   if (TiTstIMonth2 < TiTstIMonth) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<) months: FAIL");
   }

   if (TiTstIMonth2 <= TiTstIMonth) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) months: FAIL");
   }

   if (TiTstIDay2 < TiTstIDay) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<) days: FAIL");
   }

   if (TiTstIDay2 <= TiTstIDay) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) days: FAIL");
   }

   // test failure modes

   if (!(TiTstIYear < TiTstIYear2)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<) years: FAIL");
   }

   if (!(TiTstIYear <= TiTstIYear2)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) years: FAIL");
   }

   if (!(TiTstIMonth < TiTstIMonth2)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<) months: FAIL");
   }

   if (!(TiTstIMonth <= TiTstIMonth2)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) months: FAIL");
   }

   if (!(TiTstIDay < TiTstIDay2)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<) days: FAIL");
   }

   if (!(TiTstIDay <= TiTstIDay2)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) days: FAIL");
   }

   // Test > operator (and > part of >= operator) for non-calendars

   if (TiRef > TiTstFS) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>): FAIL");
   }

   if (TiRef >= TiTstFS) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=): FAIL");
   }

   if (!(TiTstFS > TiRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>): FAIL");
   }

   if (!(TiTstFS >= TiRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=): FAIL");
   }

   // Test > operator (and > part of >= operator) for calendar interval

   if (TiTstIYear > TiTstIYear2) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>) years: FAIL");
   }

   if (TiTstIYear >= TiTstIYear2) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) years: FAIL");
   }

   if (TiTstIMonth > TiTstIMonth2) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>) months: FAIL");
   }

   if (TiTstIMonth >= TiTstIMonth2) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) months: FAIL");
   }

   if (TiTstIDay > TiTstIDay2) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>) days: FAIL");
   }

   if (TiTstIDay >= TiTstIDay2) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) days: FAIL");
   }

   // test failure modes

   if (!(TiTstIYear2 > TiTstIYear)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>) years: FAIL");
   }

   if (!(TiTstIYear2 >= TiTstIYear)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) years: FAIL");
   }

   if (!(TiTstIMonth2 > TiTstIMonth)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>) months: FAIL");
   }

   if (!(TiTstIMonth2 >= TiTstIMonth)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) months: FAIL");
   }

   if (!(TiTstIDay2 > TiTstIDay)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>) days: FAIL");
   }

   if (!(TiTstIDay2 >= TiTstIDay)) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) days: FAIL");
   }

   // Test equivalence part of comparisons

   TiTstFS      = TiRef;
   TiTstIYear2  = TiTstIYear;
   TiTstIMonth2 = TiTstIMonth;
   TiTstIDay2   = TiTstIDay;

   if (TiTstFS <= TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=): FAIL");
   }

   if (TiTstFS >= TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=): FAIL");
   }

   if (TiTstIYear2 <= TiTstIYear) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) years: FAIL");
   }

   if (TiTstIYear2 >= TiTstIYear) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) years: FAIL");
   }

   if (TiTstIMonth2 <= TiTstIMonth) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) months: FAIL");
   }

   if (TiTstIMonth2 >= TiTstIMonth) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) months: FAIL");
   }

   if (TiTstIDay2 <= TiTstIDay) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(<=) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(<=) days: FAIL");
   }

   if (TiTstIDay2 >= TiTstIDay) {
      LOG_INFO("TimeMgrTest/TimeInterval: operator(>=) days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: operator(>=) days: FAIL");
   }

   // Test addition operator for interval

   TiTstFS      = TiRef + TiRef;
   TiTstIYear2  = TiTstIYear + TiTstIYear;
   TiTstIMonth2 = TiTstIMonth + TiTstIMonth;
   TiTstIDay2   = TiTstIDay + TiTstIDay;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 10 && NTst == 3 && DTst == 4) {
      LOG_INFO("TimeMgrTest/TimeInterval: addition: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: addition: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: addition years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: addition years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: addition months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: addition months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: addition days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: addition days: FAIL");
   }

   // Test increment

   TiTstFS += TiRef;
   TiTstIYear2 += TiTstIYear;
   TiTstIMonth2 += TiTstIMonth;
   TiTstIDay2 += TiTstIDay;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 16 && NTst == 1 && DTst == 8) {
      LOG_INFO("TimeMgrTest/TimeInterval: increment: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: increment: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == (IRef + IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: increment years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: increment years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == (IRef + IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: increment months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: increment months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == (IRef + IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: increment days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: increment days: FAIL");
   }

   // Test subtraction

   TiTstFS      = TiTstFS - TiRef;
   TiTstIYear2  = TiTstIYear2 - TiTstIYear;
   TiTstIMonth2 = TiTstIMonth2 - TiTstIMonth;
   TiTstIDay2   = TiTstIDay2 - TiTstIDay;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 10 && NTst == 3 && DTst == 4) {
      LOG_INFO("TimeMgrTest/TimeInterval: subtraction: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: subtraction: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: subtraction years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: subtraction years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: subtraction months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: subtraction months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: subtraction days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: subtraction days: FAIL");
   }

   // Test decrement

   TiTstFS -= TiRef;
   TiTstIYear2 -= TiTstIYear;
   TiTstIMonth2 -= TiTstIMonth;
   TiTstIDay2 -= TiTstIDay;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 5 && NTst == 3 && DTst == 8) {
      LOG_INFO("TimeMgrTest/TimeInterval: decrement: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: decrement: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: decrement years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: decrement years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: decrement months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: decrement months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: decrement days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: decrement days: FAIL");
   }

   // Test multiply by integer scalar

   TiTstFS      = TiRef * 3;
   TiTstIYear2  = TiTstIYear * 3;
   TiTstIMonth2 = TiTstIMonth * 3;
   TiTstIDay2   = TiTstIDay * 3;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 16 && NTst == 1 && DTst == 8) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply days: FAIL");
   }

   // Test multiply by integer scalar commutative version

   TiTstFS      = 3 * TiRef;
   TiTstIYear2  = 3 * TiTstIYear;
   TiTstIMonth2 = 3 * TiTstIMonth;
   TiTstIDay2   = 3 * TiTstIDay;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 16 && NTst == 1 && DTst == 8) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm int multiply: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm int multiply: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm int multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm int multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm int multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm int multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm int multiply days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm int multiply days: FAIL");
   }

   // Test multiply by integer scalar in place

   TiTstFS *= 3;
   TiTstIYear2 *= 3;
   TiTstIMonth2 *= 3;
   TiTstIDay2 *= 3;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 48 && NTst == 3 && DTst == 8) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply in place: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply in place: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 9 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply in place years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply in place years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 9 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply in place months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply in place months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 9 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply in place days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply in place days: FAIL");
   }

   // Test multiply by real scalar

   TiTstFS      = TiRef * 3.25;
   TiTstIYear2  = TiTstIYear * 3.25;
   TiTstIMonth2 = TiTstIMonth * 3.25;
   TiTstIDay2   = TiTstIDay * 3.25;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 17 && NTst == 15 && DTst == 32) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply days: FAIL");
   }

   // Test multiply by real scalar commutative version

   TiTstFS      = 3.25 * TiRef;
   TiTstIYear2  = 3.25 * TiTstIYear;
   TiTstIMonth2 = 3.25 * TiTstIMonth;
   TiTstIDay2   = 3.25 * TiTstIDay;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 17 && NTst == 15 && DTst == 32) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm real multiply: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm real multiply: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm real multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm real multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm real multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm real multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm real multiply days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm real multiply days: FAIL");
   }

   // Test multiply by real scalar in place

   TiTstFS *= 3.25;
   TiTstIYear2 *= 3.25;
   TiTstIMonth2 *= 3.25;
   TiTstIDay2 *= 3.25;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 56 && NTst == 99 && DTst == 128) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply in place: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply in place: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 52) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply in place years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply in place years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 52) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply in place months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply in place "
                "months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 52) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply in place days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply in place days: FAIL");
   }

   // Test divide by integer functions

   TiTstFS      = TiRef / 3;
   TiTstIYear2  = TiTstIYear / 3;
   TiTstIMonth2 = TiTstIMonth / 3;
   TiTstIDay2   = TiTstIDay / 3;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 1 && NTst == 19 && DTst == 24) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide days: FAIL");
   }

   // Test divide by integer scalar in place

   TiTstFS /= 3;
   TiTstIYear2 /= 3;
   TiTstIMonth2 /= 3;
   TiTstIDay2 /= 3;

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == 0 && NTst == 43 && DTst == 72) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide in place: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide in place: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == 0) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide in place years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide in place years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == 0) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide in place months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide in place months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == 0) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide in place days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide in place days: FAIL");
   }

   // Test negative absolute value first to change sign

   TiTstFS      = TimeInterval::negAbsValue(TiRef);
   TiTstIYear2  = TimeInterval::negAbsValue(TiTstIYear);
   TiTstIMonth2 = TimeInterval::negAbsValue(TiTstIMonth);
   TiTstIDay2   = TimeInterval::negAbsValue(TiTstIDay);

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == -WRef && NTst == -NRef && DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: negAbsValue: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: negAbsValue: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == -IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: negAbsValue years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: negAbsValue years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == -IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: negAbsValue months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: negAbsValue months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == -IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: negAbsValue days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: negAbsValue days: FAIL");
   }

   // Test is positive function with negative results

   if (!(TiTstFS.isPositive()) && !(TiTstIYear2.isPositive()) &&
       !(TiTstIDay2.isPositive())) {
      LOG_INFO("TimeMgrTest/TimeInterval: not isPositive: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: not isPositive: FAIL");
   }

   // Test absolute value by changing the above to abs value

   TimeInterval TiTst3;
   TimeInterval TiTstIYear3;
   TimeInterval TiTstIMonth3;
   TimeInterval TiTstIDay3;

   TiTst3       = TimeInterval::absValue(TiTstFS);
   TiTstIYear3  = TimeInterval::absValue(TiTstIYear2);
   TiTstIMonth3 = TimeInterval::absValue(TiTstIMonth2);
   TiTstIDay3   = TimeInterval::absValue(TiTstIDay2);

   Err1 = TiTst3.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == WRef && NTst == NRef && DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: absValue: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: absValue: FAIL");
   }

   Err1 = TiTstIYear3.get(ITst, TimeUnits::Years);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: absValue years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: absValue years: FAIL");
   }

   Err1 = TiTstIMonth3.get(ITst, TimeUnits::Months);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: absValue months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: absValue months: FAIL");
   }

   Err1 = TiTstIDay3.get(ITst, TimeUnits::Days);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: absValue days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: absValue days: FAIL");
   }

   // Test is positive function with positive results

   if (TiTst3.isPositive() && TiTstIYear3.isPositive() &&
       TiTstIMonth3.isPositive() && TiTstIDay3.isPositive()) {
      LOG_INFO("TimeMgrTest/TimeInterval: isPositive: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: isPositive: FAIL");
   }

   // Test time interval constructor from string
   std::string TiStr = "1001_11:23:45.375";
   TimeInterval TiTstStr(TiStr);
   RRef = 86527425.375;
   Err1 = TiTstStr.get(RTst, TimeUnits::Seconds);

   if (Err1 == 0 && RTst == RRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: string constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: string constructor: FAIL");
   }

   return ErrAll;

} // end testTimeInterval

//------------------------------------------------------------------------------
// TimeInstant test

int testTimeInstant(void) {

   LOG_INFO("TimeMgrTest: TimeInstant tests ---------------------------------");

   // Initialize error codes
   I4 Err1{0};
   I4 Err2{0};
   I4 Err3{0};
   I4 ErrAll{0};

   // Use default constructor to create first (empty) instant
   TimeInstant TiEmpty;

   // Create some reference elapsed times based on 4 most likely calendars
   // (Gregorian, NoLeap, 360 day and no-calendar). And use reference date
   // of July 4, 2019 at 3:16:23.25.

   I8 YearRef{2019};
   I8 MonthRef{7};
   I8 DayRef{4};
   I8 HourRef{15};
   I8 MinuteRef{16};
   I8 WRef{23};
   I8 NRef{1};
   I8 DRef{4};
   R8 RRef{23.25};

   //--- Gregorian calendar tests
   // Construct an instant for the most likely Gregorian use case
   // Then test using get functions

   Calendar::reset(); // reset calendar
   Calendar::init("Gregorian");

   TimeInstant TiGreg(YearRef, MonthRef, DayRef, HourRef, MinuteRef, RRef);

   I8 YearChk    = 0;
   I8 MonthChk   = 0;
   I8 DayChk     = 0;
   I8 HourChk    = 0;
   I8 MinuteChk  = 0;
   R8 RSecondChk = 0.0;

   Err1 = TiGreg.get(YearChk, MonthChk, DayChk, HourChk, MinuteChk, RSecondChk);
   if (Err1 == 0 && YearChk == YearRef && MonthChk == MonthRef &&
       DayChk == DayRef && HourChk == HourRef && MinuteChk == MinuteRef &&
       abs(RSecondChk - RRef) < 1.e-15) {
      LOG_INFO("TimeMgrTest/TimeInstant: constructor YMDHMS(real): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: constructor YMDHMS(real): FAIL");
   }

   // Now use set function to create an identical instant

   TimeInstant TiGreg2;
   Err1 = TiGreg2.set(YearRef, MonthRef, DayRef, HourRef, MinuteRef, RRef);

   YearChk    = 0;
   MonthChk   = 0;
   DayChk     = 0;
   HourChk    = 0;
   MinuteChk  = 0;
   RSecondChk = 0.0;

   Err2 =
       TiGreg2.get(YearChk, MonthChk, DayChk, HourChk, MinuteChk, RSecondChk);
   if (Err1 == 0 && Err2 == 0 && YearChk == YearRef && MonthChk == MonthRef &&
       DayChk == DayRef && HourChk == HourRef && MinuteChk == MinuteRef &&
       abs(RSecondChk - RRef) < 1.e-15) {
      LOG_INFO("TimeMgrTest/TimeInstant: get/set YMDHMS(real): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: get/set YMDHMS(real): FAIL");
   }

   // Can now also check equivalence and equivalence part of >=, <=

   if (TiGreg == TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(==): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(==): FAIL");
   }

   if (TiGreg >= TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(>=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(>=): FAIL");
   }

   if (TiGreg <= TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(<=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(<=): FAIL");
   }

   // Now must test remaining operators involving time intervals

   // Test increment/decrement by time interval
   TimeInterval IntervalSec(1, TimeUnits::Seconds);

   TiGreg2             = TiGreg + IntervalSec;
   TimeInstant TiGreg3 = TiGreg2 - IntervalSec;

   YearChk   = 0;
   MonthChk  = 0;
   DayChk    = 0;
   HourChk   = 0;
   MinuteChk = 0;
   I8 WChk   = 0;
   I8 NChk   = 0;
   I8 DChk   = 0;
   Err1 = TiGreg2.get(YearChk, MonthChk, DayChk, HourChk, MinuteChk, WChk, NChk,
                      DChk);

   if (Err1 == 0 && YearChk == YearRef && MonthChk == MonthRef &&
       DayChk == DayRef && HourChk == HourRef && MinuteChk == MinuteRef &&
       NChk == NRef && DChk == DRef && WChk == WRef + 1) {
      LOG_INFO("TimeMgrTest/TimeInstant: addition second interval: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: addition second interval: FAIL");
   }

   if (TiGreg3 == TiGreg) {
      LOG_INFO("TimeMgrTest/TimeInstant: subtraction second interval: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: subtraction second interval: FAIL");
   }

   // Test a 5-year integration in several units (year, month, day, hour,
   // minute). Include a nominal leap year.

   TimeInterval IntervalYear5(5, TimeUnits::Years);
   TimeInterval IntervalYear(1, TimeUnits::Years);
   TimeInterval IntervalMonth(2, TimeUnits::Months);
   TimeInterval IntervalDay(1, TimeUnits::Days);
   TimeInterval IntervalHour(2, TimeUnits::Hours);
   TimeInterval IntervalMinute(20, TimeUnits::Minutes);

   // Add the five year interval to create a final target date
   TiGreg2 = TiGreg + IntervalYear5;

   // Test intervals for Gregorian calendars
   TimeInstant TiFinal = TiGreg;
   for (int N = 1; N <= 5; ++N) {
      TiFinal += IntervalYear;
   }
   if (TiFinal == TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: Gregorian annual integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: Gregorian annual integration: FAIL");
   }

   TiFinal = TiGreg;
   for (int N = 1; N <= 30; ++N) {
      TiFinal += IntervalMonth;
   }
   if (TiFinal == TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: Gregorian monthly integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: Gregorian monthly integration: FAIL");
   }

   TiFinal = TiGreg;
   for (int N = 1; N <= 365 * 5 + 2; ++N) { // period includes 2 leap years
      TiFinal += IntervalDay;
   }
   if (TiFinal == TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: Gregorian daily integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: Gregorian daily integration: FAIL");
   }

   TiFinal = TiGreg;
   for (int N = 1; N <= 12 * (365 * 5 + 2); ++N) {
      TiFinal += IntervalHour;
   }
   if (TiFinal == TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: Gregorian hourly integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: Gregorian hourly integration: FAIL");
   }

   TiFinal = TiGreg;
   for (int N = 1; N <= (3 * 24) * (365 * 5 + 2); ++N) {
      TiFinal += IntervalMinute;
   }
   if (TiFinal == TiGreg2) {
      LOG_INFO("TimeMgrTest/TimeInstant: Gregorian minute integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: Gregorian minute integration: FAIL");
   }

   // Test time string generator and constructor from string

   std::string StrDateRef = "2019-07-04_15:16:23.2500";
   std::string StrDateChk = TiGreg.getString(4, 4, "_");

   if (StrDateChk == StrDateRef) {
      LOG_INFO("TimeMgrTest/TimeInstant: getString: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: getString: FAIL");
   }

   TimeInstant TiFromString(StrDateChk);
   if (TiFromString == TiGreg) {
      LOG_INFO("TimeMgrTest/TimeInstant from string: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant from string: FAIL");
   }

   //--- Test instants in no-leap calendar
   Calendar::reset(); // reset calendar
   Calendar::init("No Leap");

   // Construct a no-leap instant using frac second interface
   TimeInstant TiNoLeap(YearRef, MonthRef, DayRef, HourRef, MinuteRef, WRef,
                        NRef, DRef);

   YearChk   = 0;
   MonthChk  = 0;
   DayChk    = 0;
   HourChk   = 0;
   MinuteChk = 0;
   WChk      = 0;
   NChk      = 0;
   DChk      = 0;

   Err1 = TiNoLeap.get(YearChk, MonthChk, DayChk, HourChk, MinuteChk, WChk,
                       NChk, DChk);
   if (Err1 == 0 && YearChk == YearRef && MonthChk == MonthRef &&
       DayChk == DayRef && HourChk == HourRef && MinuteChk == MinuteRef &&
       WChk == WRef && NChk == NRef && DChk == DRef) {
      LOG_INFO("TimeMgrTest/TimeInstant: constructor YMDHMS(frac): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: constructor YMDHMS(frac): FAIL");
   }

   // Now use get/set interface to set a slightly earlier instant

   TimeInstant TiNoLeap2;
   Err2 = TiNoLeap2.set(YearRef, MonthRef, DayRef, HourRef, MinuteRef, WRef - 1,
                        NRef, DRef);

   YearChk   = 0;
   MonthChk  = 0;
   DayChk    = 0;
   HourChk   = 0;
   MinuteChk = 0;
   WChk      = 0;
   NChk      = 0;
   DChk      = 0;
   Err3 = TiNoLeap2.get(YearChk, MonthChk, DayChk, HourChk, MinuteChk, WChk,
                        NChk, DChk);

   if (Err2 == 0 && Err3 == 0 && YearChk == YearRef && MonthChk == MonthRef &&
       DayChk == DayRef && HourChk == HourRef && MinuteChk == MinuteRef &&
       WChk == WRef - 1 && NChk == NRef && DChk == DRef) {
      LOG_INFO("TimeMgrTest/TimeInstant: get/set by YMDHMS(frac): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: get/set by YMDHMS(frac): FAIL");
   }

   // Can use these to test a few more operators like non-equivalence, <

   // Non-equiv for different time instant in same calendar
   if (TiNoLeap != TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(!=) for "
               "different time instant: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(!=) for "
                "different time instant: FAIL");
   }

   // Test forms of > operator
   if (TiNoLeap >= TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(>=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(>=): FAIL");
   }

   if (TiNoLeap > TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(>): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(>): FAIL");
   }

   // Test forms of < operator
   if (TiNoLeap2 <= TiNoLeap) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(<=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(<=): FAIL");
   }

   if (TiNoLeap2 < TiNoLeap) {
      LOG_INFO("TimeMgrTest/TimeInstant: operator(<): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: operator(<): FAIL");
   }

   // NoLeap and NoLeap2 above differ by one second, so create
   // a 1-second interval using the difference between the two instants
   // and check with the previously defined 1-second interval

   TimeInterval IntervalSec2 = TiNoLeap - TiNoLeap2;

   if (IntervalSec == IntervalSec2) {
      LOG_INFO("TimeMgrTest/TimeInstant: create interval from "
               "diff of instants: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: create interval from "
                "diff of instants: FAIL");
   }

   // Now test addition, subtraction for seconds

   // Test increment and decrement using the second interval and noLeap
   TimeInstant TiNoLeap3 = TiNoLeap2;

   TiNoLeap2 += IntervalSec;
   if (TiNoLeap2 == TiNoLeap) {
      LOG_INFO("TimeMgrTest/TimeInstant: increment by second interval: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: increment by second interval: FAIL");
   }

   TiNoLeap2 -= IntervalSec;
   if (TiNoLeap2 == TiNoLeap3) {
      LOG_INFO("TimeMgrTest/TimeInstant: decrement by second interval: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: decrement by second interval: FAIL");
   }

   // Test a 5-year integration in several units (year, month, day, hour,
   // minute). Include a nominal leap year.

   TimeInterval IntYear5NL(5, TimeUnits::Years);
   TimeInterval IntYearNL(1, TimeUnits::Years);
   TimeInterval IntMonthNL(2, TimeUnits::Months);
   TimeInterval IntDayNL(1, TimeUnits::Days);
   TimeInterval IntHourNL(2, TimeUnits::Hours);
   TimeInterval IntMinuteNL(20, TimeUnits::Minutes);
   // Add the five year interval to create a final target for each calendar
   TiNoLeap2 = TiNoLeap + IntYear5NL;

   // Test integration for each interval
   TiFinal = TiNoLeap;
   for (int N = 1; N <= 5; ++N) {
      TiFinal += IntYearNL;
   }
   if (TiFinal == TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: NoLeap annual integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: NoLeap annual integration: FAIL");
   }

   TiFinal = TiNoLeap;
   for (int N = 1; N <= 30; ++N) {
      TiFinal += IntMonthNL;
   }
   if (TiFinal == TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: NoLeap monthly integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: NoLeap monthly integration: FAIL");
   }

   TiFinal = TiNoLeap;
   for (int N = 1; N <= 365 * 5; ++N) {
      TiFinal += IntDayNL;
   }
   if (TiFinal == TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: NoLeap daily integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: NoLeap daily integration: FAIL");
   }

   TiFinal = TiNoLeap;
   for (int N = 1; N <= 12 * 365 * 5; ++N) {
      TiFinal += IntHourNL;
   }
   if (TiFinal == TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: NoLeap hourly integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: NoLeap hourly integration: FAIL");
   }

   TiFinal = TiNoLeap;
   for (int N = 1; N <= (3 * 24) * (365 * 5); ++N) {
      TiFinal += IntMinuteNL;
   }
   if (TiFinal == TiNoLeap2) {
      LOG_INFO("TimeMgrTest/TimeInstant: NoLeap minute integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: NoLeap minute integration: FAIL");
   }

   //--- Test instants in 360 day calendar
   Calendar::reset(); // reset calendar
   Calendar::init("360 Day");

   // Construct a 360Day instant using frac second interface
   TimeInstant Ti360Day(YearRef, MonthRef, DayRef, HourRef, MinuteRef, WRef,
                        NRef, DRef);

   YearChk   = 0;
   MonthChk  = 0;
   DayChk    = 0;
   HourChk   = 0;
   MinuteChk = 0;
   WChk      = 0;
   NChk      = 0;
   DChk      = 0;
   Err1      = Ti360Day.get(YearChk, MonthChk, DayChk, HourChk, MinuteChk, WChk,
                            NChk, DChk);

   if (Err1 == 0 && YearChk == YearRef && MonthChk == MonthRef &&
       DayChk == DayRef && HourChk == HourRef && MinuteChk == MinuteRef &&
       WChk == WRef && NChk == NRef && DChk == DRef) {
      LOG_INFO("TimeMgrTest/TimeInstant: construct 360Day: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: construct 360Day: FAIL");
   }

   // Test a 5-year integration in several units (year, month, day, hour,
   // minute).

   TimeInterval IntYear5360(5, TimeUnits::Years);
   TimeInterval IntYear360(1, TimeUnits::Years);
   TimeInterval IntMonth360(2, TimeUnits::Months);
   TimeInterval IntDay360(1, TimeUnits::Days);
   TimeInterval IntHour360(2, TimeUnits::Hours);
   TimeInterval IntMinute360(20, TimeUnits::Minutes);
   // Add the five year interval to create a final target
   TimeInstant Ti360Day2 = Ti360Day + IntYear5360;

   // Test intervals for 360Day calendars
   TiFinal = Ti360Day;
   for (int N = 1; N <= 5; ++N) {
      TiFinal += IntYear360;
   }
   if (TiFinal == Ti360Day2) {
      LOG_INFO("TimeMgrTest/TimeInstant: 360Day annual integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: 360Day annual integration: FAIL");
   }

   TiFinal = Ti360Day;
   for (int N = 1; N <= 30; ++N) {
      TiFinal += IntMonth360;
   }
   if (TiFinal == Ti360Day2) {
      LOG_INFO("TimeMgrTest/TimeInstant: 360Day monthly integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: 360Day monthly integration: FAIL");
   }

   TiFinal = Ti360Day;
   for (int N = 1; N <= 360 * 5; ++N) { // period includes 2 leap years
      TiFinal += IntDay360;
   }
   if (TiFinal == Ti360Day2) {
      LOG_INFO("TimeMgrTest/TimeInstant: 360Day daily integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: 360Day daily integration: FAIL");
   }

   TiFinal = Ti360Day;
   for (int N = 1; N <= 12 * 360 * 5; ++N) {
      TiFinal += IntHour360;
   }
   if (TiFinal == Ti360Day2) {
      LOG_INFO("TimeMgrTest/TimeInstant: 360Day hourly integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: 360Day hourly integration: FAIL");
   }

   TiFinal = Ti360Day;
   for (int N = 1; N <= (3 * 24) * (360 * 5); ++N) {
      TiFinal += IntMinute360;
   }
   if (TiFinal == Ti360Day2) {
      LOG_INFO("TimeMgrTest/TimeInstant: 360Day minute integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: 360Day minute integration: FAIL");
   }

   //--- Test instants in no calendar option
   Calendar::reset(); // reset calendar
   Calendar::init("No Calendar");

   // Create 2 no-calendar time instants using different elapsed time
   // constructors

   YearChk   = 0;
   MonthChk  = 0;
   DayChk    = 0;
   HourChk   = 0;
   MinuteChk = 0;
   TimeInstant TiNone(YearChk, MonthChk, DayChk, HourChk, MinuteChk, WRef, NRef,
                      DRef);

   TimeInstant TiNone2(YearChk, MonthChk, DayChk, HourChk, MinuteChk, RRef);

   if (TiNone == TiNone2) {
      LOG_INFO("TimeMgrTest/TimeInstant: No-calendar time constructors: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: No-calendar time constructors: FAIL");
   }

   // Test a five-year integration with only units that make sense when there
   // is no calendar.

   TimeInterval IntHourNone(2, TimeUnits::Hours);
   TimeInterval IntMinuteNone(20, TimeUnits::Minutes);
   TimeInterval IntSeconds5yr(86400 * 365 * 5, TimeUnits::Seconds);

   // Step forward in time
   // In this case, annual, monthly or daily intervals are meaningless
   TiFinal = TiNone;
   TiNone2 = TiNone + IntSeconds5yr;

   for (int N = 1; N <= 12 * 365 * 5; ++N) {
      TiFinal += IntHourNone;
   }
   if (TiFinal == TiNone2) {
      LOG_INFO("TimeMgrTest/TimeInstant: No-calendar hourly integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: No-calendar "
                "hourly integration: FAIL");
   }

   TiFinal = TiNone;
   for (int N = 1; N <= (3 * 24) * (365 * 5); ++N) {
      TiFinal += IntMinuteNone;
   }
   if (TiFinal == TiNone2) {
      LOG_INFO("TimeMgrTest/TimeInstant: No-calendar minute integration: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInstant: No-calendar "
                "minute integration: FAIL");
   }

   //--- Finished time instant tests, switch back to Gregorian calendar
   Calendar::reset(); // reset calendar
   Calendar::init("Gregorian");
   return ErrAll;

} // end testTimeInstant

//------------------------------------------------------------------------------
// Alarm test

int testAlarm(void) {

   LOG_INFO("TimeMgrTest: Alarm tests ---------------------------------------");

   // Initialize error codes
   I4 Err1{0};
   I4 Err2{0};
   I4 Err3{0};
   I4 ErrAll{0};

   // For various time intervals, we create alarms at relevant times
   // and step through time to trigger the alarm.

   // Do all testing in Gregorian calendar

   Calendar::reset();
   Calendar::init("Gregorian");

   // Create a zero start time and generic start time
   TimeInstant Time0(2000, 1, 1, 0, 0, 0.0);
   TimeInstant StartTime(2019, 8, 15, 14, 25, 23.25);

   // Test a year-based alarm and periodic alarm using a start time and
   // monthly interval

   TimeInstant TimeNewYear2020(2020, 1, 1, 0, 0, 0.0);
   Alarm AlarmNewYear2020("New Year 2020", TimeNewYear2020);

   TimeInterval IntervalAnnual(1, TimeUnits::Years);
   TimeInterval IntervalMonthly(1, TimeUnits::Months);
   Alarm AlarmEveryYear("Every Year", IntervalAnnual, Time0);

   TimeInstant CurTime = StartTime;

   // Test retrieval of the time interval
   const TimeInterval *NewInterval = AlarmEveryYear.getInterval();
   if (*NewInterval == IntervalAnnual) {
      LOG_INFO("TimeMgrTest/Alarm: retrieve interval: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Alarm: retrieve interval: FAIL");
   }

   // Test initial retrieval of the previous ring time (should be time0)
   TimeInstant PriorRingTimeTest    = Time0;
   const TimeInstant *PriorRingTime = AlarmEveryYear.getRingTimePrev();
   if (*PriorRingTime == PriorRingTimeTest) {
      LOG_INFO("TimeMgrTest/Alarm: retrieve initial ring time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Alarm: retrieve initial ring time: FAIL");
   }

   // quick test of update status function
   Err1 = AlarmNewYear2020.updateStatus(CurTime);
   Err2 = AlarmEveryYear.updateStatus(CurTime);
   if (Err1 == 0 && Err2 == 0) {
      LOG_INFO("TimeMgrTest/Alarm: update status: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Alarm: update status: FAIL");
   }

   // reset interval timer to make sure first ring time is in future
   Err1 = AlarmEveryYear.reset(CurTime);
   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Alarm: initial reset: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Alarm: initial reset: FAIL");
   }
   // reset the previous ring time to proper value corresponding to current time
   TimeInstant TestTime = Time0;
   while (TestTime < CurTime) {
      PriorRingTimeTest = TestTime;
      TestTime += IntervalAnnual;
   }

   // now integrate forward for 18 months
   for (int N = 1; N <= 18; ++N) {
      // increment time in monthly intervals
      CurTime += IntervalMonthly;

      // update alarm state based on current time
      Err1 = AlarmNewYear2020.updateStatus(CurTime);
      Err2 = AlarmEveryYear.updateStatus(CurTime);
      if (Err1 != 0 || Err2 != 0) {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Alarm: update annual alarms: FAIL");
      }

      // Test whether one-time alarm should be ringing or not
      if (N == 5) {
         if (AlarmNewYear2020.isRinging()) {
            LOG_INFO("TimeMgrTest/Alarm: one-time annual alarm: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time annual alarm: FAIL");
         }
         Err1 = AlarmNewYear2020.stop();
         if (Err1 == 0) {
            LOG_INFO("TimeMgrTest/Alarm: one-time annual alarm stop: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time annual alarm stop: FAIL");
         }
      } else {
         if (AlarmNewYear2020.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time annual alarm "
                      "should not be ringing: FAIL");
         }
      }

      // Test retrieval of previous ring time
      const TimeInstant *PriorRingTime = AlarmEveryYear.getRingTimePrev();
      if (*PriorRingTime != PriorRingTimeTest) {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Alarm: retrieve prior ring time: FAIL");
      }

      // Test whether interval alarm should be ringing or not
      if (N == 5 || N == 17) {
         if (AlarmEveryYear.isRinging()) {
            LOG_INFO("TimeMgrTest/Alarm: periodic annual alarm: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic annual alarm: FAIL");
         }
         Err1 = AlarmEveryYear.reset(CurTime);
         if (Err1 == 0) {
            LOG_INFO("TimeMgrTest/Alarm: periodic annual alarm reset: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic annual alarm reset: FAIL");
         }
         PriorRingTimeTest += IntervalAnnual;
      } else {
         if (AlarmEveryYear.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic annual alarm "
                      "should not be ringing: FAIL");
         }
      }
   }

   // Test a month-based alarm and periodic alarm using a start time and
   // daily interval (pick a leap year just to catch an edge case)

   TimeInstant Time2020Mar1(2020, 3, 1, 0, 0, 0.0);
   Alarm Alarm2020Mar1("2020-03-01", Time2020Mar1);

   Alarm AlarmEveryMonth("Every Month", IntervalMonthly, Time0);
   TimeInterval IntervalDaily(1, TimeUnits::Days);

   CurTime = StartTime; // start time is 2019-08-15_14:25:23.25
   Err1    = AlarmEveryMonth.reset(CurTime); // ensure first alarm in future

   Err3 = 1; // use to limit messages
   for (int N = 1; N <= 365; ++N) {
      // increment time in daily intervals
      CurTime += IntervalDaily;

      // update alarm state based on current time
      Err1 = Alarm2020Mar1.updateStatus(CurTime);
      Err2 = AlarmEveryMonth.updateStatus(CurTime);
      if (Err1 != 0 || Err2 != 0) {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Alarm: update monthly alarms: FAIL");
      }
      // Test whether one-time alarm should be ringing or not
      if (N == 199) {
         if (Alarm2020Mar1.isRinging()) {
            LOG_INFO("TimeMgrTest/Alarm: one-time monthly alarm: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time monthly alarm: FAIL");
         }
         Err1 = Alarm2020Mar1.stop();
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time monthly alarm stop: FAIL");
         }
      } else {
         if (Alarm2020Mar1.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time monthly alarm "
                      "should not be ringing: FAIL");
         }
      }

      // Test whether interval alarm should be ringing or not
      if (N == 17 || N == 47 || N == 78 || N == 108 || N == 139 || N == 170 ||
          N == 199 || N == 230 || N == 260 || N == 291 || N == 321 ||
          N == 352) {
         if (AlarmEveryMonth.isRinging()) {
            if (Err3 == 1) // only print first instance for pass
               LOG_INFO("TimeMgrTest/Alarm: periodic monthly alarm: PASS");
            Err3 = 0;
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic monthly alarm: FAIL");
         }
         Err1 = AlarmEveryMonth.reset(CurTime);
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic monthly alarm reset: FAIL");
         }
      } else {
         if (AlarmEveryMonth.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic monthly alarm "
                      "should not be ringing: FAIL");
         }
      }
   }

   // Test a day-based alarm and periodic alarm using a start time and
   // hourly interval

   TimeInstant Time2019Aug20(2019, 8, 20, 0, 0, 0.0);
   Alarm Alarm2019Aug20("2020-08-20", Time2019Aug20);

   Alarm AlarmEveryDay("Every Day", IntervalDaily, Time0);
   TimeInterval IntervalHourly(1, TimeUnits::Hours);

   CurTime = StartTime; // start time is 2019-08-15_14:25:23.25
   Err1    = AlarmEveryDay.reset(CurTime); // ensure first alarm in future

   Err3 = 1; // use to limit output
   for (int N = 1; N <= 240; ++N) {
      // increment time in hourly intervals
      CurTime += IntervalHourly;

      // update alarm state based on current time
      Err1 = Alarm2019Aug20.updateStatus(CurTime);
      Err2 = AlarmEveryDay.updateStatus(CurTime);
      if (Err1 != 0 || Err2 != 0) {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Alarm: update daily alarms: FAIL");
      }

      // Test whether one-time alarm should be ringing or not
      if (N == 106) {
         if (Alarm2019Aug20.isRinging()) {
            LOG_INFO("TimeMgrTest/Alarm: one-time daily alarm: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time daily alarm: FAIL");
         }
         Err1 = Alarm2019Aug20.stop();
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time daily alarm stop: FAIL");
         }
      } else {
         if (Alarm2019Aug20.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time daily alarm "
                      "should not be ringing: FAIL");
         }
      }

      // Test whether interval alarm should be ringing or not
      if (N == 10 || N == 34 || N == 58 || N == 82 || N == 106 || N == 130 ||
          N == 154 || N == 178 || N == 202 || N == 226 || N == 250) {
         if (AlarmEveryDay.isRinging()) {
            if (Err3 == 1) // only print first instance for pass
               LOG_INFO("TimeMgrTest/Alarm: periodic daily alarm: PASS");
            Err3 = 0;
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic daily alarm: FAIL");
         }
         Err1 = AlarmEveryDay.reset(CurTime);
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic daily alarm reset: FAIL");
         }
      } else {
         if (AlarmEveryMonth.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic monthly alarm "
                      "should not be ringing: FAIL");
         }
      }
   }

   // Test an hour-based alarm and periodic alarm using a start time and
   // minute interval

   TimeInstant Time9am(2019, 8, 16, 9, 0, 0.0);
   Alarm Alarm9am("2019-08-16_0900", Time9am);

   Alarm AlarmEveryHour("Every Hour", IntervalHourly, Time0);
   TimeInterval IntervalMinute(1, TimeUnits::Minutes);

   CurTime = StartTime; // start time is 2019-08-15_14:25:23.25
   Err1    = AlarmEveryHour.reset(CurTime); // ensure first alarm in future

   Err3 = 1; // use to limit number of output lines
   for (int N = 1; N <= 2880; ++N) {
      // increment time in minute intervals
      CurTime += IntervalMinute;

      // update alarm state based on current time
      Err1 = Alarm9am.updateStatus(CurTime);
      Err2 = AlarmEveryHour.updateStatus(CurTime);
      if (Err1 != 0 || Err2 != 0) {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Alarm: update hourly alarms: FAIL");
      }

      // Test whether one-time alarm should be ringing or not
      if (N == 1115) {
         if (Alarm9am.isRinging()) {
            LOG_INFO("TimeMgrTest/Alarm: one-time hourly alarm: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time hourly alarm: FAIL");
         }
         Err1 = Alarm9am.stop();
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time hourly alarm stop: FAIL");
         }
      } else {
         if (Alarm9am.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time hourly alarm "
                      "should not be ringing: FAIL");
         }
      }

      // Test whether interval alarm should be ringing or not
      if ((N - 35) % 60 == 0) {
         if (AlarmEveryHour.isRinging()) {
            if (Err3 == 1) // only print success on first instance
               LOG_INFO("TimeMgrTest/Alarm: periodic hourly alarm: PASS");
            Err3 = 0;
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic hourly alarm: FAIL");
         }
         Err1 = AlarmEveryHour.reset(CurTime);
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic hourly alarm reset: FAIL");
         }
      } else {
         if (AlarmEveryHour.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic hourly alarm "
                      "should not be ringing: FAIL");
         }
      }
   }

   // Test a 6-hour periodic alarm using a start time and
   // hourly interval to test a non-unit time interval

   TimeInterval Interval6Hour(6, TimeUnits::Hours);
   Alarm AlarmEvery6Hour("Every 6 Hours", Interval6Hour, Time0);

   CurTime = StartTime; // start time is 2019-08-15_14:25:23.25
   Err1    = AlarmEvery6Hour.reset(CurTime); // ensure first alarm in future

   Err3 = 1; // limit output
   for (int N = 1; N <= 120; ++N) {
      // increment time in hourly intervals
      CurTime += IntervalHourly;

      // update alarm state based on current time
      Err2 = AlarmEvery6Hour.updateStatus(CurTime);
      if (Err2 != 0) {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Alarm: update 6-hourly alarms: FAIL");
      }

      // Test whether interval alarm should be ringing or not
      if ((N - 4) % 6 == 0) {
         if (AlarmEvery6Hour.isRinging()) {
            if (Err3 == 1) // only print first instance of pass
               LOG_INFO("TimeMgrTest/Alarm: periodic 6-hourly alarm: PASS");
            Err3 = 0;
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic 6-hourly alarm: FAIL");
         }
         Err1 = AlarmEvery6Hour.reset(CurTime);
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic 6-hourly alarm reset: FAIL");
         }
      } else {
         if (AlarmEvery6Hour.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic 6-hourly alarm "
                      "should not be ringing: FAIL");
         }
      }
   }

   // Test a minute-based alarm and periodic 20-minute alarm using a start
   // time and second interval

   TimeInstant Time30min(2019, 8, 15, 14, 55, 23.25);
   Alarm Alarm30min("30min", Time30min);

   TimeInterval Interval20min(20, TimeUnits::Minutes);
   Alarm AlarmEvery20min("Every 20 minutes", Interval20min, Time0);
   TimeInterval IntervalSecond(1, TimeUnits::Seconds);

   CurTime = StartTime; // start time is 2019-08-15_14:25:23.25
   Err1    = AlarmEvery20min.reset(CurTime); // ensure first alarm in future

   Err3 = 1;                          // to limit output
   for (int N = 1; N <= 10800; ++N) { // integrate for 3 hours
      // increment time in second intervals
      CurTime += IntervalSecond;

      // update alarm state based on current time
      Err1 = Alarm30min.updateStatus(CurTime);
      Err2 = AlarmEvery20min.updateStatus(CurTime);
      if (Err1 != 0 && Err2 != 0) {
         ++ErrAll;
         LOG_ERROR("TimeMgrTest/Alarm: update minute alarms: FAIL");
      }

      // Test whether one-time alarm should be ringing or not
      if (N == 1800) {
         if (Alarm30min.isRinging()) {
            LOG_INFO("TimeMgrTest/Alarm: one-time minute alarm: PASS");
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time minute alarm: FAIL");
         }
         Err1 = Alarm30min.stop();
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time minute alarm stop: FAIL");
         }
      } else {
         if (Alarm30min.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: one-time minute alarm "
                      "should not be ringing: FAIL");
         }
      }

      // Test whether interval alarm should be ringing or not
      if ((N - 877) % 1200 == 0) {
         if (AlarmEvery20min.isRinging()) {
            if (Err3 == 1) // only print first instance of pass
               LOG_INFO("TimeMgrTest/Alarm: periodic minute alarm: PASS");
            Err3 = 0;
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic minute alarm: FAIL");
         }
         Err1 = AlarmEvery20min.reset(CurTime);
         if (Err1 != 0) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic minute alarm reset: FAIL");
         }
      } else {
         if (AlarmEvery20min.isRinging()) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Alarm: periodic minute alarm "
                      "should not be ringing: FAIL");
         }
      }
   }

   return ErrAll;

} // end testAlarm

//------------------------------------------------------------------------------
// Clock test

int testClock(void) {

   LOG_INFO("TimeMgrTest: Clock tests ---------------------------------------");

   // Initialize error codes
   I4 Err1{0};
   I4 Err2{0};
   I4 Err3{0};
   I4 ErrAll{0};

   // For various time intervals, we create alarms at relevant times
   // and step through time to trigger the alarm.

   // Do all testing in Gregorian calendar
   Calendar::reset();
   Calendar::init("Gregorian");

   // Create an initial model clock with a 2000 start time and
   // one hour time step.

   TimeInstant Time0(2000, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);

   Clock ModelClock(Time0, TimeStep);

   // Test some basic retrieval functions

   TimeInstant TimeCheck; // init to default values
   TimeInterval StepCheck;

   TimeCheck = ModelClock.getStartTime();
   StepCheck = ModelClock.getTimeStep();

   if (TimeCheck == Time0) {
      LOG_INFO("TimeMgrTest/Clock: get start time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: get start time: FAIL");
   }

   if (StepCheck == TimeStep) {
      LOG_INFO("TimeMgrTest/Clock: get time step: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: get time step: FAIL");
   }

   // Define a number of periodic and one-time alarms

   TimeInstant TimeNewYear2020(2020, 1, 1, 0, 0, 0.0);
   TimeInstant Time2020Mar1(2020, 3, 1, 0, 0, 0.0);
   TimeInstant Time2019Aug20(2019, 8, 20, 0, 0, 0.0);

   Alarm AlarmNewYear2020("New Year 2020", TimeNewYear2020);
   Alarm Alarm2020Mar1("2020-03-01", Time2020Mar1);
   Alarm Alarm2019Aug20("2020-08-20", Time2019Aug20);

   TimeInterval IntervalAnnual(1, TimeUnits::Years);
   TimeInterval IntervalMonthly(1, TimeUnits::Months);
   TimeInterval IntervalDaily(1, TimeUnits::Days);
   TimeInterval IntervalHourly(1, TimeUnits::Hours);
   TimeInterval Interval6Hour(6, TimeUnits::Hours);
   TimeInterval Interval20min(20, TimeUnits::Minutes);

   Alarm AlarmEveryYear("Every Year", IntervalAnnual, Time0);
   Alarm AlarmEveryMonth("Every Month", IntervalMonthly, Time0);
   Alarm AlarmEveryDay("Every Day", IntervalDaily, Time0);
   Alarm AlarmEveryHour("Every Hour", IntervalHourly, Time0);
   Alarm AlarmEvery6Hour("Every 6 Hours", Interval6Hour, Time0);
   Alarm AlarmEvery20min("Every 20 minutes", Interval20min, Time0);

   // Test adding alarms to clock

   Err1 = ModelClock.attachAlarm(&AlarmNewYear2020);
   Err2 = ModelClock.attachAlarm(&Alarm2020Mar1);
   Err3 = ModelClock.attachAlarm(&Alarm2019Aug20);

   if (Err1 == 0 && Err2 == 0 && Err3 == 0) {
      LOG_INFO("TimeMgrTest/Clock: attach one-time alarms: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: attach one-time alarms: FAIL");
   }

   Err1 = ModelClock.attachAlarm(&AlarmEveryYear);
   Err2 = ModelClock.attachAlarm(&AlarmEveryMonth);
   Err3 = ModelClock.attachAlarm(&AlarmEveryDay);

   if (Err1 == 0 && Err2 == 0 && Err3 == 0) {
      LOG_INFO("TimeMgrTest/Clock: attach periodic alarms 1: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: attach periodic alarms 1: FAIL");
   }

   Err1 = ModelClock.attachAlarm(&AlarmEveryHour);
   Err2 = ModelClock.attachAlarm(&AlarmEvery6Hour);
   Err3 = ModelClock.attachAlarm(&AlarmEvery20min);

   if (Err1 == 0 && Err2 == 0 && Err3 == 0) {
      LOG_INFO("TimeMgrTest/Clock: attach periodic alarms 2: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: attach periodic alarms 2: FAIL");
   }

   // Test changing the time step

   Err1 = ModelClock.changeTimeStep(Interval20min);

   StepCheck = ModelClock.getTimeStep();
   if (StepCheck == Interval20min && Err1 == 0) {
      LOG_INFO("TimeMgrTest/Clock: change time step: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: change time step: FAIL");
   }

   // Test setting new current time and retrieving current, previous,
   // and next times.

   TimeInstant CurrTime(2019, 1, 1, 0, 0, 0.0);
   TimeInstant PrevTime(2018, 12, 31, 23, 40, 0.0);
   TimeInstant NextTime(2019, 1, 1, 0, 20, 0.0);

   Err1      = ModelClock.setCurrentTime(CurrTime);
   TimeCheck = ModelClock.getCurrentTime();

   if (TimeCheck == CurrTime && Err1 == 0) {
      LOG_INFO("TimeMgrTest/Clock: set/get current time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: set/get current time: FAIL");
   }

   TimeCheck = ModelClock.getPreviousTime();
   if (TimeCheck == PrevTime) {
      LOG_INFO("TimeMgrTest/Clock: set/get previous time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: set/get previous time: FAIL");
   }

   TimeCheck = ModelClock.getNextTime();
   if (TimeCheck == NextTime) {
      LOG_INFO("TimeMgrTest/Clock: set/get next time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Clock: set/get next time: FAIL");
   }

   // Update periodic alarms to new current time

   Err1 = AlarmEveryYear.reset(CurrTime);
   Err1 = AlarmEveryMonth.reset(CurrTime);
   Err1 = AlarmEveryDay.reset(CurrTime);
   Err1 = AlarmEveryHour.reset(CurrTime);
   Err1 = AlarmEvery6Hour.reset(CurrTime);
   Err1 = AlarmEvery20min.reset(CurrTime);

   // Test the integration of a model clock by advancing forward
   // in time and checking alarms. Integrate forward 2 years with a
   // 20 min timestep.

   TimeInstant StopTime(2021, 1, 1, 0, 0, 0.0);
   bool FirstStep{true};
   bool RingCheck{false};
   I8 Year{0};
   I8 Month{0};
   I8 Day{0};
   I8 Hour{0};
   I8 Minute{0};
   R8 Second{0.0};

   while (CurrTime <= StopTime) {

      Err1 = ModelClock.advance(); // advance one time step

      if (Err1 == 0) {
         if (FirstStep)
            LOG_INFO("TimeMgrTest/Clock: advance: PASS");
      } else {
         ++ErrAll;
         break;
         LOG_ERROR("TimeMgrTest/Clock: advance: FAIL");
      }

      // retrieve current time for both loop cycling and tests below
      CurrTime = ModelClock.getCurrentTime();

      // check various one-time alarms and stop if needed
      RingCheck = AlarmNewYear2020.isRinging();
      if (CurrTime == TimeNewYear2020) {
         if (RingCheck) {
            LOG_INFO("TimeMgrTest/Clock: alarm NewYear2020: PASS");
            Err1 = AlarmNewYear2020.stop();
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm NewYear2020: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm NewYear2020: FAIL");
         }
      }

      RingCheck = Alarm2020Mar1.isRinging();
      if (CurrTime == Time2020Mar1) {
         if (RingCheck) {
            LOG_INFO("TimeMgrTest/Clock: alarm 2020Mar1: PASS");
            Err1 = Alarm2020Mar1.stop();
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm 2020Mar1: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm 2020Mar1: FAIL");
         }
      }

      RingCheck = Alarm2019Aug20.isRinging();
      if (CurrTime == Time2019Aug20) {
         if (RingCheck) {
            LOG_INFO("TimeMgrTest/Clock: alarm 2019Aug20: PASS");
            Err1 = Alarm2019Aug20.stop();
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm 2019Aug20: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm 2019Aug20: FAIL");
         }
      }

      // check 20-min alarm should always be ringing
      RingCheck = AlarmEvery20min.isRinging();
      if (RingCheck) {
         if (FirstStep)
            LOG_INFO("TimeMgrTest/Clock: alarm Every20min: PASS");
         Err1 = AlarmEvery20min.reset(CurrTime);
      } else {
         ++ErrAll;
         LOG_INFO("TimeMgrTest/Clock: alarm Every20min: FAIL");
      }

      // Extract year, month, day, hour, min, seconds from current
      // time to check other periodic alarms

      Err1 = CurrTime.get(Year, Month, Day, Hour, Minute, Second);

      // Check annual alarm
      RingCheck = AlarmEveryYear.isRinging();
      if (Month == 1 && Day == 1 && Hour == 0 && Minute == 0 && Second == 0.0) {
         if (RingCheck) {
            // success but avoid printing excessive PASS output
            Err1 = AlarmEveryYear.reset(CurrTime);
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryYear 1: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryYear 2: FAIL");
         }
      }

      // Check monthly alarm
      RingCheck = AlarmEveryMonth.isRinging();
      if (Day == 1 && Hour == 0 && Minute == 0 && Second == 0.0) {
         if (RingCheck) {
            // success but avoid printing excessive PASS output
            Err1 = AlarmEveryMonth.reset(CurrTime);
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryMonth 1: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryMonth 2: FAIL");
         }
      }

      // Check daily alarm
      RingCheck = AlarmEveryDay.isRinging();
      if (Hour == 0 && Minute == 0 && Second == 0.0) {
         if (RingCheck) {
            // success but avoid printing excessive PASS output
            Err1 = AlarmEveryDay.reset(CurrTime);
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryDay 1: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryDay 2: FAIL");
         }
      }

      // Check hourly alarm
      RingCheck = AlarmEveryHour.isRinging();
      if (Minute == 0 && Second == 0.0) {
         if (RingCheck) {
            // success but avoid printing excessive PASS output
            Err1 = AlarmEveryHour.reset(CurrTime);
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryHour 1: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm EveryHour 2: FAIL");
         }
      }

      // Check 6-hour alarm
      RingCheck = AlarmEvery6Hour.isRinging();
      if (Hour % 6 == 0 && Minute == 0 && Second == 0.0) {
         if (RingCheck) {
            // success but avoid printing excessive PASS output
            Err1 = AlarmEvery6Hour.reset(CurrTime);
         } else {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm Every6Hour 1: FAIL");
         }
      } else {
         if (RingCheck) {
            ++ErrAll;
            LOG_ERROR("TimeMgrTest/Clock: alarm Every6Hour 2: FAIL");
         }
      }

      // reset flag for stuff on first step
      if (FirstStep)
         FirstStep = false;
   }

   return ErrAll;

} // end testClock

//------------------------------------------------------------------------------
// The test driver.

int main(int argc, char *argv[]) {

   I4 Err{0};
   I4 TotErr{0};

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);

   // Initialize the Machine Environment and retrieve the default environment
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();

   // Initialize the Logging system
   initLogging(DefEnv);

   Err = testTimeFrac();
   TotErr += Err;

   Err = testCalendar();
   TotErr += Err;

   Err = testTimeInterval();
   TotErr += Err;

   Err = testTimeInstant();
   TotErr += Err;

   Err = testAlarm();
   TotErr += Err;

   Err = testClock();
   TotErr += Err;

   MPI_Finalize();

   if (TotErr == 0) {
      LOG_INFO("TimeMgrTest: Successful completion");
   } else {
      LOG_INFO("TimeMgrTest: Failed");
   }

   if (TotErr >= 256)
      TotErr = 255;

   return TotErr;

} // end of main
//===-----------------------------------------------------------------------===/
