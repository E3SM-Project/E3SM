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
// Calendar test

int testCalendar(void) {

   LOG_INFO("TimeMgrTest: Calendar tests ------------------------------------");

   // Initialize error codes
   OMEGA::I4 Err1{0};
   OMEGA::I4 Err2{0};
   OMEGA::I4 ErrAll{0};

   // Test default constructor.
   // Also tests the get routine.

   OMEGA::Calendar CalEmpty;

   OMEGA::CalendarKind Kind0 = OMEGA::CalendarNoCalendar;
   OMEGA::I4 ID0             = 1;
   std::string Name0(" ");
   OMEGA::I4 DaysPerMonth0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   OMEGA::I4 MonthsPerYear0  = 12;
   OMEGA::I4 SecondsPerDay0  = 0;
   OMEGA::I4 SecondsPerYear0 = 0;
   OMEGA::I4 DaysPerYear0    = 0;

   OMEGA::CalendarKind Kind1 = OMEGA::CalendarGregorian;
   OMEGA::I4 ID1             = 1;
   std::string Name1("junk");
   OMEGA::I4 DaysPerMonth1[] = {99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99};
   OMEGA::I4 MonthsPerYear1  = 999;
   OMEGA::I4 SecondsPerDay1  = 999;
   OMEGA::I4 SecondsPerYear1 = 999;
   OMEGA::I4 DaysPerYear1    = 999;

   Err1 = CalEmpty.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                       &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Err1 != 0 || Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err2 = 1;

   for (int I = 0; I < MonthsPerYear0; I++) {
      if (DaysPerMonth0[I] != DaysPerMonth1[I])
         Err2 = 1;
   }

   if (Err2 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: default constructor and get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: default constructor or get: FAIL");
   }

   // Test rename function
   Err1  = CalEmpty.rename("junk");
   Name1 = "junk";

   Err2 = CalEmpty.get(nullptr, &Name0, nullptr, nullptr, nullptr, nullptr,
                       nullptr, nullptr);

   if (Err1 == 0 && Err2 == 0 && Name1 == Name0) {
      LOG_INFO("TimeMgrTest/Calendar: rename: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: rename: FAIL");
   }

   // Test custom calendar

   Kind0             = OMEGA::CalendarCustom;
   ID0               = 2;
   Name0             = "Custom";
   DaysPerMonth0[0]  = 10;
   DaysPerMonth0[1]  = 10;
   DaysPerMonth0[2]  = 10;
   DaysPerMonth0[3]  = 10;
   DaysPerMonth0[4]  = 10;
   DaysPerMonth0[5]  = 10;
   DaysPerMonth0[6]  = 10;
   DaysPerMonth0[7]  = 10;
   DaysPerMonth0[8]  = 10;
   DaysPerMonth0[9]  = 10;
   DaysPerMonth0[10] = 10;
   DaysPerMonth0[11] = 14;
   MonthsPerYear0    = 12;
   SecondsPerDay0    = 100;
   SecondsPerYear0   = 12400;
   DaysPerYear0      = 124;

   Kind1             = OMEGA::CalendarCustom;
   ID1               = 0;
   Name1             = "junk";
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

   OMEGA::Calendar CalCustom(Name0, DaysPerMonth0, SecondsPerDay0,
                             SecondsPerYear0, DaysPerYear0);

   Err1 = CalCustom.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                        &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Err1 != 0 || Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0) {
      Err2 = 1;
   } else {
      Err2 = 0;
   }

   if (Err2 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: custom constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: custom constructor: FAIL");
   }

   // Test copy constructor

   ID0 = 3;

   OMEGA::Calendar CalCopy(CalCustom);

   Err1 = CalCopy.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                      &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Err1 != 0 || Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err2 = 1;
   else
      Err2 = 0;

   if (Err2 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: copy constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: copy constructor: FAIL");
   }

   // Test equivalence and non-equivalence

   if (CalCustom == CalCopy) {
      LOG_INFO("TimeMgrTest/Calendar: operator(==): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: operator(==): FAIL");
   }

   if (CalCustom != CalEmpty) {
      LOG_INFO("TimeMgrTest/Calendar: operator(!=): PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: operator(!=): FAIL");
   }

   // Test calendar construction for Gregorian

   Kind0             = OMEGA::CalendarGregorian;
   ID0               = 4;
   Name0             = "Gregorian";
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

   Kind1             = OMEGA::CalendarNoCalendar;
   ID1               = 0;
   Name1             = "junk";
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

   OMEGA::Calendar CalGreg("Gregorian", OMEGA::CalendarGregorian);

   Err1 = CalGreg.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                      &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Gregorian constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Gregorian constructor: FAIL");
   }

   // No leap calendar is identical

   Kind1             = OMEGA::CalendarNoCalendar;
   ID1               = 0;
   Name1             = "junk";
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

   ID0   = 5;
   Kind0 = OMEGA::CalendarNoLeap;
   Name0 = "Noleap";

   OMEGA::Calendar CalNoLeap("Noleap", OMEGA::CalendarNoLeap);

   Err1 = CalNoLeap.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                        &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: No leap constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: No leap constructor: FAIL");
   }

   // Straight Julian calendar the same

   Kind1             = OMEGA::CalendarNoCalendar;
   ID1               = 0;
   Name1             = "junk";
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

   ID0   = 6;
   Kind0 = OMEGA::CalendarJulian;
   Name0 = "Julian";

   OMEGA::Calendar CalJulian("Julian", OMEGA::CalendarJulian);

   Err1 = CalJulian.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                        &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Julian constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Julian constructor: FAIL");
   }

   // Test calendar construction for 360 day

   Kind0             = OMEGA::Calendar360Day;
   ID0               = 7;
   Name0             = "360Day";
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

   Kind1             = OMEGA::CalendarNoCalendar;
   ID1               = 0;
   Name1             = "junk";
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

   OMEGA::Calendar Cal360Day("360Day", OMEGA::Calendar360Day);

   Err1 = Cal360Day.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                        &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: 360 Day constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: 360 Day constructor: FAIL");
   }

   // Test calendar construction for Julian day

   Kind0             = OMEGA::CalendarJulianDay;
   ID0               = 8;
   Name0             = "JulianDay";
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

   Kind1             = OMEGA::CalendarNoCalendar;
   ID1               = 0;
   Name1             = "junk";
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

   OMEGA::Calendar CalJulianDay("JulianDay", OMEGA::CalendarJulianDay);

   Err1 = CalJulianDay.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                           &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Julian Day constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Julian Day constructor: FAIL");
   }

   // Modified Julian day identical

   Kind0 = OMEGA::CalendarModJulianDay;
   ID0   = 9;
   Name0 = "ModJulianDay";

   Kind1             = OMEGA::CalendarNoCalendar;
   ID1               = 0;
   Name1             = "junk";
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

   OMEGA::Calendar CalModJulianDay("ModJulianDay", OMEGA::CalendarModJulianDay);

   Err1 =
       CalModJulianDay.get(&ID1, &Name1, &Kind1, DaysPerMonth1, &MonthsPerYear1,
                           &SecondsPerDay1, &SecondsPerYear1, &DaysPerYear1);

   if (Kind1 != Kind0 || ID1 != ID0 || Name1 != Name0 ||
       MonthsPerYear1 != MonthsPerYear0 || SecondsPerDay1 != SecondsPerDay0 ||
       SecondsPerYear1 != SecondsPerYear0 || DaysPerYear1 != DaysPerYear0)
      Err1 = 1;

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: Modified Julian Day constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: Modified Julian Day constructor: FAIL");
   }

   // Test leap year functions for different calendars

   OMEGA::I8 TstYear = 1981;

   // Calendar with no leap year
   if (!CalJulianDay.isLeapYear(TstYear, Err1)) {
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

   // Leap year calendar but not a leap year
   if (!CalGreg.isLeapYear(TstYear, Err1)) {
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

   // Leap year calendar normal leap years
   TstYear = 1984;
   if (CalGreg.isLeapYear(TstYear, Err1)) {
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
   if (CalJulian.isLeapYear(TstYear, Err1)) {
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
   if (!CalNoLeap.isLeapYear(TstYear, Err1)) {
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

   // Special Gregorian leap year exceptions
   TstYear = 1900;
   if (!CalGreg.isLeapYear(TstYear, Err1)) {
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
   if (CalGreg.isLeapYear(TstYear, Err1)) {
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
   OMEGA::TimeFrac TstTime(210480449184, 1, 4);
   TstYear = 1957;
   OMEGA::I8 TstMonth{10};
   OMEGA::I8 TstDay{4};
   OMEGA::I8 TstHour{19};
   OMEGA::I8 TstMinute{26};
   OMEGA::I8 TstSecondW{24};
   OMEGA::I8 TstSecondN{1};
   OMEGA::I8 TstSecondD{4};

   OMEGA::TimeFrac ChkTime(0, 0, 1);
   OMEGA::I8 ChkYear{0};
   OMEGA::I8 ChkMonth{0};
   OMEGA::I8 ChkDay{0};
   OMEGA::I8 ChkHour{0};
   OMEGA::I8 ChkMinute{0};
   OMEGA::I8 ChkSecondW{0};
   OMEGA::I8 ChkSecondN{0};
   OMEGA::I8 ChkSecondD{1};

   ChkTime =
       CalGreg.getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                              TstSecondW, TstSecondN, TstSecondD);

   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert Gregorian date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert Gregorian date to "
                "elapsed time: FAIL");
   }

   Err1 = CalGreg.getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
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

   // Use same date for no leap calendar
   OMEGA::I8 TmpSeconds =
       (1957 * (OMEGA::I8)86400 * 365) +
       (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 4 - 1) * (OMEGA::I8)86400 +
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
       CalNoLeap.getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert NoLeap date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert NoLeap date to "
                "elapsed time: FAIL");
   }

   Err1 = CalNoLeap.getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
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

   // Use same date for 360-day
   TmpSeconds = (1957 * (OMEGA::I8)86400 * 360) +
                (9 * 30 + 4 - 1) * (OMEGA::I8)86400 + 19 * 3600 + 26 * 60 + 24;
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
       Cal360Day.getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert 360-day date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert 360-day date to "
                "elapsed time: FAIL");
   }

   Err1 = Cal360Day.getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
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

   // Modify same date to fit in custom calendar (1957-10-4, 00:01:24.25)
   TmpSeconds = (1957 * (OMEGA::I8)12400) + (9 * 10 + 4 - 1) * 100 + 60 + 24;
   Err1       = TstTime.set(TmpSeconds, 1, 4);
   TstYear    = 1957;
   TstMonth   = 10;
   TstDay     = 4;
   TstHour    = 0;
   TstMinute  = 1;
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
       CalCustom.getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert custom date to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert custom date to "
                "elapsed time: FAIL");
   }

   Err1 = CalCustom.getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
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

   // Julian Day - test Julian Day of 2436116.31 (plus 1/4 second)
   TmpSeconds = (2436116 - 1) * (OMEGA::I8)86400 + 7 * 3600 + 26 * 60 + 24;
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

   ChkTime = CalJulianDay.getElapsedTime(TstYear, TstMonth, TstDay, TstHour,
                                         TstMinute, TstSecondW, TstSecondN,
                                         TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert Julian day to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert Julian day to "
                "elapsed time: FAIL");
   }

   Err1 =
       CalJulianDay.getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
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

   // Mod Julian Day - test Mod Julian Day of 2436116.31 (plus 1/4 second)
   TmpSeconds = (2436116 - 1) * (OMEGA::I8)86400 + 7 * 3600 + 26 * 60 + 24;
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

   ChkTime = CalModJulianDay.getElapsedTime(TstYear, TstMonth, TstDay, TstHour,
                                            TstMinute, TstSecondW, TstSecondN,
                                            TstSecondD);
   if (ChkTime == TstTime) {
      LOG_INFO("TimeMgrTest/Calendar: convert Mod Julian day to "
               "elapsed time: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: convert Mod Julian day to "
                "elapsed time: FAIL");
   }

   Err1 = CalModJulianDay.getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay,
                                      ChkHour, ChkMinute, ChkSecondW,
                                      ChkSecondN, ChkSecondD);
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
       CalJulian.getElapsedTime(TstYear, TstMonth, TstDay, TstHour, TstMinute,
                                TstSecondW, TstSecondN, TstSecondD);

   Err1 = CalJulian.getDateTime(ChkTime, ChkYear, ChkMonth, ChkDay, ChkHour,
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

   // Test calendar date increment function

   // Check normal year increment/decrement

   TstYear  = 1983;
   TstMonth = 6;
   TstDay   = 15;
   ChkYear  = 1985;
   ChkMonth = 6;
   ChkDay   = 15;

   Err1 = CalGreg.incrementDate(2, OMEGA::TimeUnits::Years, TstYear, TstMonth,
                                TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by year: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date by year: FAIL");
   }

   TstYear = 1983;
   Err1 = CalNoLeap.incrementDate(2, OMEGA::TimeUnits::Years, TstYear, TstMonth,
                                  TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by year: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by year: FAIL");
   }

   TstYear = 1983;
   ChkYear = 1981;
   Err1    = CalNoLeap.incrementDate(-2, OMEGA::TimeUnits::Years, TstYear,
                                     TstMonth, TstDay);
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

   Err1 = CalGreg.incrementDate(2, OMEGA::TimeUnits::Months, TstYear, TstMonth,
                                TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by month: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date "
                "by month: FAIL");
   }

   TstMonth = 6;
   Err1     = CalNoLeap.incrementDate(2, OMEGA::TimeUnits::Months, TstYear,
                                      TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by month: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by month: FAIL");
   }

   TstMonth = 6;
   ChkMonth = 4;
   Err1     = CalNoLeap.incrementDate(-2, OMEGA::TimeUnits::Months, TstYear,
                                      TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement NoLeap date by month: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement NoLeap date by month: FAIL");
   }

   // Check year rollover for longer month intervals

   TstYear  = 1984;
   TstMonth = 10;
   TstDay   = 15;
   ChkYear  = 1986;
   ChkMonth = 4;
   ChkDay   = 15;

   Err1 = CalGreg.incrementDate(18, OMEGA::TimeUnits::Months, TstYear, TstMonth,
                                TstDay);
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

   Err1 = CalGreg.incrementDate(-18, OMEGA::TimeUnits::Months, TstYear,
                                TstMonth, TstDay);
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

   Err1 = CalGreg.incrementDate(1, OMEGA::TimeUnits::Months, TstYear, TstMonth,
                                TstDay);
   if (Err1 != 0) {
      LOG_INFO("TimeMgrTest/Calendar: increment catch bad day range: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment catch bad day range: FAIL");
   }

   // Test normal daily increments/decrements (include a leap day for Gregorian)

   TstYear  = 1984;
   TstMonth = 2;
   TstDay   = 25;
   ChkYear  = 1984;
   ChkMonth = 3;
   ChkDay   = 6;

   Err1 = CalGreg.incrementDate(10, OMEGA::TimeUnits::Days, TstYear, TstMonth,
                                TstDay);
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
   TstMonth = 2;
   TstDay   = 25;
   ChkYear  = 1984;
   ChkMonth = 3;
   ChkDay   = 7;

   Err1 = CalNoLeap.incrementDate(10, OMEGA::TimeUnits::Days, TstYear, TstMonth,
                                  TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by 10 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by 10 days: FAIL");
   }

   TstYear  = 1984;
   TstMonth = 3;
   TstDay   = 6;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 = CalGreg.incrementDate(-10, OMEGA::TimeUnits::Days, TstYear, TstMonth,
                                TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement Gregorian date by "
               "10 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement Gregorian date by "
                "10 days: FAIL");
   }

   TstYear  = 1984;
   TstMonth = 3;
   TstDay   = 7;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 = CalNoLeap.incrementDate(-10, OMEGA::TimeUnits::Days, TstYear,
                                  TstMonth, TstDay);
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
   ChkMonth = 3;
   ChkDay   = 31;

   Err1 = CalGreg.incrementDate(400, OMEGA::TimeUnits::Days, TstYear, TstMonth,
                                TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment Gregorian date by "
               "400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment Gregorian date by "
                "400 days: FAIL");
   }

   TstYear  = 1984;
   TstMonth = 2;
   TstDay   = 25;
   ChkYear  = 1985;
   ChkMonth = 4;
   ChkDay   = 1;

   Err1 = CalNoLeap.incrementDate(400, OMEGA::TimeUnits::Days, TstYear,
                                  TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: increment NoLeap date by 400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: increment NoLeap date by "
                "400 days: FAIL");
   }

   TstYear  = 1985;
   TstMonth = 3;
   TstDay   = 31;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 = CalGreg.incrementDate(-400, OMEGA::TimeUnits::Days, TstYear, TstMonth,
                                TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement Gregorian date by "
               "400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement Gregorian date by "
                "400 days: FAIL");
   }

   TstYear  = 1985;
   TstMonth = 4;
   TstDay   = 1;
   ChkYear  = 1984;
   ChkMonth = 2;
   ChkDay   = 25;

   Err1 = CalNoLeap.incrementDate(-400, OMEGA::TimeUnits::Days, TstYear,
                                  TstMonth, TstDay);
   if (Err1 == 0 && TstYear == ChkYear && TstMonth == ChkMonth &&
       TstDay == ChkDay) {
      LOG_INFO("TimeMgrTest/Calendar: decrement NoLeap date by 400 days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: decrement NoLeap date by "
                "400 days: FAIL");
   }

   // Test validate

   Err1 = CalCustom.validate();

   if (Err1 == 0) {
      LOG_INFO("TimeMgrTest/Calendar: validate: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/Calendar: validate: FAIL");
   }

   return ErrAll;

} // end testCalendar

//------------------------------------------------------------------------------
// TimeInterval test

int testTimeInterval(void) {

   LOG_INFO("TimeMgrTest: TimeInterval tests --------------------------------");

   // Initialize error codes
   OMEGA::I4 Err1{0};
   OMEGA::I4 Err2{0};
   OMEGA::I4 ErrAll{0};

   // Initialize some reference values for the fractional
   // representation of a TimeInterval 5 3/8 seconds.
   OMEGA::I8 WRef{5};
   OMEGA::I8 NRef{3};
   OMEGA::I8 DRef{8};
   OMEGA::R8 RRef{5.375};
   OMEGA::I8 IRef{WRef};

   OMEGA::I8 WTst{5};
   OMEGA::I8 NTst{3};
   OMEGA::I8 DTst{8};
   OMEGA::R8 RTst{5.375};
   OMEGA::I8 ITst{WTst};

   // Test default constructor to create a reference fraction
   // Also implicitly tests one form of the get routine.

   OMEGA::TimeInterval TiRef;

   // Test get function for a fractional second representation
   Err1 = TiRef.get(WTst, NTst, DTst);

   if (Err1 == 0 && WTst == 0 && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: default constructor and get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: default constructor and get: FAIL");
   }

   // Test constructor from fractional seconds

   OMEGA::TimeInterval TiTstFS(WRef, NRef, DRef);

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
   OMEGA::TimeInterval TiTstRSec(RRef, OMEGA::TimeUnits::Seconds);

   if (TiTstRSec == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real seconds constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real seconds constructor: FAIL");
   }

   // Test time interval constructor from real minutes
   OMEGA::TimeInterval TiTstRMin((RRef / SECONDS_PER_MINUTE),
                                 OMEGA::TimeUnits::Minutes);
   if (TiTstRMin == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real minutes constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real minutes constructor: FAIL");
   }

   // Test time interval constructor from real hours
   OMEGA::TimeInterval TiTstRHour((RRef / SECONDS_PER_HOUR),
                                  OMEGA::TimeUnits::Hours);
   if (TiTstRHour == TiRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real hours constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real hours constructor: FAIL");
   }

   // Test time interval constructor from real days
   OMEGA::TimeInterval TiTstRDay(RRef, OMEGA::TimeUnits::Days);

   Err1 = TiTstRDay.get(RTst, OMEGA::TimeUnits::Days);
   if (Err1 == 0 && RTst == 5.0) {
      LOG_INFO("TimeMgrTest/TimeInterval: real days constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real days constructor: FAIL");
   }

   // Test time interval constructor from real months
   OMEGA::TimeInterval TiTstRMonth(RRef, OMEGA::TimeUnits::Months);

   Err1 = TiTstRMonth.get(RTst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && RTst == 5.0) {
      LOG_INFO("TimeMgrTest/TimeInterval: real months constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real months constructor: FAIL");
   }

   // Test time interval constructor from real years
   OMEGA::TimeInterval TiTstRYear(RRef, OMEGA::TimeUnits::Years);

   Err1 = TiTstRYear.get(RTst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && RTst == 5.0) {
      LOG_INFO("TimeMgrTest/TimeInterval: real years constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real years constructor: FAIL");
   }

   // Test time interval constructor from integer seconds
   OMEGA::TimeInterval TiTstISec(IRef, OMEGA::TimeUnits::Seconds);

   Err1 = TiTstISec.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == WRef && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer seconds constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer seconds constructor: FAIL");
   }

   // Test time interval constructor from integer minutes
   OMEGA::TimeInterval TiTstIMin(IRef, OMEGA::TimeUnits::Minutes);

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
   OMEGA::TimeInterval TiTstIHour(IRef, OMEGA::TimeUnits::Hours);

   Err1 = TiTstIHour.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == WRef * SECONDS_PER_HOUR && NTst == 0 && DTst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer hours constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer hours constructor: FAIL");
   }

   // Test calendar-based time interval constructor in years
   OMEGA::TimeInterval TiTstIYear(IRef, OMEGA::TimeUnits::Years);

   Err1 = TiTstIYear.get(ITst, OMEGA::TimeUnits::Years);

   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: year constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: year constructor: FAIL");
   }

   // Test calendar-based time interval constructor in months
   OMEGA::TimeInterval TiTstIMonth(IRef, OMEGA::TimeUnits::Months);

   Err1 = TiTstIMonth.get(ITst, OMEGA::TimeUnits::Months);

   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: month constructor: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: month constructor: FAIL");
   }

   // Test calendar-based time interval constructor in days
   OMEGA::TimeInterval TiTstIDay(IRef, OMEGA::TimeUnits::Days);

   Err1 = TiTstIDay.get(ITst, OMEGA::TimeUnits::Days);

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

   Err1 = TiTstRSec.set(RRef, OMEGA::TimeUnits::Seconds);
   Err2 = TiTstRSec.get(RTst, OMEGA::TimeUnits::Seconds);

   if (Err1 == 0 && Err2 == 0 && RTst == RRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real seconds put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real seconds put/get: FAIL");
   }

   // Test integer seconds put/get

   Err1 = TiTstISec.set(IRef, OMEGA::TimeUnits::Seconds);
   Err2 = TiTstISec.get(ITst, OMEGA::TimeUnits::Seconds);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer seconds put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer seconds put/get: FAIL");
   }

   // Test real hours put/get

   Err1 = TiTstRSec.set(RRef, OMEGA::TimeUnits::Hours);
   Err2 = TiTstRSec.get(RTst, OMEGA::TimeUnits::Hours);

   if (Err1 == 0 && Err2 == 0 && RTst == RRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real hours put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real hours put/get: FAIL");
   }

   // Test integer hours put/get

   Err1 = TiTstISec.set(IRef, OMEGA::TimeUnits::Hours);
   Err2 = TiTstISec.get(ITst, OMEGA::TimeUnits::Hours);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer hours put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer hours put/get: FAIL");
   }

   // Test real minutes put/get

   Err1 = TiTstRSec.set(RRef, OMEGA::TimeUnits::Minutes);
   Err2 = TiTstRSec.get(RTst, OMEGA::TimeUnits::Minutes);

   if (Err1 == 0 && Err2 == 0 && RTst == RRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: real minutes put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real minutes put/get: FAIL");
   }

   // Test integer minutes put/get

   Err1 = TiTstISec.set(IRef, OMEGA::TimeUnits::Minutes);
   Err2 = TiTstISec.get(ITst, OMEGA::TimeUnits::Minutes);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer minutes put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer minutes put/get: FAIL");
   }

   // Test integer years put/get

   Err1 = TiTstISec.set(IRef + 1, OMEGA::TimeUnits::Years);
   Err2 = TiTstISec.get(ITst, OMEGA::TimeUnits::Years);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef + 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer years put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer years put/get: FAIL");
   }

   // Test integer months put/get

   Err1 = TiTstISec.set(IRef + 1, OMEGA::TimeUnits::Months);
   Err2 = TiTstISec.get(ITst, OMEGA::TimeUnits::Months);

   if (Err1 == 0 && Err2 == 0 && ITst == IRef + 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: integer months put/get: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: integer months put/get: FAIL");
   }

   // Test integer days put/get

   Err1 = TiTstISec.set(IRef + 1, OMEGA::TimeUnits::Days);
   Err2 = TiTstISec.get(ITst, OMEGA::TimeUnits::Days);

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

   Err1 = TiTstIYear.set(IRef, OMEGA::TimeUnits::Years);
   Err1 = TiTstIMonth.set(IRef, OMEGA::TimeUnits::Months);
   Err1 = TiTstIDay.set(IRef, OMEGA::TimeUnits::Days);

   OMEGA::TimeInterval TiTstIYear2(IRef - 1, OMEGA::TimeUnits::Years);
   OMEGA::TimeInterval TiTstIMonth2(IRef - 1, OMEGA::TimeUnits::Months);
   OMEGA::TimeInterval TiTstIDay2(IRef - 1, OMEGA::TimeUnits::Days);

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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: addition years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: addition years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: addition months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: addition months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == (IRef + IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: increment years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: increment years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == (IRef + IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: increment months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: increment months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: subtraction years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: subtraction years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == (IRef + IRef)) {
      LOG_INFO("TimeMgrTest/TimeInterval: subtraction months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: subtraction months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: decrement years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: decrement years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: decrement months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: decrement months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm int multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm int multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 3 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm int multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm int multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 9 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply in place years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply in place years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 9 * IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: int multiply in place months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int multiply in place months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm real multiply years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm real multiply years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 16) {
      LOG_INFO("TimeMgrTest/TimeInterval: comm real multiply months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: comm real multiply months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 52) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply in place years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply in place years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 52) {
      LOG_INFO("TimeMgrTest/TimeInterval: real multiply in place months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: real multiply in place "
                "months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 1) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == 0) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide in place years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide in place years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == 0) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide in place months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide in place months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
   if (Err1 == 0 && ITst == 0) {
      LOG_INFO("TimeMgrTest/TimeInterval: int divide in place days: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: int divide in place days: FAIL");
   }

   // Test negative absolute value first to change sign

   TiTstFS      = OMEGA::TimeInterval::negAbsValue(TiRef);
   TiTstIYear2  = OMEGA::TimeInterval::negAbsValue(TiTstIYear);
   TiTstIMonth2 = OMEGA::TimeInterval::negAbsValue(TiTstIMonth);
   TiTstIDay2   = OMEGA::TimeInterval::negAbsValue(TiTstIDay);

   Err1 = TiTstFS.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == -WRef && NTst == -NRef && DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: negAbsValue: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: negAbsValue: FAIL");
   }

   Err1 = TiTstIYear2.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == -IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: negAbsValue years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: negAbsValue years: FAIL");
   }

   Err1 = TiTstIMonth2.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == -IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: negAbsValue months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: negAbsValue months: FAIL");
   }

   Err1 = TiTstIDay2.get(ITst, OMEGA::TimeUnits::Days);
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

   OMEGA::TimeInterval TiTst3;
   OMEGA::TimeInterval TiTstIYear3;
   OMEGA::TimeInterval TiTstIMonth3;
   OMEGA::TimeInterval TiTstIDay3;

   TiTst3       = OMEGA::TimeInterval::absValue(TiTstFS);
   TiTstIYear3  = OMEGA::TimeInterval::absValue(TiTstIYear2);
   TiTstIMonth3 = OMEGA::TimeInterval::absValue(TiTstIMonth2);
   TiTstIDay3   = OMEGA::TimeInterval::absValue(TiTstIDay2);

   Err1 = TiTst3.get(WTst, NTst, DTst);
   if (Err1 == 0 && WTst == WRef && NTst == NRef && DTst == DRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: absValue: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: absValue: FAIL");
   }

   Err1 = TiTstIYear3.get(ITst, OMEGA::TimeUnits::Years);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: absValue years: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: absValue years: FAIL");
   }

   Err1 = TiTstIMonth3.get(ITst, OMEGA::TimeUnits::Months);
   if (Err1 == 0 && ITst == IRef) {
      LOG_INFO("TimeMgrTest/TimeInterval: absValue months: PASS");
   } else {
      ++ErrAll;
      LOG_ERROR("TimeMgrTest/TimeInterval: absValue months: FAIL");
   }

   Err1 = TiTstIDay3.get(ITst, OMEGA::TimeUnits::Days);
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

   return ErrAll;

} // end testTimeInterval

//------------------------------------------------------------------------------
// The test driver.

int main(int argc, char *argv[]) {

   OMEGA::I4 Err{0};
   OMEGA::I4 TotErr{0};

   Err = testTimeFrac();
   TotErr += Err;

   Err = testCalendar();
   TotErr += Err;

   Err = testTimeInterval();
   TotErr += Err;

   if (TotErr == 0) {
      LOG_INFO("TimeMgrTest: Successful completion");
   } else {
      LOG_INFO("TimeMgrTest: Failed");
   }

} // end of main
//===-----------------------------------------------------------------------===/
