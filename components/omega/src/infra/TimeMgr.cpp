//===-- infra/TimeMgr.cpp - time manager ------------------------*- C++ -*-===//
//
// The OMEGA time manager tracks simulation time during a forward-integration
// in a manner to avoid accumulated roundoff in simulations that can extend
// over millions of time steps, is compatible with numerous calendar systems,
// and manages alarms for triggering events at precise times. This time manager
// module consists of six classes:
// 1. The TimeFrac class is the base time representation for a number of
//    seconds as an integer fraction
// 2. The Calendar class holds important information for supported
//    calendar kinds
// 3. The TimeInstant class represents a point in time on a given calendar
// 4. The TimeInterval class represents the amount of time between different
//    points in time
// 5. The Alarm class allows the user to trigger events with single-use or
//    periodic alarms
// 6. The Clock class tracks time during an integration and helps manage
//    attached alarms
//
// Based in part on ESMF time manager, so copyright attached here:
// Earth System Modeling Framework
// Copyright 2002-2024, University Corporation for Atmospheric Research,
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics
// Laboratory, University of Michigan, National Centers for Environmental
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
// NASA Goddard Space Flight Center.
//
//===----------------------------------------------------------------------===//

#include "TimeMgr.h"
#include "Logging.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <memory>

// max/min macros if they don't already exist
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

namespace OMEGA {

//------------------------------------------------------------------------------
// Utility function to extract time units from a string
TimeUnits TimeUnitsFromString(
    const std::string TimeUnitString ///< string describing time units
) {

   // First convert to lower case for easier comparison
   std::string CompStr = TimeUnitString;
   std::transform(CompStr.begin(), CompStr.end(), CompStr.begin(),
                  [](unsigned char C) { return std::tolower(C); });

   // Now convert the string to the associated enum TimeUnit

   if (CompStr == "year" or CompStr == "years" or CompStr == "nyears") {
      return TimeUnits::Years;
   } else if (CompStr == "month" or CompStr == "months" or
              CompStr == "nmonths") {
      return TimeUnits::Months;
   } else if (CompStr == "day" or CompStr == "days" or CompStr == "ndays") {
      return TimeUnits::Days;
   } else if (CompStr == "hour" or CompStr == "hours" or CompStr == "nhours") {
      return TimeUnits::Hours;
   } else if (CompStr == "minute" or CompStr == "minutes" or
              CompStr == "nminutes") {
      return TimeUnits::Minutes;
   } else if (CompStr == "second" or CompStr == "seconds" or
              CompStr == "nseconds") {
      return TimeUnits::Seconds;
   } else {
      return TimeUnits::None;
   }
};

//------------------------------------------------------------------------------
// TimeFrac definitions
//------------------------------------------------------------------------------

// Utility functions that don't operate on object
//------------------------------------------------------------------------------
// TimeFracGCD - Compute greatest common denominator for base time fractions
// Uses Euclid's algorithm to determine the Greatest Common Divisor of
// A and B.

I8 TimeFracGCD(const I8 A,   // [in] first  denominator
               const I8 B) { // [in] second denominator

   // deal with a zero input
   if (A == 0 && B != 0)
      return llabs(B);
   else if (A != 0 && B == 0)
      return llabs(A);
   else if (A == 0 && B == 0)
      return 1;

   I8 Abs_A = llabs(A);
   I8 Abs_B = llabs(B);
   I8 Large = MAX(Abs_A, Abs_B);
   I8 Small = MIN(Abs_A, Abs_B);
   I8 Remainder;

   // initial remainder
   Remainder = Large % Small;

   // divide smaller number by previous remainder until remainder goes to 0
   while (Remainder != 0) {
      Large     = Small;
      Small     = Remainder;
      Remainder = Large % Small;
   }

   // the GCD is the last non-zero remainder or the passed-in
   // smallest number
   return Small;

} // end TimeFracGCD

//------------------------------------------------------------------------------
// TimeFracLCM -  Compute least common multiple for base time fractions
// Uses GCD to determine the Least Common Multiple of A and B.
// LCM = (A * B) / GCD(A,B)

I8 TimeFracLCM(const I8 A,   // first  num for which LCM needed
               const I8 B) { // second num for which LCM needed

   I8 GCD = TimeFracGCD(A, B);
   // Check GCD does not return zero
   if (GCD == 0) {
      LOG_ERROR("TimeMgr: TimeFracLCM TimeFracGCD returned 0.");
      return 0;
   }

   return llabs((A / GCD) * B); // avoid (a * b) directly to prevent
                                //   overflow when a and b are large;
                                //   return absolute value

} // end TimeFracLCM

// TimeFrac accessors
//------------------------------------------------------------------------------
// TimeFrac::set - sets base time components directly

I4 TimeFrac::set(I8 W,   // [in] whole integer seconds
                 I8 N,   // [in] fractional seconds, numerator
                 I8 D) { // [in] fractional seconds, denominator

   // Initialize error code to success
   I4 Err{0};

   // Whole and Numer must be either both positive or both negative (product
   // non-neg), Denom must be always positive and >= 1. If conditions are met,
   // set TimeFrac to input values, else output an error message

   if ((W * N >= 0) && D >= 1) {
      Whole = W;          // Sets whole number to input value
      Numer = N;          // Sets numerator    to input value
      Denom = D;          // Sets denominator  to input value
      Err   = simplify(); // Ensure fraction in simplest form
   } else {
      LOG_ERROR("TimeMgr: Invalid input for TimeFrac::set, product of W and N "
                "must be non-negative and D must be >= 1.");

      Err = 1;
   }

   return Err;

} // end TimeFrac::set

//------------------------------------------------------------------------------
// TimeFrac::setHMS - sets a base time by converting from hour, min, sec

I4 TimeFrac::setHMS(I4 Hours,     // [in] integer hours
                    I4 Minutes,   // [in] integer minutes
                    I4 Seconds) { // [in] integer seconds

   // Initialize error code as success
   I4 Err{0};

   // set fractional component to zero
   Numer = 0;
   Denom = 1;

   // convert all inputs to seconds and set final whole seconds to sum
   Whole = (I8)Seconds + (I8)Hours * SECONDS_PER_HOUR +
           (I8)Minutes * SECONDS_PER_MINUTE;

   Err = simplify(); // Ensure fraction in simplest form
   if (Err != 0)
      LOG_ERROR("TimeMgr: Error simplifying final fraction");

   // return
   return Err;

} // end TimeFrac::setHMS

//------------------------------------------------------------------------------
// TimeFrac::setSeconds - sets a base time by converting real seconds

I4 TimeFrac::setSeconds(R8 Seconds) { // [in] floating point seconds

   // Initialize error code as success
   I4 Err{0};

   // Ensure real input is in reasonable fractional range
   // due to limitation of I8 representation.
   R8 Rabs = fabs(Seconds);
   if ((Rabs > 0.0 && Rabs < 1e-17) || Rabs > 1e18) {
      Err = 1;
      LOG_ERROR("TimeMgr: Input out of range, abs value larger than 1e18 or "
                "smaller than 1e-17");
      return Err;
   }

   // Must convert real seconds to native integer fraction.
   // This is performed using a method of continued fractions (CF) and a
   // standard variation of Euclid's GCD algorithm.

   // save sign to restore later
   I4 Sign = (Seconds < 0) ? -1 : 1;

   // Initialize algorithm
   Whole     = 0;
   Numer     = 0;
   Denom     = 1;
   R8 Target = Rabs;
   if (Target == 0.0)
      return Err; // if input is 0.0, we're done

   // Strip off any whole number part (w) first to avoid 64-bit overflow in
   // the algorithm if given a large value.
   if (Target >= 1.0) {
      I8 W = (I8)Rabs;
      Target -= (R8)W;
      Whole = Sign * W;
      if (Target < 1e-17)
         return Err; // if input is an integer, we're done
   }

   // Target precision is relative to given input (not target, due to possible
   // loss of significant digits from the subraction of whole).  Precision is
   // 1/10 of 15th (DBL_DIG) decimal digit (least significant), which is
   // about half a decimal digit more than in double precison, or within half
   // of DBL_EPSILON (2.2e-16) of input magnitude (< 0.5*fabs(2.2e-16 * input)).
   // Anything more stringent can cause instability in the algorithm such as
   // infinite looping due to 64-bit overflow or computation of unreasonably
   // large n and d values.  Test cases: rin = (8+21/23)-8, rin = 9.1 - 9.

   R8 P = pow(10.0, -(DBL_DIG - (I4)log10(Rabs)));

   R8 R         = Target;
   I8 Nprevprev = 0;
   I8 Nprev     = 1;
   I8 Dprevprev = 1;
   I8 Dprev     = 0;

   // Iterate until double precision representation to the number of decimal
   // significant digits (at least DBL_DIG) is achieved.
   // Worst case is 37 iterations for the golden ratio
   // phi = (1 + sqrt(5))/2; the "most irrational number".

   I8 A;
   I8 N; // numerator
   I8 D; // denominator
   R8 F = 0.0;
   I4 I = 0;
   do {
      // Compute next convergent N/D
      A = (I8)R;
      N = A * Nprev + Nprevprev;
      D = A * Dprev + Dprevprev;
      if (fabs((R8)N / (R8)D - Target) < P)
         break; // done when double precision reached

      F = R - (R8)A;
      // Prevent divide-by-zero (exceed do-able precision) Will this condition
      // ever be met, while primary target termination condition above is not?
      // I.e., perhaps this line is not necessary, and the line above and the
      // line below can be combined into one, eliminating the variable f.
      if (F < 1e-17)
         break;

      // Prepare for next iteration
      R         = 1.0 / F;
      Nprevprev = Nprev;
      Nprev     = N;
      Dprevprev = Dprev;
      Dprev     = D;
   } while (true);

   // done iterating, set final value (whole number already set above)
   Numer = N * Sign;
   Denom = D;

   Err = simplify(); // Ensure fraction in simplest form
   if (Err != 0)
      LOG_ERROR("TimeMgr: Error simplifying final fraction");
   return Err; // Return error code

} // end TimeFrac::setSeconds

//------------------------------------------------------------------------------
// TimeFrac::setMinutes - sets a base time by converting real minutes

I4 TimeFrac::setMinutes(R8 minutes) { // [in] floating point minutes

   // Initialize error code as success
   I4 Err{0};

   // Convert to real seconds and use the seconds routine
   R8 seconds = minutes * SECONDS_PER_MINUTE;
   Err        = setSeconds(seconds);

   // Check for errors and return
   if (Err != 0)
      LOG_ERROR("TimeMgr: Error calling TimeFrac::setSeconds");
   return Err;

} // end TimeFrac::setMinutes

//------------------------------------------------------------------------------
// TimeFrac::setHours - sets a base time by converting real hours

I4 TimeFrac::setHours(R8 hours) { // [in] floating point hours

   // Initialize error code as success
   I4 Err{0};

   // Convert to real seconds and use the seconds routine
   R8 seconds = hours * SECONDS_PER_HOUR;
   Err        = setSeconds(seconds);

   // Check for errors and return
   if (Err != 0)
      LOG_ERROR("TimeMgr: Error calling TimeFrac::setSeconds");
   return Err;

} // end TimeFrac::setHours

//------------------------------------------------------------------------------
// TimeFrac::setWhole - Set whole seconds separately

I4 TimeFrac::setWhole(const I8 W) { // [in] Whole number of seconds

   Whole = W; // Sets whole number to input value
   return 0;

} // end TimeFrac::setWhole

//------------------------------------------------------------------------------
// TimeFrac::setNumer - Set numerator of fractional seconds separately

I4 TimeFrac::setNumer(const I8 N) { // [in] Numerator of fractional seconds

   Numer = N; // Sets numerator to input value
   return 0;

} // end TimeFrac::setNumer

//------------------------------------------------------------------------------
// TimeFrac::setDenom - Set denominator of fractional seconds separately

I4 TimeFrac::setDenom(const I8 D) { // [in] Denominator of fraction

   Denom = D; // Sets denominator to input value
   return 0;

} // end TimeFrac::setDenom

//------------------------------------------------------------------------------
// TimeFrac::get - Retrieve base time in native fractional seconds

I4 TimeFrac::get(I8 &W,         // [out] whole seconds
                 I8 &N,         // [out] fractional second numerator
                 I8 &D) const { // [out] fractional second denominator

   // retrieve components
   W = Whole;
   N = Numer;
   D = Denom;

   return 0;

} // end TimeFrac::get

//------------------------------------------------------------------------------
// TimeFrac::getHMS - Gets base time in integer hours, minutes, seconds
// Note that this version truncates the fractional second (rounds toward zero).

I4 TimeFrac::getHMS(I4 &Hours,           // [out] integer hours
                    I4 &Minutes,         // [out] integer minutes
                    I4 &Seconds) const { // [out] integer seconds

   // Initialize error code to success
   I4 Err{0};

   // make local copy for manipulation
   TimeFrac RemainingTime = *this;
   RemainingTime.simplify(); // ensure maximum whole seconds

   // Simce result must be integer, we drop the fractional seconds
   I8 RemainingSeconds = RemainingTime.Whole;

   // extract hours and check whether in range
   I8 TmpHours = RemainingSeconds / SECONDS_PER_HOUR;

   if (TmpHours < INT_MIN || TmpHours > INT_MAX) {
      Err = 1;
      LOG_ERROR("TimeMgr: value exceeds machine limits for I4");
      return Err;
   }

   // hours from remaining time
   RemainingSeconds %= SECONDS_PER_HOUR;

   // now extract minutes and remove from remaining time
   I8 TmpMinutes = RemainingSeconds / SECONDS_PER_MINUTE;
   RemainingSeconds %= SECONDS_PER_MINUTE;

   // convert results for return
   Hours   = (I4)TmpHours;
   Minutes = (I4)TmpMinutes;
   Seconds = (I4)RemainingSeconds;

   // return error code
   return Err;

} // end TimeFrac::getHMS

//------------------------------------------------------------------------------
// TimeFrac::getSeconds - Retrieves a base time and converts to a real number of
// seconds

R8 TimeFrac::getSeconds(void) const { // \result Time in real seconds

   // check for divide-by-zero
   I4 Err{0};
   if (Denom == 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: encountered 0 denominator.");
      return 0.0;
   }

   // convert integer fractional seconds to real result
   return (R8)Whole + (R8)Numer / (R8)Denom;

} // end TimeFrac::getSeconds

//------------------------------------------------------------------------------
// TimeFrac::getHours - Retrieves a base time and converts to a real number of
// hours

R8 TimeFrac::getHours(void) const { // \result Time in real hours

   // check for divide-by-zero
   I4 Err{0};
   if (Denom == 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: encountered 0 denominator.");
      return 0.0;
   }

   // make local copy for manipulation
   TimeFrac TmpTime = *this;

   // avoid error introduced by floating point divide or multiply by
   // performing a fractional division in integer arithmetic first,
   // then converting to floating point.
   TmpTime /= SECONDS_PER_HOUR;

   // convert integer fractional seconds to real result
   return (R8)TmpTime.Whole + (R8)TmpTime.Numer / (R8)TmpTime.Denom;

} // end TimeFrac::getHours

//------------------------------------------------------------------------------
// TimeFrac::getMinutes - Retrieves a base time and converts to a real number of
// minutes

R8 TimeFrac::getMinutes(void) const { // \result Time in real minutes

   // check for divide-by-zero
   int Err{0};
   if (Denom == 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: encountered 0 denominator.");
      return 0.0;
   }

   // make local copy for manipulation
   TimeFrac TmpTime = *this;

   // avoid error introduced by floating point divide or multiply by
   // performing a fractional division in integer arithmetic first,
   // then converting to floating point.
   TmpTime /= SECONDS_PER_MINUTE;

   // convert integer fractional seconds to real result
   return (R8)TmpTime.Whole + (R8)TmpTime.Numer / (R8)TmpTime.Denom;

} // end TimeFrac::getMinutes

//------------------------------------------------------------------------------
// TimeFrac::getWhole - Retrieve the whole seconds component of base time

I8 TimeFrac::getWhole(void) const { return Whole; }

//------------------------------------------------------------------------------
// TimeFrac::getNumer - Retrieve the numerator of the fractional component of
// base time

I8 TimeFrac::getNumer(void) const { return Numer; }

//------------------------------------------------------------------------------
// TimeFrac::getDenom - Retrieve the denominator of the fractional component of
// base time

I8 TimeFrac::getDenom(void) const { return Denom; }

// TimeFrac constructors/destructors
//------------------------------------------------------------------------------
// TimeFrac - default constructor defines a zero time fraction

TimeFrac::TimeFrac(void)
    : Whole(0), Numer(0), Denom(1) { // Denom = 1 to prevent divide-by-zero
} // end TimeFrac default constructor

//------------------------------------------------------------------------------
// TimeFrac - Copy constructor

TimeFrac::TimeFrac(const TimeFrac &TimeToCopy) // [in] TimeFrac to copy
    : Whole(TimeToCopy.Whole), Numer(TimeToCopy.Numer),
      Denom(TimeToCopy.Denom) {} // end TimeFrac copy constructor

//------------------------------------------------------------------------------
// TimeFrac - Constructs TimeFrac by fractional second components

TimeFrac::TimeFrac(I8 W,   // [in] integer whole seconds
                   I8 N,   // [in] fractional seconds, numerator
                   I8 D) { // [in] fractional seconds, denominator

   // Calls the set function to set TimeFrac to input values
   I4 Err = this->set(W, N, D);

   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeFrac constructor returned error from component "
                "set routine");

} // end TimeFrac constructor by components

//------------------------------------------------------------------------------
// TimeFrac - Constructs TimeFrac from a real number of seconds

TimeFrac::TimeFrac(const R8 Seconds) { // [in] Time in real seconds

   // Calls the setSecond function to convert seconds to an integer
   // fraction for TimeFrac
   I4 Err = this->setSeconds(Seconds);

   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeFrac constructor (real seconds) Returned error "
                "from component set routine");

} // end TimeFrac constructor from real seconds

//----------------------------------------------------------------------
// ~TimeFrac - destructor for base time

TimeFrac::~TimeFrac(void) { // No operation - nothing to clean up
} // end ~TimeFrac

// TimeFrac operators
//------------------------------------------------------------------------------
// TimeFrac::operator(==) - Equivalence comparison operator for TimeFrac

bool TimeFrac::operator==(
    const TimeFrac &Time) const { // [in] TimeFrac to compare

   // make local copies; don't change the originals.
   TimeFrac Ttmp1 = *this;
   TimeFrac Ttmp2 = Time;

   // ensure proper base time
   if (Ttmp1.simplify() != 0 || Ttmp2.simplify() != 0) {
      LOG_ERROR("TimeMgr: TimeFrac(==) simplfy error");
      return false;
   }

   // put both fractions on the same denominator, then compare
   I8 LCM = TimeFracLCM(Ttmp1.Denom, Ttmp2.Denom);
   return Ttmp1.Whole == Ttmp2.Whole && Ttmp1.Numer * (LCM / Ttmp1.Denom) ==
                                            Ttmp2.Numer * (LCM / Ttmp2.Denom);

} // end TimeFrac::operator(==)

//------------------------------------------------------------------------------
// TimeFrac::operator(!=)
// Non-equivalence comparison operator for TimeFrac

bool TimeFrac::operator!=(
    const TimeFrac &Time) const { // [in] TimeFrac to compare

   // make local copies; don't change the originals.
   TimeFrac Ttmp1 = *this;
   TimeFrac Ttmp2 = Time;

   // ensure proper base time
   if (Ttmp1.simplify() != 0 || Ttmp2.simplify() != 0) {
      LOG_ERROR("TimeMgr: TimeFrac(!=) simplfy error");
      return false;
   }

   // put both fractions on the same denominator, then compare
   I8 LCM = TimeFracLCM(Ttmp1.Denom, Ttmp2.Denom);
   return Ttmp1.Whole != Ttmp2.Whole || Ttmp1.Numer * (LCM / Ttmp1.Denom) !=
                                            Ttmp2.Numer * (LCM / Ttmp2.Denom);

} // end TimeFrac::operator(!=)

//------------------------------------------------------------------------------
// TimeFrac::operator(<)
// Less than comparison operator for TimeFrac

bool TimeFrac::operator<(
    const TimeFrac &Time) const { // [in] TimeFrac to compare

   // make local copies; don't change the originals.
   TimeFrac Ttmp1 = *this;
   TimeFrac Ttmp2 = Time;

   // ensure proper fractions; check for divide-by-zero
   if (Ttmp1.simplify() != 0 || Ttmp2.simplify() != 0) {
      LOG_ERROR("TimeMgr: TimeFrac(<) simplfy error");
      return false;
   }

   // ignore fractional part if whole parts are different
   if (Ttmp1.Whole != Ttmp2.Whole) {
      return Ttmp1.Whole < Ttmp2.Whole;
   } else { // must look at fractional part
      // put both fractions on the same denominator, then compare
      I8 LCM = TimeFracLCM(Ttmp1.Denom, Ttmp2.Denom);
      return Ttmp1.Numer * (LCM / Ttmp1.Denom) <
             Ttmp2.Numer * (LCM / Ttmp2.Denom);
   }

} // end TimeFrac::operator(<)

//------------------------------------------------------------------------------
// TimeFrac::operator(>)
// Greater than comparison operator for TimeFrac

bool TimeFrac::operator>(
    const TimeFrac &Time) const { // [in] TimeFrac to compare

   // make local copies so we do not change the originals.
   TimeFrac Ttmp1 = *this;
   TimeFrac Ttmp2 = Time;

   // ensure proper fractions; check for divide-by-zero
   if (Ttmp1.simplify() != 0 || Ttmp2.simplify() != 0) {
      LOG_ERROR("TimeMgr: TimeFrac(<) simplfy error");
      return false;
   }

   // ignore fractional part if whole parts are different
   if (Ttmp1.Whole != Ttmp2.Whole) {
      return Ttmp1.Whole > Ttmp2.Whole;
   }
   // must look at fractional part
   else {
      // put both fractions on the same denominator, then compare
      I8 LCM = TimeFracLCM(Ttmp1.Denom, Ttmp2.Denom);
      return Ttmp1.Numer * (LCM / Ttmp1.Denom) >
             Ttmp2.Numer * (LCM / Ttmp2.Denom);
   }

} // end TimeFrac::operator(>)

//------------------------------------------------------------------------------
// TimeFrac::operator(<=)
// Less than or equal comparison operator for TimeFrac

bool TimeFrac::operator<=(
    const TimeFrac &Time) const { // [in] TimeFrac to compare

   // reuse < and == operators defined above
   return *this < Time || *this == Time;

} // end TimeFrac::operator(<=)

//------------------------------------------------------------------------------
// TimeFrac::operator(>=)
// Greater than or equal comparison operator for TimeFrac

bool TimeFrac::operator>=(
    const TimeFrac &Time) const { // [in] TimeFrac to compare

   // reuse > and == operators defined above
   return *this > Time || *this == Time;

} // end TimeFrac::operator(>=)

//------------------------------------------------------------------------------
// TimeFrac::operator(+)
// Addition operator for TimeFrac

TimeFrac
TimeFrac::operator+(const TimeFrac &Time) const { // [in] TimeFrac to add

   // check for divide-by-zero
   if (Denom == 0 || Time.Denom == 0) {
      LOG_ERROR("TimeMgr: TimeFrac(+) encountered 0 denominator");
      return TimeFrac(0, 0, 1);
   }

   TimeFrac Sum;

   // fractional part addition
   Sum.Denom = TimeFracLCM(Denom, Time.Denom);
   Sum.Numer =
       Numer * (Sum.Denom / Denom) + Time.Numer * (Sum.Denom / Time.Denom);

   // whole part addition
   Sum.Whole = Whole + Time.Whole;

   // ensure simplified form
   Sum.simplify();

   return Sum;

} // end TimeFrac::operator(+)

//------------------------------------------------------------------------------
// TimeFrac::operator(-)
// Subtraction operator for TimeFrac

TimeFrac
TimeFrac::operator-(const TimeFrac &Time) const { // [in] TimeFrac to subtract

   // check for divide-by-zero
   if (Denom == 0 || Time.Denom == 0) {
      LOG_ERROR("TimeMgr: TimeFrac(-) encountered 0 denominator");
      return TimeFrac(0, 0, 1);
   }

   TimeFrac Diff;

   // fractional part subtraction
   // must convert to same denominator
   Diff.Denom = TimeFracLCM(Denom, Time.Denom);
   Diff.Numer =
       Numer * (Diff.Denom / Denom) - Time.Numer * (Diff.Denom / Time.Denom);

   // whole part subtraction
   Diff.Whole = Whole - Time.Whole;

   // ensure simplified form
   Diff.simplify();

   return Diff;

} // end TimeFrac::operator(-)

//------------------------------------------------------------------------------
// TimeFrac::operator(+=)
// Increment operator for TimeFrac

TimeFrac &
TimeFrac::operator+=(const TimeFrac &Time) { // [in] TimeFrac increment

   // reuse (+) operator defined above
   *this = *this + Time;
   return *this;

} // end TimeFrac::operator(+=)

//------------------------------------------------------------------------------
// TimeFrac::operator(-=)
// Decrement operator for TimeFrac

TimeFrac &
TimeFrac::operator-=(const TimeFrac &Time) { // [in] TimeFrac decrement

   // reuse (-) operator defined above
   *this = *this - Time;
   return *this;

} // end TimeFrac::operator(-=)

//------------------------------------------------------------------------------
// TimeFrac::operator(*)
// Multiplication by integer scalar

TimeFrac
TimeFrac::operator*(const I4 Multiplier) const { // [in] integer multiplier

   TimeFrac Product;

   // fractional part multiplication.
   Product.Numer = Numer * Multiplier;
   Product.Denom = Denom;

   // whole part multiplication
   Product.Whole = Whole * Multiplier;

   // ensure simplified form
   Product.simplify();

   return Product;

} // end TimeFrac::operator(*)

//------------------------------------------------------------------------------
// TimeFrac::operator(*=)
// Multiplication in place by integer scalar

TimeFrac &TimeFrac::operator*=(const I4 Multiplier) { // [in] integer multiplier

   // just reuse (*) operator defined above
   *this = *this * Multiplier;
   return *this;

} // end TimeFrac::operator(*=)

//------------------------------------------------------------------------------
// TimeFrac::operator(*)
// Multiplication by real scalar

TimeFrac
TimeFrac::operator*(const R8 Multiplier) const { // [in] real multiplier

   TimeFrac Product; // initialized to zero

   // convert real scalar to a base time using conversion in setSeconds
   TimeFrac Mult;
   Mult.setSeconds(Multiplier);

   // now multiply the two fractions
   Product.Denom = Denom * Mult.Denom;
   Product.Numer =
       (Whole * Denom + Numer) * (Mult.Whole * Mult.Denom + Mult.Numer);

   // ensure simplified form
   Product.simplify();

   return Product;

} // end TimeFrac::operator(*)

//------------------------------------------------------------------------------
// TimeFrac::operator(*=)
// Multiplication in place by real scalar

TimeFrac &TimeFrac::operator*=(const R8 Multiplier) { // [in] real multiplier

   // reuse (*) operator defined above
   *this = *this * Multiplier;
   return *this;

} // end TimeFrac::operator(*=)

//------------------------------------------------------------------------------
// TimeFrac::operator(/)
// Divide TimeFrac by integer scalar

TimeFrac TimeFrac::operator/(I4 Divisor) const { // [in] integer divisor

   // check for divide-by-zero
   if (Divisor == 0) {
      LOG_ERROR("TimeMgr: TimeFrac(/) attempt to divide by 0");
      return TimeFrac(0, 0, 1);
   }

   TimeFrac Quotient;
   I8 Remainder;
   I8 Denominator;

   // Fractional part division. To avoid overflows/underflows with large
   // denominators, we do not just blindly multiply denominator.
   // Instead, divide numerator and add back any remainder.
   Quotient.Numer = Numer / (I8)Divisor;
   Quotient.Denom = Denom;
   Remainder      = Numer % (I8)Divisor;

   // check remainder and add back
   if (Remainder != 0) {
      // upper bounds check of (d * divisor)
      Denominator = Denom * (I8)Divisor; // must do it!

      // add back
      Quotient += TimeFrac(0, Remainder, Denominator);
   }

   // whole part division; add back any remainder
   Quotient.Whole = Whole / (I8)Divisor;
   Remainder      = Whole % (I8)Divisor;
   if (Remainder != 0) {
      Quotient += TimeFrac(0, Remainder, (I8)Divisor);
   }

   // ensure simplified form
   Quotient.simplify();

   return Quotient;

} // end TimeFrac::operator(/)

//------------------------------------------------------------------------------
// TimeFrac::operator(/=)
// Divide TimeFrac in place by integer scalar

TimeFrac &TimeFrac::operator/=(I4 divisor) { // [in] integer divisor

   // just reuse (/) operator defined above
   *this = *this / divisor;
   return *this;

} // end TimeFrac::operator(/=)

//------------------------------------------------------------------------------
// TimeFrac::operator(/)
// Divide two TimeFracs and return a double precision result

R8 TimeFrac::operator/(
    const TimeFrac &Time) const { // [in] TimeFrac to divide by

   // check for divide-by-zero
   if (Denom == 0 || Time.Denom == 0 ||
       Time.Whole * Time.Denom + Time.Numer == 0) {
      LOG_ERROR("TimeMgr: TimeFrac(/) attempt to divide by 0");
      return 0.0;
   }

   R8 Quotient = (Whole + (R8)Numer / (R8)Denom) /
                 (Time.Whole + (R8)Time.Numer / (R8)Time.Denom);

   return Quotient;

} // end TimeFrac::operator(/)

//------------------------------------------------------------------------------
// TimeFrac::operator(%)
// Computes the modulus of two fractions

TimeFrac
TimeFrac::operator%(const TimeFrac &Time) const { // [in] TimeFrac to modulo by

   TimeFrac Remainder;

   // if no fractional part, just modulus the whole parts
   if (Numer == 0 && Time.Numer == 0) {
      // check for divide-by-zero
      if (Time.Whole != 0) {
         Remainder.Numer = Whole % Time.Whole;
      } else {
         LOG_ERROR("TimeMgr: TimeFrac modulus attempt to divide by 0");
         return TimeFrac(0, 0, 1);
      }
   }
   // otherwise, perform Time modulus
   else {
      // check for divide-by-zero
      if (Denom == 0 || Time.Denom == 0) {
         LOG_ERROR("TimeMgr: TimeFrac modulus divide by 0");
         return TimeFrac(0, 0, 1);
      }

      I8 LCM = TimeFracLCM(Denom, Time.Denom);

      // convert *this Time and given Time into improper form with a
      // common denominator and then compute remainder
      Remainder.Denom =
          (Time.Whole * Time.Denom + Time.Numer) * (LCM / Time.Denom);
      if (Remainder.Denom != 0) {
         Remainder.Numer =
             ((Whole * Denom + Numer) * (LCM / Denom)) % Remainder.Denom;
      } else {
         LOG_ERROR("TimeMgr: TimeFrac modulus divide by 0");
         return TimeFrac(0, 0, 1);
      }
   }

   // ensure simplified form
   Remainder.simplify();

   return Remainder;

} // end TimeFrac::operator(%)

//------------------------------------------------------------------------------
// TimeFrac::(%=)
// Compute modulus of two fractions in place

TimeFrac &
TimeFrac::operator%=(const TimeFrac &Time) { // [in] TimeFrac to modulo by

   // just reuse (%) operator defined above!
   *this = *this % Time;
   return *this;

} // end TimeFrac::operator(%=)

//------------------------------------------------------------------------------
// TimeFrac::operator(=)
// Assignment operator for TimeFrac

TimeFrac &TimeFrac::operator=(const TimeFrac &Time) { // [in] TimeFrac to assign

   if (&Time != this) {
      this->Whole = Time.Whole;
      this->Numer = Time.Numer;
      this->Denom = Time.Denom;
   }

   return *this;

} // end TimeFrac::operator=

// TimeFrac utility methods
//------------------------------------------------------------------------------
// TimeFrac::convert
// Converts a base time fraction to new denominator

I4 TimeFrac::convert(I8 Denominator) {

   I4 Err{0};

   // Check for invalid denominator
   if (Denominator == 0 || Denom == 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: TimeFrac::convert encountered 0 denominator");
      return Err;
   }

   I8 Conversion = Whole * Denominator + (Numer * Denominator) / Denom;

   // set new values
   Whole = 0;
   Numer = Conversion;
   Denom = Denominator;

   // don't simplify; leave on specified denominator.

   return Err; // return error code

} // end TimeFrac::convert

//------------------------------------------------------------------------------
// TimeFrac::simplify
// Ensure proper fraction and reduce to lowest denominator
// If fraction >= 1, add to whole part, and adjust fraction to
// remainder. Then reduce to lowest denominator.

I4 TimeFrac::simplify(void) {

   // Initialize error code
   I4 Err{0};

   // check for divide-by-zero
   if (Denom == 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: TimeFrac::simplify encountered 0 denominator");
      return Err;
   }

   // normalize to proper fraction - if fraction > 1 add to the
   // whole number component to bring fraction under 1
   I8 W;
   if (llabs((W = Numer / Denom)) >= 1) {
      Whole += W;
      Numer %= Denom;
   }

   // ensure whole and fraction parts are same sign

   // if whole is positive and fraction is negative
   if (Whole > 0 && ((Numer < 0 && Denom > 0) || (Denom < 0 && Numer > 0))) {
      Whole--;        // subtract one from whole number
      Numer += Denom; //  and add it to the fraction part

      // else if whole is negative and fraction is positive
   } else if ((Whole < 0 && (Numer > 0 && Denom > 0)) ||
              (Denom < 0 && Numer < 0)) {
      Whole++;        // add one to whole number
      Numer -= Denom; // and subtract it from the fraction part
   }

   // normalize fraction sign
   if (Denom < 0) {
      Denom *= -1; // change signs
      Numer *= -1;
   }

   // reduce to lowest denominator

   I8 GCD = TimeFracGCD(Numer, Denom);
   // Check GCD does not return zero
   if (GCD == 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: TimeFrac::simplify TimeFracGCD returned 0");
      return Err;
   }

   Numer /= GCD;
   Denom /= GCD;

   return Err; // return error code

} // end TimeFrac::simplify

//------------------------------------------------------------------------------
// Calendar definitions
//------------------------------------------------------------------------------

// initialize static calendar instance
std::unique_ptr<Calendar> Calendar::OmegaCal = nullptr;

// Calendar accessors
//------------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Calendar::get - retrieves pointer to model calendar
Calendar *Calendar::get() {
   if (isDefined()) {
      return Calendar::OmegaCal.get();
   } else {
      LOG_CRITICAL("Attempt to retrieve undefined calendar");
   }
} // end retrieve calendar pointer

//------------------------------------------------------------------------------
// Calendar::get - retrievals for each calendar property

CalendarKind Calendar::getKind() {

   CalendarKind OutKind = CalendarUnknown;
   // check if calendar defined
   if (isDefined()) {

      // get calendar type and check for valid value
      OutKind = Calendar::OmegaCal->CalKind;
      if (OutKind == CalendarUnknown)
         LOG_CRITICAL("Cannot retrieve calendar kind"
                      " - Calendar undefined or not initialized");

   } else {
      LOG_CRITICAL("Cannot retrieve calendar kind - Calendar not initialized");
   }

   return OutKind;
}

std::vector<I4> Calendar::getDaysPerMonth() {

   // check if calendar defined
   if (isDefined) {
      // also check for valid vector
      if (Calendar::OmegaCal->DaysPerMonth.size() > 0) {
         return Calendar::OmegaCal->DaysPerMonth;
      } else {
         LOG_CRITICAL("Invalid DaysPerMonth vector");
      }
   } else {
      LOG_CRITICAL("Cannot retrieve calendar props - Calendar not initialized");
   }
}

I4 Calendar::getMonthsPerYear() {

   // check if calendar defined
   if (isDefined()) {
      return Calendar::OmegaCal->MonthsPerYear;
   } else {
      LOG_CRITICAL("Cannot retrieve MonthsPerYear - Calendar not initialized");
      return 0;
   }
}

I4 Calendar::getSecondsPerDay() {

   // check if calendar defined
   if (isDefined()) {
      return Calendar::OmegaCal->SecondsPerDay;
   } else {
      LOG_CRITICAL("Cannot retrieve SecondsPerDay - Calendar not initialized");
      return 0;
   }
}

I4 Calendar::getSecondsPerYear() {

   // check if calendar defined
   if (isDefined()) {
      return Calendar::OmegaCal->SecondsPerYear;
   } else {
      LOG_CRITICAL("Cannot retrieve SecondsPerYear - Calendar not initialized");
      return 0;
   }
}

I4 Calendar::getDaysPerYear() {
   // check if calendar defined
   if (isDefined()) {
      return Calendar::OmegaCal->DaysPerYear;
   } else {
      LOG_CRITICAL("Cannot retrieve DaysPerYear - Calendar not initialized");
      return 0;
   }
}

// Calendar constructors/destructors
//------------------------------------------------------------------------------
// Calendar::init - creates the calendar from a supported calendar option
void Calendar::init(
    std::string CalendarKindStr // [in] string name for calendar type
) {

   // If calendar has already been set, throw a critical error
   if (isDefined()) {
      LOG_CRITICAL("Omega calendar already set, can not be changed");
      return;
   }

   // Convert string name for type of calendar to a supported kind
   CalendarKind CalKindLoc = CalendarUnknown;
   for (int I = 0; I < NUM_SUPPORTED_CALENDARS; ++I) {
      if (CalendarKindStr == CalendarKindName[I]) {
         CalKindLoc = (CalendarKind)(I);
         break;
      }
   }

   // Call constructor and assign the calendar if it is valid
   if (CalKindLoc == CalendarUnknown) {
      LOG_CRITICAL("Attempt to create calendar of unknown name: {}",
                   CalendarKindStr);
      Calendar::OmegaCal = nullptr;
   } else {
      Calendar::OmegaCal = std::unique_ptr<Calendar>(new Calendar(CalKindLoc));
   }

   LOG_INFO("Omega using calendar type: {}", CalendarKindStr);
}

//------------------------------------------------------------------------------
// Calendar::init - creates the calendar from a custom calendar

void Calendar::init(std::vector<I4> &InDaysPerMonth, // [in] vector days/month
                    int InSecondsPerDay,             // [in] seconds per day
                    int InSecondsPerYear,            // [in] seconds per year
                    int InDaysPerYear                // [in] days per year
) {

   // If calendar has already been set, throw a critical error
   if (isDefined()) {
      LOG_CRITICAL("Omega calendar already set, can not be changed");
      return;
   }

   // Call constructor and assign the single instance
   // instance = std::unique_ptr<Singleton>(new Singleton());
   Calendar::OmegaCal = std::unique_ptr<Calendar>(new Calendar(
       InDaysPerMonth, InSecondsPerDay, InSecondsPerYear, InDaysPerYear));

   LOG_INFO("Omega using a custom calendar");

} // end custom calendar create

//------------------------------------------------------------------------------
// Calendar::Calendar - constructor creates standard calendar by kind
// This routine creates a standard calendar from a selection of
// standard supported calendars (e.g. Gregorian, Julian, no-leap).
// It is the preferred method for creating a calendar.

Calendar::Calendar(CalendarKind InCalKind) { // [in] calendar type

   // set calendar kind
   CalKind     = InCalKind;
   CalKindName = CalendarKindName[CalKind];

   // months per year fixed
   MonthsPerYear = MONTHS_PER_YEAR;

   // set remaining variables based on type (kind) of calendar
   DaysPerMonth.resize(MonthsPerYear);
   switch (CalKind) {
   // All these cases have similar months, days per month
   case CalendarGregorian:
   case CalendarJulian:
   case CalendarNoLeap:
      // note that February by default is still 28 days
      // leap years are treated as exceptions later
      DaysPerMonth[0]  = 31;
      DaysPerMonth[1]  = 28;
      DaysPerMonth[2]  = 31;
      DaysPerMonth[3]  = 30;
      DaysPerMonth[4]  = 31;
      DaysPerMonth[5]  = 30;
      DaysPerMonth[6]  = 31;
      DaysPerMonth[7]  = 31;
      DaysPerMonth[8]  = 30;
      DaysPerMonth[9]  = 31;
      DaysPerMonth[10] = 30;
      DaysPerMonth[11] = 31;
      SecondsPerDay    = SECONDS_PER_DAY;
      DaysPerYear      = 365;
      SecondsPerYear   = SECONDS_PER_DAY * DaysPerYear;
      break;

   case CalendarJulianDay:
   case CalendarModJulianDay:
      // Julian day calendars have no month or year, just
      // days since reference time
      for (I4 I = 0; I < MonthsPerYear; I++)
         DaysPerMonth[I] = 0;
      SecondsPerDay  = SECONDS_PER_DAY;
      SecondsPerYear = 0;
      DaysPerYear    = 0;
      break;

   case Calendar360Day:
      // 12 months of 30 days each
      for (I4 I = 0; I < MonthsPerYear; I++)
         DaysPerMonth[I] = 30;
      SecondsPerDay  = SECONDS_PER_DAY;
      DaysPerYear    = 360;
      SecondsPerYear = SECONDS_PER_DAY * DaysPerYear;
      break;

   case CalendarNoCalendar:
      // no calendar needed, only uses seconds since time zero
      for (I4 I = 0; I < MonthsPerYear; I++)
         DaysPerMonth[I] = 0;
      SecondsPerDay  = 0;
      SecondsPerYear = 0;
      DaysPerYear    = 0;
      break;

   case CalendarCustom:
      // custom requires more info, so separate interface should be used
      LOG_ERROR(
          "TimeMgr: Must use custom constructor for custom calendar type");
      break;

   default:
      // unknown calendar kind, invalidate and return with error
      CalKind     = CalendarUnknown;
      CalKindName = CalendarKindName[CalKind];
      LOG_ERROR("Attempt to create a calendar of unknown kind");
      break;

   } // end switch on CalKind

} // end Calendar::Calendar standard constructor

//------------------------------------------------------------------------------
// Calendar::Calendar - constructs a custom calendar
// In some cases, the standard calendars are inappropriate, so this
// custom constructor permits users to define their own calendar.

Calendar::Calendar(std::vector<I4> &InDaysPerMonth, // [in] array of days/month
                   int InSecondsPerDay,             // [in] seconds per day
                   int InSecondsPerYear,            // [in] seconds per year
                   int InDaysPerYear) {             // [in] days per year

   // define calendar as custom
   CalKind     = CalendarCustom;
   CalKindName = CalendarKindName[CalKind];

   // months per year is fixed
   MonthsPerYear = MONTHS_PER_YEAR;

   // set remaining calendar vars from input values
   DaysPerMonth.resize(MonthsPerYear);
   for (I4 I = 0; I < MonthsPerYear; I++)
      DaysPerMonth[I] = InDaysPerMonth[I];
   SecondsPerDay  = InSecondsPerDay;
   SecondsPerYear = InSecondsPerYear;
   DaysPerYear    = InDaysPerYear;

} // end Calendar::Calendar custom constructor

//------------------------------------------------------------------------------
// Calendar::~Calendar - destructor for calendar

Calendar::~Calendar(void) {}

//------------------------------------------------------------------------------
// Calendar::reset - Destroys the single instance
// This should only be used during testing as it will invalidate most
// time instants and other time manager functions
void Calendar::reset() {

   LOG_WARN("Removing Omega calendar");
   LOG_WARN("This invalidates all prior time manager values"
            "and should only be used in testing");
   Calendar::OmegaCal.reset(); // Resets the calendar pointer
}

// Calendar operators
//------------------------------------------------------------------------------
// Calendar::isDefined - checks whether a calendar has been initialized
bool Calendar::isDefined() {
   if (Calendar::OmegaCal) {
      return true;
   } else {
      return false;
   }
}

//------------------------------------------------------------------------------
// Calendar::validate - validate Calendar state
// Validates the Calendar state, checking for improper values.

I4 Calendar::validate() const {

   // initialize error code
   I4 Err{0};

   // Check if calendar kind was invalid
   if (this->CalKind == CalendarUnknown) {
      Err = 1;
      LOG_ERROR("TimeMgr: Calendar::validate unknown calendar kind");
      return Err;
   }

   // Check type of calendar
   if (this->CalKind < 1 || this->CalKind > NUM_SUPPORTED_CALENDARS) {
      Err = 1;
      LOG_ERROR("TimeMgr: Calendar::validate calendar kind out of range");
      return Err;
   }

   // Check seconds per day
   if (this->SecondsPerDay < 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: Calendar::validate invalid seconds per day");
      return Err;
   }

   // Check seconds per year
   if (this->SecondsPerYear < 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: Calendar::validate invalid seconds per year");
      return Err;
   }

   // Check days per year
   if (this->DaysPerYear < 0) {
      Err = 1;
      LOG_ERROR("TimeMgr: Calendar::validate invalid days per year");
      return Err;
   }

   // Check days per month
   if (this->MonthsPerYear > 0 && this->MonthsPerYear <= MONTHS_PER_YEAR) {
      I4 CountDaysPerYr{0};
      for (I4 I = 0; I < this->MonthsPerYear; I++) {
         CountDaysPerYr += this->DaysPerMonth[I];
         if (this->DaysPerMonth[I] < 0) {
            Err = 1;
            LOG_ERROR("TimeMgr: Calendar::validate invalid days per month "
                      "DaysPerMonth[{}] = {}",
                      I, this->DaysPerMonth[I]);
            return Err;
         }
      }

      // check whether sum of days per month is same as days per year
      if (this->DaysPerYear != CountDaysPerYr) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::validate DaysPerYear != sum of "
                   "DaysPerMonth[]");
         return Err;
      }
   }

   // all done, return error code
   return Err;

} // end Calendar::validate

//------------------------------------------------------------------------------
// Calendar::isLeapYear - Determine if given year is a leap year
// Logical function that returns true if input year is a leap year.

bool Calendar::isLeapYear(I8 Year // [in] a calendar year
) {

   // Check calendar existence
   if (!isDefined()) {
      LOG_CRITICAL("Attempt to check leap year - calendar not defined");
      return false;
   }

   // now check for leap year depending on calendar kind
   switch (Calendar::OmegaCal->CalKind) {
   case CalendarGregorian:
      // leap year is divisible by 400 or divisible by 4 and not 100.
      return (Year % 400 == 0) || ((Year % 4 == 0) && (Year % 100 != 0));
      break;

   case CalendarJulian:
      // leap year is divisible by 4.
      return Year % 4 == 0;
      break;

   default:
      // all other calendars don't have leap years.
      return false;
      break;
   } // end switch CalKind

} // end Calendar::isLeapYear

//------------------------------------------------------------------------------
// Calendar::getElapsedTime - return time since ref time for a given date/time
// Retrieves the total elapsed time (seconds, in TimeFrac form) since the
// calendar reference time associated with a particular calendar date and
// time of day.

TimeFrac Calendar::getElapsedTime(
    const I8 Year,   ///< [in] calendar year
    const I8 Month,  ///< [in] calendar month
    const I8 Day,    ///< [in] calendar day
    const I8 Hour,   ///< [in] time of day-hour
    const I8 Minute, ///< [in] time of day-min
    const I8 Whole,  ///< [in] time of day-whole seconds
    const I8 Numer,  ///< [in] time of day-frac secs (numerator)
    const I8 Denom   ///< [in] time of day-frac secs (denom)
) {

   // initialize error code, the basetime result, and common temps
   I4 Err{0};
   TimeFrac Result(0, 0, 1);
   I8 ResultWhole{0}; // whole seconds for result
   I8 DayOfYear{0};   // day since beginning of year
   I8 JD{0};          // Julian Day used for many conversions
   I8 HourTmp{0};     // For half-day JD conversions
   I4 FebDays{0};     // For tracking leap days

   // Check that calendar has been defined
   if (!isDefined()) {
      LOG_CRITICAL("Cannot get elapsed time - calendar not defined");
      return Result;
   }

   // Branch on calendar type for actual calculation
   switch (Calendar::OmegaCal->CalKind) {

   // convert Gregorian Date to time since reference of noon 3/1/-4800
   case CalendarGregorian: {
      // Validate inputs
      if (Year < -4800 || Month < 1 || Month > 12 || Day < 1) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid Gregorian date");
         break;
      }
      // invalid before 3/1/-4800
      if (Year == -4800 && Month < 3) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid date prior to "
                   "reference Gregorian date");
         break;
      }
      // check day of the month for any month except February
      if (Month != 2 && Day > Calendar::OmegaCal->DaysPerMonth[Month - 1]) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid day of month for "
                   "Gregorian calendar");
         break;
      }
      // if February, take leap year into account before checking
      //   day of the month
      if (Month == 2) {
         FebDays = Calendar::OmegaCal->DaysPerMonth[1];
         if (isLeapYear(Year))
            FebDays += 1;
         if (Day > FebDays) {
            Err = 1;
            LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid day of "
                      "February for Gregorian calendar");
            break;
         }
      }
      // Check valid hour, minute, seconds
      if (Hour < 0 || Hour > 24 || Minute < 0 || Minute > 60 || Whole < 0 ||
          Whole > SECONDS_PER_MINUTE) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid time of day "
                   "(hour, min, sec) for Gregorian calendar");
         break;
      }
      // convert Gregorian date to Julian days
      // using Fliegel, van Flandern algorithm
      // Gregorian date (year, month, day) => Julian days (JD)
      I4 Temp = (Month - 14) / 12;
      JD      = (1461 * (Year + 4800 + Temp)) / 4 +
           (367 * (Month - 2 - 12 * Temp)) / 12 -
           (3 * ((Year + 4900 + Temp) / 100)) / 4 + Day - 32075;

      // Julian Day starts at noon, so correct for the half day
      HourTmp = Hour - 12;
      if (HourTmp < 0) {
         HourTmp += 24;
         JD -= 1;
      }

      // Finally, convert JD to seconds and add hours, minutes, secs
      ResultWhole = JD * Calendar::OmegaCal->SecondsPerDay +
                    HourTmp * SECONDS_PER_HOUR + Minute * SECONDS_PER_MINUTE +
                    Whole;
      Err = Result.set(ResultWhole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime Gregorian - error "
                   "setting final result");
      break;
   }

   // The Julian <-> Julian day conversion algorithm is from D.A. Hatcher,
   // Q.JlR. astr. Soc. (1984) 25, 53-55.  It is valid from 3/1/-4712 forward.
   case CalendarJulian: {
      // Validate inputs
      if (Year < -4712 || Month < 1 || Month > 12 || Day < 1) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid Julian date");
         break;
      }
      if (Year == -4712 && Month < 3) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid date prior to "
                   "reference Julian date");
         break;
      }
      // check day of the month for any month except February
      if (Month != 2 && Day > Calendar::OmegaCal->DaysPerMonth[Month - 1]) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid day of month for "
                   "Julian calendar");
         break;
      }
      // if February, take leap year into account before checking
      //   day of the month
      if (Month == 2) {
         FebDays = Calendar::OmegaCal->DaysPerMonth[1];
         if (isLeapYear(Year))
            FebDays += 1;
         if (Day > FebDays) {
            Err = 1;
            LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid day of "
                      "February for Julian calendar");
            break;
         }
      }
      // Check valid hour, minute, seconds
      if (Hour < 0 || Hour > 24 || Minute < 0 || Minute > 60 || Whole < 0 ||
          Whole > SECONDS_PER_MINUTE) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid time of day "
                   "(hour, min, sec) for Julian calendar");
         break;
      }

      // Convert date to Julian day number
      I8 YPrime = Year - ((12 - Month) / 10);
      I4 MPrime = (Month + 9) % 12;
      I8 Y      = (I8)(365.25 * (YPrime + 4712));
      I4 D      = (I4)((30.6 * MPrime) + 0.5);
      JD        = Y + D + Day + 59;

      // Julian Day starts at noon, so correct for the half day
      HourTmp = Hour - 12;
      if (HourTmp < 0) {
         HourTmp += 24;
         JD -= 1;
      }

      // Finally, conver JD to seconds and add hours, minutes, secs
      ResultWhole = JD * Calendar::OmegaCal->SecondsPerDay +
                    HourTmp * SECONDS_PER_HOUR + Minute * SECONDS_PER_MINUTE +
                    Whole;
      Err = Result.set(ResultWhole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime Julian - error setting "
                   "final result");
      break;
   }

   // For Julian Day calendars, just counting days from ref time
   case CalendarJulianDay:
   case CalendarModJulianDay: {
      // Validate inputs
      if (Year != 0 || Month != 0 || Day < 1) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime Invalid date for Julian "
                   "Day calendars");
         break;
      }
      // Check valid hour, minute, seconds
      if (Hour < 0 || Hour > 24 || Minute < 0 || Minute > 60 || Whole < 0 ||
          Whole > SECONDS_PER_MINUTE) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid time of day "
                   "(hour, min, sec) for Julian calendar");
      }
      // Elapsed time just seconds per day plus time of day
      // For regular Julian Day, the start of the day is at noon on
      // other calendars - does not matter for internal calendar
      // operations, so no correction and just noted here
      ResultWhole = (Day - 1) * Calendar::OmegaCal->SecondsPerDay +
                    Hour * SECONDS_PER_HOUR + Minute * SECONDS_PER_MINUTE +
                    Whole;
      // Now set final result and add fractional second
      Err = Result.set(ResultWhole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime Julian Day calendars - "
                   "error setting final result");
      break;
   }

   case CalendarNoLeap:
   case Calendar360Day:
   case CalendarCustom: {
      // Validate inputs
      if (Year < 0 || Month < 1 || Month > Calendar::OmegaCal->MonthsPerYear ||
          Day < 1) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid date for "
                   "fixed-length calendars");
         break;
      }
      // check day of the month for any month except February
      if (Day > Calendar::OmegaCal->DaysPerMonth[Month - 1]) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid date of month "
                   "for fixed-length calendars");
         break;
      }
      // Check valid hour, minute, seconds
      if (Hour < 0 || Hour > 24 || Minute < 0 || Minute > 60 || Whole < 0 ||
          Whole > SECONDS_PER_MINUTE) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid time of day "
                   "(hour, min, sec) for fixed-length calendars");
         break;
      }

      // Conversion straightforward for fixed-length calendars
      ResultWhole = Year * Calendar::OmegaCal->SecondsPerYear; // at beg of year
      for (I4 Imonth = 0; Imonth < Month - 1; Imonth++)
         ResultWhole += Calendar::OmegaCal->DaysPerMonth[Imonth] *
                        Calendar::OmegaCal->SecondsPerDay; // add months
      ResultWhole += (Day - 1) * Calendar::OmegaCal->SecondsPerDay; // add days
      ResultWhole +=
          Hour * SECONDS_PER_HOUR + Minute * SECONDS_PER_MINUTE + Whole;
      // Now set final result and add fractional second
      Err = Result.set(ResultWhole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime fixed-length calendars - "
                   "error setting final result");
      break;
   }

   // No calendar option should only be tracking seconds, so check for
   // inappropriate entries for other fields and simply return the
   // input seconds field
   case CalendarNoCalendar: {
      // Validate inputs
      if (Year != 0 || Month != 0 || Day != 0 || Hour != 0 || Minute != 0)
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime NoCalendar option should "
                   "not be tracking date or time of day");
      Err = Result.set(Whole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime NoCalendar calendars - "
                   "error setting final result");
      break;
   }

   // Undefined calendars
   case CalendarUnknown:
   default: {
      LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid calendar");
      break;
   }
   } // end switch CalKind

   return Result;

} // end Calendar::getElapsedTime(

//------------------------------------------------------------------------------
// Calendar::getDateTime - Determine date and time of day given elapsed time
// Determines the calendar date and time of day, given an elapsed time since
// the calendar reference time.

I4 Calendar::getDateTime(
    const TimeFrac ElapsedTime, //< [in] time (secs) from ref time
    I8 &Year,                   //< [out] calendar year
    I8 &Month,                  //< [out] calendar month
    I8 &Day,                    //< [out] calendar day
    I8 &Hour,                   //< [out] time of day-hours
    I8 &Minute,                 //< [out] time of day-minutes
    I8 &Whole,                  //< [out] time of day-whole seconds
    I8 &Numer,                  //< [out] time of day-frac secs (numerator)
    I8 &Denom                   //< [out] time of day-frac secs (denom)
) {

   // initialize error code to success and make sure input elapsed
   // time is reduced to its simplest form
   I4 Error{0};

   // Check existence of calendar
   if (!isDefined()) {
      LOG_CRITICAL("Attempt to getDateTime on non-existent calendar");
      return 1;
   }

   TimeFrac RevisedTime = ElapsedTime;
   Error                = RevisedTime.simplify();
   if (Error != 0) {
      LOG_ERROR("TimeMgr: Calendar::getDateTime error simplifying time "
                "representation");
      return Error;
   }

   // retrieve whole and fractional seconds to manipulate
   Error = RevisedTime.get(Whole, Numer, Denom);
   if (Error != 0) {
      LOG_ERROR("TimeMgr: Calendar::getDateTime error retrieving fractional "
                "second components");
      return Error;
   }

   // initialize common counters used by several calendars
   I4 DayInYear{0};
   I4 CountDays{0};
   I4 PrevDays{0};
   I8 JD{0};

   // branch to appropriate calculation based on calendar type
   switch (OmegaCal->CalKind) {
   case CalendarGregorian: {
      // Reference time same as Julian Days
      // Convert whole seconds to Julian Day (number of days)
      // and remove that from elapsed seconds
      JD = Whole / Calendar::OmegaCal->SecondsPerDay;
      Whole -= (JD * Calendar::OmegaCal->SecondsPerDay);

      // Now convert the remaining seconds to hour, minute, seconds
      Hour = Whole / SECONDS_PER_HOUR;
      Whole -= Hour * SECONDS_PER_HOUR;
      Minute = Whole / SECONDS_PER_MINUTE;
      Whole -= Minute * SECONDS_PER_MINUTE;

      // Since Julian date uses a noon time for start of day, must
      // correct to a midnight start and modify day count if necessary
      Hour += 12;
      if (Hour >= 24) {
         Hour -= 24;
         JD += 1;
      }

      // convert Julian days to Gregorian date
      // Julian days (jdays) => Gregorian date (yy, mm, dd)

      I8 TempL = JD + 68569;
      I8 TempN = (4 * TempL) / 146097;
      TempL    = TempL - (146097 * TempN + 3) / 4;
      I8 TempI = (4000 * (TempL + 1)) / 1461001;
      TempL    = TempL - (1461 * TempI) / 4 + 31;
      I8 TempJ = (80 * TempL) / 2447;

      Day   = TempL - (2447 * TempJ) / 80;
      TempL = TempJ / 11;
      Month = TempJ + 2 - (12 * TempL);
      Year  = 100 * (TempN - 49) + TempI + TempL;

      // All done with Gregorian calendar
      break;
   }

   // The Julian <-> Julian day conversion algorithm is from D.A. Hatcher,
   // Q.JlR. astr. Soc. (1984) 25, 53-55.  It is valid from 3/1/-4712 forward.
   case CalendarJulian: {
      // Reference time same as Julian Days
      // Convert whole seconds to Julian Day (number of days)
      // and remove that from elapsed seconds
      JD = Whole / Calendar::OmegaCal->SecondsPerDay;
      Whole -= JD * Calendar::OmegaCal->SecondsPerDay;

      // Now convert the remaining seconds to hour, minute, seconds
      Hour = Whole / SECONDS_PER_HOUR;
      Whole -= Hour * SECONDS_PER_HOUR;
      Minute = Whole / SECONDS_PER_MINUTE;
      Whole -= Minute * SECONDS_PER_MINUTE;

      // Since Julian date uses a noon time for start of day, must
      // correct to a midnight start and modify day count if necessary
      Hour += 12;
      if (Hour >= 24) {
         Hour -= 24;
         JD += 1;
      }

      // Julian day (D) => Julian date (yy, mm, dd)
      Year      = (I8)((R8)JD / 365.25) - 4712;
      I8 DPrime = (I8)fmod((JD - 59.25), 365.25);
      Month     = ((I8)((DPrime + 0.5) / 30.6) + 2) % 12 + 1;
      Day       = (I8)(fmod((DPrime + 0.5), 30.6)) + 1;

      break;
   }

   case CalendarJulianDay:
   case CalendarModJulianDay: {
      // Since Julian day calendars just count days since ref time
      // year and month not defined
      // Convert whole seconds to Julian Day (number of days)
      // and remove that from elapsed seconds
      JD = Whole / Calendar::OmegaCal->SecondsPerDay;
      Whole -= JD * Calendar::OmegaCal->SecondsPerDay;
      Year  = 0;
      Month = 0;
      Day   = JD + 1; // correct for 1-based counting

      // Now convert the remaining seconds to hour, minute, seconds
      // Note that for JD calendar we retain the convention of a noon
      // start for the day and do not correct the hour.
      Hour = Whole / SECONDS_PER_HOUR;
      Whole -= Hour * SECONDS_PER_HOUR;
      Minute = Whole / SECONDS_PER_MINUTE;
      Whole -= Minute * SECONDS_PER_MINUTE;
      break;
   }

   // for calendars with fixed number of days and year 0000 reference
   case CalendarNoLeap:
   case Calendar360Day:
   case CalendarCustom: {
      Year = Whole / Calendar::OmegaCal->SecondsPerYear;      // determine year
      Whole -= Calendar::OmegaCal->SecondsPerYear * Year;     // holds remainder
      DayInYear = Whole / Calendar::OmegaCal->SecondsPerDay;  // day in year
      Whole -= Calendar::OmegaCal->SecondsPerDay * DayInYear; // holds remainder
      DayInYear += 1; // correct for 1-based counting
      // find month, day
      CountDays = 0;
      PrevDays  = 0;
      for (I4 I = 0; I < Calendar::OmegaCal->MonthsPerYear; I++) {
         CountDays += Calendar::OmegaCal->DaysPerMonth[I];
         if (DayInYear <= CountDays) { // day is in this month
            Month = I + 1;
            Day   = DayInYear - PrevDays;
            break;
         }
         PrevDays += Calendar::OmegaCal->DaysPerMonth[I];
      }
      // keep peeling off for hour, minute
      Hour = Whole / SECONDS_PER_HOUR;
      Whole -= Hour * SECONDS_PER_HOUR;
      Minute = Whole / SECONDS_PER_MINUTE;
      Whole -= Minute * SECONDS_PER_MINUTE;
      break;
   }

   case CalendarNoCalendar: {
      // If no calendar is used, we only track elapsed time so simply
      // return elapsed time in current form and zero all other elements.
      Year   = 0;
      Month  = 0;
      Day    = 0;
      Hour   = 0;
      Minute = 0;
      break;
   }

   case CalendarUnknown:
   default: {
      Error = 1;
      LOG_ERROR("TimeMgr: Calendar::getDateTime invalid calendar");
      break;
   }
   } // end switch calKindFlag

   return Error;

} // end Calendar::getDateTime

//------------------------------------------------------------------------------
// Calendar::incrementDate - Increments (decrements) date by a given interval
// Increments (or decrements if negative) the date by a given integer time
// interval. Only calendar intervals are supported.

I4 Calendar::incrementDate(
    const I8 Interval,     //< [in] interval to advance time
    const TimeUnits Units, //< [in] units for input interval
    I8 &Year,              //< [in,out] calendar year for time to be advanced
    I8 &Month,             //< [in,out] calendar month for time to be advanced
    I8 &Day                //< [in,out] calendar day for time to be advanced
) {

   // initialize error code
   I4 Err{0};

   // Check calendar exists
   if (!isDefined()) {
      LOG_CRITICAL("Attempt to increment non-existent calendar");
      return 1;
   }

   // Increment date based on supported interval types
   switch (Units) {
   // For year intervals, mostly just increment/decrement the year
   // while keeping all other units fixed, though check for leap day error
   case TimeUnits::Years: {
      switch (Calendar::OmegaCal->CalKind) {
      case CalendarGregorian:
      case CalendarNoLeap:
      case CalendarJulian:
      case Calendar360Day:
      case CalendarCustom: {
         Year += Interval;
         if (!(isLeapYear(Year)) && Month == 2 && Day == 29) {
            Err = 1;
            LOG_ERROR("TimeMgr: Calendar::incrementDate day out of range for "
                      "new year increment");
         }
         break;
      }
      // all other options return an error since year undefined
      case CalendarJulianDay:
      case CalendarModJulianDay:
      case CalendarNoCalendar:
      case CalendarUnknown:
      default: {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::incrementDate invalid calendar for year "
                   "interval");
         break;
      }
      } // end switch on calendar type
      break;
   }

   // For monthly intervals, check for year rollovers
   case TimeUnits::Months: {
      switch (Calendar::OmegaCal->CalKind) {
      // For most calendars, add the monthly interval and adjust year
      // and day accordingly.
      case CalendarGregorian:
      case CalendarNoLeap:
      case CalendarJulian:
      case Calendar360Day:
      case CalendarCustom: {
         I8 TmpMonth = Month + Interval;
         // correct the year if the interval pushes beyond the end of year
         // use a while loop in case interval extends across multiple years
         while (TmpMonth > Calendar::OmegaCal->MonthsPerYear) {
            Year += 1;
            TmpMonth -= Calendar::OmegaCal->MonthsPerYear;
         }
         // For negative intervals, check the opposite direction
         while (TmpMonth < 1) {
            Year -= 1;
            TmpMonth += Calendar::OmegaCal->MonthsPerYear;
         }
         // Set month and check that the day of the month is valid
         Month      = TmpMonth;
         I8 MaxDays = Calendar::OmegaCal->DaysPerMonth[Month - 1];
         if (isLeapYear(Year) && Month == 2)
            MaxDays += 1;
         if (Day < 1 || Day > MaxDays) {
            Err = 1;
            LOG_ERROR("TimeMgr: Calendar::incrementDate day outside range for "
                      "new month");
         }
         break;
      }
      // all other options return an error since month undefined
      case CalendarJulianDay:
      case CalendarModJulianDay:
      case CalendarNoCalendar:
      case CalendarUnknown:
      default: {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::incrementDate invalid calendar for "
                   "month interval");
         break;
      }
      } // end switch on calendar type
      break;
   }

   // For day intervals, increment days and adjust year, month
   case TimeUnits::Days: {
      switch (Calendar::OmegaCal->CalKind) {
      case CalendarGregorian:
      case CalendarNoLeap:
      case CalendarJulian:
      case Calendar360Day:
      case CalendarCustom: {
         if (Interval > 0) { // positive interval - step forward
            I4 DayMax = Calendar::OmegaCal->DaysPerMonth[Month - 1];
            if (Month == 2 && isLeapYear(Year))
               DayMax += 1;
            for (I4 I = 1; I <= Interval; ++I) {
               Day += 1;
               if (Day > DayMax) { // rollover month boundary
                  Day = 1;
                  Month += 1;
                  // rollover year boundary
                  if (Month > Calendar::OmegaCal->MonthsPerYear) {
                     Month = 1;
                     Year += 1;
                  }
                  DayMax = Calendar::OmegaCal->DaysPerMonth[Month - 1];
                  if (Month == 2 && isLeapYear(Year))
                     DayMax += 1;
               }
            }
         } else { // negative interval - step backward
            for (I4 I = 1; I <= llabs(Interval); ++I) {
               Day -= 1;
               if (Day < 1) {
                  Month -= 1;
                  if (Month < 1) {
                     Month = Calendar::OmegaCal->MonthsPerYear;
                     Year -= 1;
                  }
                  Day = Calendar::OmegaCal->DaysPerMonth[Month - 1];
                  if (Month == 2 && isLeapYear(Year))
                     Day += 1;
               }
            }
         }
         break;
      }
      // For day-only calendars, just increment/decrement the day
      case CalendarJulianDay:
      case CalendarModJulianDay: {
         Day += Interval;
         break;
      }
      // all other options return an error
      case CalendarNoCalendar:
      case CalendarUnknown:
      default: {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::incrementDate invalid calendar for day "
                   "interval");
         break;
      }
      } // end switch on calendar type
      break;
   }

   // All other intervals return an error
   case TimeUnits::None:
   case TimeUnits::Seconds:
   case TimeUnits::Minutes:
   case TimeUnits::Hours:
   default: {
      Err = 1;
      LOG_ERROR("TimeMgr: Calendar::incrementDateTime only calendar intervals "
                "supported");
      break;
   }
   } // end switch on interval type

   // All done, return the results
   return Err;

} // end Calendar::incrementDateTime

//------------------------------------------------------------------------------
// TimeInterval definitions
//------------------------------------------------------------------------------

// TimeInterval constructors/destructors
//------------------------------------------------------------------------------
// TimeInterval::TimeInterval - default constructor
// This is a default constructor that creates a zero, non-calendar interval.

TimeInterval::TimeInterval(void) {

   // Create a base time interval of zero
   Interval = TimeFrac(0, 0, 1);
   Units    = TimeUnits::Seconds;

   // Default is not a calendar interval
   IsCalendar  = false;
   CalInterval = 0;

} // end TimeInterval::TimeInterval

//------------------------------------------------------------------------------
// TimeInterval::TimeInterval
// This is a time interval constructor that creates a non-calendar interval from
// the integer fraction representation (full + numerator + denominator).

TimeInterval::TimeInterval(const I8 Whole, // [in] whole integer seconds
                           const I8 Numer, // [in] fractional seconds numerator
                           const I8 Denom // [in] fractional seconds denominator
) {

   // Create non-calendar TimeFrac interval using TimeFrac constructor
   Interval = TimeFrac(Whole, Numer, Denom);
   Units    = TimeUnits::Seconds;

   // Default is not a calendar interval
   IsCalendar  = false;
   CalInterval = 0;

} // end TimeInterval::TimeInterval frac second constructor

//------------------------------------------------------------------------------
// TimeInterval::TimeInterval
// Constructs a time interval from an integer length and specified unit of time.

TimeInterval::TimeInterval(
    const I4 InLength,      // [in] Time interval length
    const TimeUnits InUnits // [in] Units of time for this interval
) {

   // Define a default time interval
   Interval    = TimeFrac(0, 0, 1);
   Units       = InUnits;
   IsCalendar  = false;
   CalInterval = 0;
   I4 Err{0};

   // Now set values based on input units
   I8 Length = InLength;
   switch (Units) {
   case TimeUnits::Seconds:
      Err = Interval.set(Length, 0, 1);
      break;
   case TimeUnits::Minutes:
      Err = Interval.setHMS(0, Length, 0);
      break;
   case TimeUnits::Hours:
      Err = Interval.setHMS(Length, 0, 0);
      break;
   case TimeUnits::Years:  // these three are all calendar
   case TimeUnits::Months: // intervals with the input length
   case TimeUnits::Days:   // units already set correctly
      IsCalendar  = true;
      CalInterval = Length;
      break;
   case TimeUnits::None:
   default:
      LOG_ERROR("TimeMgr: invalid time units for TimeInterval constructor");
   }

} // end TimeInterval::TimeInterval I4 length/unit constructor

//------------------------------------------------------------------------------
// TimeInterval::TimeInterval
// Constructs a time interval from an integer length and specified unit of time.

TimeInterval::TimeInterval(
    const I8 Length,        // [in] Time interval length
    const TimeUnits InUnits // [in] Units of time for this interval
) {

   // Define a default time interval
   Interval    = TimeFrac(0, 0, 1);
   Units       = InUnits;
   IsCalendar  = false;
   CalInterval = 0;
   I4 Err{0};

   // Now set values based on input units
   switch (InUnits) {
   case TimeUnits::Seconds:
      Err = Interval.set(Length, 0, 1);
      break;
   case TimeUnits::Minutes:
      Err = Interval.setHMS(0, Length, 0);
      break;
   case TimeUnits::Hours:
      Err = Interval.setHMS(Length, 0, 0);
      break;
   case TimeUnits::Years:  // these three are all calendar
   case TimeUnits::Months: // intervals with the input length
   case TimeUnits::Days:   // units already set correctly
      IsCalendar  = true;
      CalInterval = Length;
      break;
   case TimeUnits::None:
   default:
      LOG_ERROR("TimeMgr: invalid time units for TimeInterval constructor");
   }

} // end TimeInterval::TimeInterval I8 length/unit constructor

//------------------------------------------------------------------------------
// TimeInterval::TimeInterval
// Constructs a time interval from a real length and specified unit of time.
// If a calendar interval is requested the real input is converted to an
// integer number of months, days or years using standard conversion.

TimeInterval::TimeInterval(
    const R8 Length,        // [in] Time interval length
    const TimeUnits InUnits // [in] Units of time for this interval
) {

   // Define a default time interval
   Interval    = TimeFrac(0, 0, 1);
   Units       = InUnits;
   IsCalendar  = false;
   CalInterval = 0;
   I4 Err{0};

   // Now set values based on input units
   switch (InUnits) {
   case TimeUnits::Seconds:
      Err = Interval.setSeconds(Length); // TimeFrac set
      break;
   case TimeUnits::Minutes:
      Err = Interval.setMinutes(Length); // TimeFrac set
      break;
   case TimeUnits::Hours:
      Err = Interval.setHours(Length); // TimeFrac set
      break;
   case TimeUnits::Years:  // these three are all calendar
   case TimeUnits::Months: // intervals with the input length
   case TimeUnits::Days:   // units already set correctly
      IsCalendar  = true;
      CalInterval = Length; // Length is coverted to an integer
      break;
   case TimeUnits::None:
   default:
      LOG_ERROR("TimeMgr: invalid time units for TimeInterval constructor");
   }

} // end TimeInterval::TimeInterval real length/unit constructor

//------------------------------------------------------------------------------
// TimeInterval::TimeInterval
// Construct a time interval from a standard string in the form
// DDDD_HH:MM:SS.SSSS where the width of DD and SS strings can be of
// arbitrary width (within reason) and the separators can be any single
// non-numeric character
TimeInterval::TimeInterval(
    std::string &TimeString // [in] string form of time interval
) {

   // Not a calendar interval
   IsCalendar  = false;
   CalInterval = 0;

   I4 Err{0};

   // Extract variables from string
   I8 Day     = 0;
   I8 Hour    = 0;
   I8 Minute  = 0;
   R8 RSecond = 0.;

   std::istringstream ss(TimeString);
   char discard;
   ss >> Day >> discard >> Hour >> discard >> Minute >> discard >> RSecond;

   R8 SecondsSet = Day * SECONDS_PER_DAY + Hour * SECONDS_PER_HOUR +
                   Minute * SECONDS_PER_MINUTE + RSecond;

   // Use the set function to define the Time Interval
   Err = this->set(SecondsSet, TimeUnits::Seconds);

   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInterval constructor from string error "
                "using using set function to construct TimeInterval.");

} // end TimeInterval::TimeInterval constructor from string

//------------------------------------------------------------------------------
// TimeInterval::~TimeInterval - destructor for time interval
// Destructor for a time interval. No allocated space so does nothing.

TimeInterval::~TimeInterval(void) { // Nothing to be done

} // end TimeInterval destructor

// TimeInterval accessors
//------------------------------------------------------------------------------
// TimeInterval::set
// Set a non-calendar time interval from the base time integer fraction
// representation (whole seconds + numerator + denominator).

I4 TimeInterval::set(const I8 Whole, // [in] whole integer seconds
                     const I8 Numer, // [in] fractional seconds numerator
                     const I8 Denom  // [in] fractional seconds denominator
) {

   // Create non-calendar TimeFrac interval using TimeFrac constructor
   Interval = TimeFrac(Whole, Numer, Denom);
   Units    = TimeUnits::Seconds;

   // Default is not a calendar interval
   IsCalendar  = false;
   CalInterval = 0;

   return 0;

} // end TimeInterval::set (fractional seconds)

//------------------------------------------------------------------------------
// TimeInterval::set
// Sets a time interval from an integer length and specified unit of time.

I4 TimeInterval::set(
    const I4 InLength,      // [in] Time interval length
    const TimeUnits InUnits // [in] Units of time for this interval
) {

   // Reset to a default time interval
   Interval    = TimeFrac(0, 0, 1);
   Units       = InUnits;
   IsCalendar  = false;
   CalInterval = 0;
   I4 Err{0};

   // Now set values based on input units
   I8 Length = InLength;
   switch (InUnits) {
   case TimeUnits::Seconds:
      Err = Interval.set(Length, 0, 1);
      break;
   case TimeUnits::Minutes:
      Err = Interval.setHMS(0, Length, 0);
      break;
   case TimeUnits::Hours:
      Err = Interval.setHMS(Length, 0, 0);
      break;
   case TimeUnits::Years:  // these three are all calendar
   case TimeUnits::Months: // intervals with the input length
   case TimeUnits::Days:   // units already set correctly
      IsCalendar  = true;
      CalInterval = Length;
      break;
   case TimeUnits::None:
   default:
      Err = 1;
      LOG_ERROR("TimeMgr: invalid time units for TimeInterval::set");
   }

   return Err;

} // end TimeInterval::set (I4 length/unit)

//------------------------------------------------------------------------------
// TimeInterval::set
// Sets a time interval from an integer length and specified unit of time.

I4 TimeInterval::set(
    const I8 Length,        // [in] Time interval length
    const TimeUnits InUnits // [in] Units of time for this interval
) {

   // Reset to a default time interval
   Interval    = TimeFrac(0, 0, 1);
   Units       = InUnits;
   IsCalendar  = false;
   CalInterval = 0;
   I4 Err{0};

   // Now set values based on input units
   switch (InUnits) {
   case TimeUnits::Seconds:
      Err = Interval.set(Length, 0, 1);
      break;
   case TimeUnits::Minutes:
      Err = Interval.setHMS(0, Length, 0);
      break;
   case TimeUnits::Hours:
      Err = Interval.setHMS(Length, 0, 0);
      break;
   case TimeUnits::Years:  // these three are all calendar
   case TimeUnits::Months: // intervals with the input length
   case TimeUnits::Days:   // units already set correctly
      IsCalendar  = true;
      CalInterval = Length;
      break;
   case TimeUnits::None:
   default:
      Err = 1;
      LOG_ERROR("TimeMgr: invalid time units for TimeInterval::set");
   }

   return Err;

} // end TimeInterval::set (I8 length/unit)

//------------------------------------------------------------------------------
// TimeInterval::set
// Sets a time interval from a real length and specified unit of time.
// Because a calendar interval only supports integer lengths, the input
// real value is converted to an integer number of months, days or years
// using standard conversions.

I4 TimeInterval::set(
    const R8 Length,        // [in] Time interval length
    const TimeUnits InUnits // [in] Units of time for this interval
) {

   // Reset to a default time interval
   Interval    = TimeFrac(0, 0, 1);
   Units       = InUnits;
   IsCalendar  = false;
   CalInterval = 0;
   I4 Err{0};

   // Now set values based on input units
   switch (InUnits) {
   case TimeUnits::Seconds:
      Err = Interval.setSeconds(Length); // TimeFrac set
      break;
   case TimeUnits::Minutes:
      Err = Interval.setMinutes(Length); // TimeFrac set
      break;
   case TimeUnits::Hours:
      Err = Interval.setHours(Length); // TimeFrac set
      break;
   case TimeUnits::Years:  // these three are all calendar
   case TimeUnits::Months: // intervals with the input length
   case TimeUnits::Days:   // units already set correctly
      IsCalendar  = true;
      CalInterval = Length; // length is coverted to an integer
      break;
   case TimeUnits::None:
   default:
      Err = 1;
      LOG_ERROR("TimeMgr: invalid time units for TimeInterval::set");
   }

   return Err;

} // end TimeInterval::set (real length/unit)

//------------------------------------------------------------------------------
// TimeInterval::get
// Retrieve a non-calendar time interval in fractional integer seconds

I4 TimeInterval::get(I8 &Whole, ///< [out] whole seconds
                     I8 &Numer, ///< [out] fractional second numerator
                     I8 &Denom  ///< [out] fractional second denominator
) const {

   I4 Err{0};

   // If this is a calendar interval, return an error
   if (IsCalendar) {
      Err = 1;
      LOG_ERROR("TimeMgr: TimeInterval::get attempt to retrieve non-calendar "
                "values from calendar interval");
   } else {
      // otherwise, just use the TimeFrac::get to extract
      Err = Interval.get(Whole, Numer, Denom);
   }

   return Err;

} // end TimeInterval::get (fractonal seconds)

//------------------------------------------------------------------------------
// TimeInterval::get (integer)
// Retrieve a time interval in requested units and integer length.
// Due to issues in unit conversion, we only allow integer retrieval to
// return in the same units that the interval was defined.

I4 TimeInterval::get(
    I8 &Length,              // requested length of interval
    const TimeUnits ReqUnits // units requested for time interval
) const {

   Length = 0;
   I4 Hours{0}; // needed for TimeFrac HMS interface
   I4 Minutes{0};
   I4 Seconds{0};
   I4 Err{0};

   // Now get values based on requested units
   // If requested units are same as the interval was defined,
   // just return length. Otherwise, attempt to convert or
   // write an error if not possible.
   switch (ReqUnits) {
   case TimeUnits::Seconds:
      Err = Interval.getHMS(Hours, Minutes, Seconds); // TimeFrac get
      if (Units == TimeUnits::Seconds) {
         Length = Seconds;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (integer) attempt to get "
                   "interval in seconds but defined in other units");
      }
      break;
   case TimeUnits::Minutes:
      Err = Interval.getHMS(Hours, Minutes, Seconds); // TimeFrac get
      if (Units == TimeUnits::Minutes) {
         Length = Minutes;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (integer) attempt to get "
                   "interval in minutes but defined in other units");
      }
      break;
   case TimeUnits::Hours:
      Err = Interval.getHMS(Hours, Minutes, Seconds); // TimeFrac get
      if (Units == TimeUnits::Hours) {
         Length = Hours;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (integer) attempt to get "
                   "interval in hours but defined in other units");
      }
      break;
   case TimeUnits::Years:
      if (Units == TimeUnits::Years) {
         Length = CalInterval;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (integer) attempt to get "
                   "interval in years but defined in other units");
      }
      break;
   case TimeUnits::Months:
      if (Units == TimeUnits::Months) {
         Length = CalInterval;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (integer) attempt to get "
                   "interval in months but defined in other units");
      }
      break;
   case TimeUnits::Days:
      if (Units == TimeUnits::Days) {
         Length = CalInterval;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (integer) attempt to get "
                   "interval in days but defined in other units");
      }
      break;
   case TimeUnits::None:
   default:
      Err = 1;
      LOG_ERROR("TimeMgr: TimeInterval::get (integer) invalid time units for "
                "time interval");
   }

   // return error code
   return Err;

} // end TimeInterval::get (integer)

//------------------------------------------------------------------------------
// TimeInterval::get (real)
// Retrieve a time interval in requested units and real length.
// For non-calendar intervals, this routine will convert units to desired
// units. Because year, month intervals change based on current date, no
// conversion is supported yet.

I4 TimeInterval::get(
    R8 &Length,              // requested length of interval
    const TimeUnits ReqUnits // units requested for time interval
) const {

   Length = 0.0;
   R8 TmpResult{0.0}; // temporary for conversion
   I4 Err{0};

   // Now get values based on requested units
   // If requested units are same as the interval was defined,
   // just return length. Otherwise perform a conversion if needed or
   // error if not supported.

   switch (ReqUnits) {
   case TimeUnits::Years:
      if (Units == TimeUnits::Years) {
         Length = CalInterval;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) attempt to get interval "
                   "in years but defined in other units");
      }
      break;
   case TimeUnits::Months:
      if (Units == TimeUnits::Months) {
         Length = CalInterval;
      } else {
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) attempt to get interval "
                   "in months but defined in other units");
      }
      break;
   case TimeUnits::Days:
      switch (Units) {
      case TimeUnits::Days:
         Length = CalInterval;
         break;
      case TimeUnits::Seconds:
      case TimeUnits::Minutes:
      case TimeUnits::Hours:
         TmpResult = Interval.getSeconds();
         Length    = TmpResult / static_cast<R8>(SECONDS_PER_DAY);
         break;
      case TimeUnits::Months:
      case TimeUnits::Years:
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) conversion of year/month "
                   "units to days not supported");
         break;
      case TimeUnits::None:
      default:
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) interval units not "
                   "defined");
      } // end switch on interval units
      break;
   case TimeUnits::Seconds:
      switch (Units) {         // result depends on units interval defined
      case TimeUnits::Seconds: // non-calendar interval simply
      case TimeUnits::Minutes: // retrieves from TimeFrac
      case TimeUnits::Hours:
         Length = Interval.getSeconds();
         break;
      case TimeUnits::Days:
         TmpResult = CalInterval;
         Length    = TmpResult * static_cast<R8>(SECONDS_PER_DAY);
         break;
      case TimeUnits::Years:  // some calendar intervals cannot
      case TimeUnits::Months: // be converted
         Length = 0.0;
         Err    = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) conversion of year/month "
                   "intervals to seconds unsupported");
         break;
      case TimeUnits::None:
      default:
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) interval units not "
                   "defined");
      } // end switch on interval units
      break;
   case TimeUnits::Minutes:
      switch (Units) {         // result depends on units interval defined
      case TimeUnits::Seconds: // non-calendar interval simply
      case TimeUnits::Minutes: // retrieves from Base Time
      case TimeUnits::Hours:
         Length = Interval.getMinutes();
         break;
      case TimeUnits::Years:  // some calendar intervals cannot
      case TimeUnits::Months: // be converted
         Length = 0.0;
         Err    = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) conversion of year/month "
                   "intervals to minutes unsupported");
         break;
      case TimeUnits::Days:
         TmpResult = CalInterval;
         Length    = TmpResult * static_cast<R8>(SECONDS_PER_DAY) /
                  static_cast<R8>(SECONDS_PER_MINUTE);
         break;
      case TimeUnits::None:
      default:
         Err = 1;
         LOG_ERROR(
             "TimeMgr: TimeInterval::get (real) interval units not defined");
      } // end switch on interval units
      break;
   case TimeUnits::Hours:
      switch (Units) {         // result depends on units interval defined
      case TimeUnits::Seconds: // non-calendar interval simply
      case TimeUnits::Minutes: // retrieves from Base Time
      case TimeUnits::Hours:
         Length = Interval.getHours();
         break;
      case TimeUnits::Years:  // some calendar intervals cannot
      case TimeUnits::Months: // be converted
         Length = 0.0;
         Err    = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) conversion of year/month "
                   "intervals to seconds unsupported");
         break;
      case TimeUnits::Days:
         TmpResult = CalInterval;
         Length    = TmpResult * static_cast<R8>(SECONDS_PER_DAY) /
                  static_cast<R8>(SECONDS_PER_HOUR);
         break;
      case TimeUnits::None:
      default:
         Err = 1;
         LOG_ERROR("TimeMgr: TimeInterval::get (real) interval units not "
                   "defined");
      } // end switch on interval units
      break;
   case TimeUnits::None:
   default:
      Err = 1;
      LOG_ERROR("TimeMgr: TimeInterval::get (real) invalid time units for time "
                "interval");
   } // end switch on requested units

   // return error code
   return Err;

} // end TimeInterval::get (real)

// TimeInterval operators
//------------------------------------------------------------------------------
// TimeInterval::operator==
// This equivalence operator returns true if the time intervals are
// equal. If time intervals are set correctly, this simply means
// checking whether all components are equal.

bool TimeInterval::operator==(
    const TimeInterval &TimeInt) const { // [in] - interval to compare

   // check calendar components if a calendar interval
   // both units and length of interval must match
   if (this->IsCalendar) {
      if (!TimeInt.IsCalendar)
         return false;
      if (this->Units != TimeInt.Units)
         return false;
      if (this->CalInterval != TimeInt.CalInterval)
         return false;
   }
   // check non-calendar interval using base time equivalence
   if (!(this->Interval == TimeInt.Interval))
      return false;

   // if we made it here, they should be equivalent
   return true;

} // end TimeInterval::operator==

//------------------------------------------------------------------------------
// TimeInterval::operator!=
// This non-equivalence operator returns true if time intervals are
// not equal. If time intervals are set correctly, this simply means
// checking whether all components are not equal.

bool TimeInterval::operator!=(
    const TimeInterval &TimeInt) const { // [in] - interval to compare

   // use the equivalence operator above
   return !(*this == TimeInt);

} // end TimeInterval::operator!=

//------------------------------------------------------------------------------
// TimeInterval::operator<
// This operator returns true if this time interval is less than
// another. Time intervals must be of the same type. That is, they
// both must be calendar or non-calendar intervals and if calendar
// interval, they must define the interval in the same units
// (years, months, or days).

bool TimeInterval::operator<(
    const TimeInterval &TimeInt) const { // [in] - interval to compare

   // calendar intervals
   if (this->IsCalendar) {
      // cannot determine result if both are not calendar intervals
      if (!(TimeInt.IsCalendar)) {
         LOG_ERROR("TimeMgr: TimeInterval::operator< attempt to compare "
                   "incompatible time intervals");
         return false;
      }
      // calendar must also compare intervals with same units
      if (this->Units == TimeInt.Units) {
         return this->CalInterval < TimeInt.CalInterval;
      } else {
         LOG_ERROR("TimeMgr: TimeInterval::operator< attempt to compare "
                   "incompatible calendar intervals");
         return false;
      }
   } else {
      // non-calendar intervals use inherited operator
      // check for incompatibility
      if (TimeInt.IsCalendar) {
         LOG_ERROR("TimeMgr: TimeInterval::operator< attempt to compare "
                   "incompatible time intervals");
         return false;
      } else {
         // use base time comparison
         return this->Interval < TimeInt.Interval;
      }
   } // end else (non-calendar)

} // end TimeInterval::operator<

//------------------------------------------------------------------------------
// TimeInterval::operator>
// This operator returns true if this time interval is greater than
// another. Time intervals must be of the same type. That is, they
// both must be calendar or non-calendar intervals and if calendar
// interval, they must define the interval in the same units
// (years, months, or days).

bool TimeInterval::operator>(
    const TimeInterval &TimeInt) const { // [in] - interval to compare

   // calendar intervals
   if (this->IsCalendar) {
      // cannot determine result if both are not calendar intervals
      if (!(TimeInt.IsCalendar)) {
         LOG_ERROR("TimeMgr: TimeInterval::operator> attempt to compare "
                   "incompatible time intervals");
         return false;
      }
      // calendar must also compare intervals with same units
      if (this->Units == TimeInt.Units) {
         return this->CalInterval > TimeInt.CalInterval;
      } else {
         LOG_ERROR("TimeMgr: TimeInterval::operator> attempt to compare "
                   "incompatible calendar intervals");
         return false;
      }
   } else {
      // non-calendar intervals use inherited operator
      // check for incompatibility
      if (TimeInt.IsCalendar) {
         LOG_ERROR("TimeInterval::operator> attempt to compare incompatible "
                   "time intervals");
         return false;
      } else {
         // use TimeFrac comparison
         return this->Interval > TimeInt.Interval;
      }
   }

} // end TimeInterval::operator>

//------------------------------------------------------------------------------
// TimeInterval::operator<=
// This operator returns true if this time interval is less than or
// equal to another. Time intervals must be of the same type. That is,
// they both must be calendar or non-calendar intervals and if calendar
// interval, they must define the interval in the same units
// (years, months, or days).

bool TimeInterval::operator<=(
    const TimeInterval &TimeInt) const { // [in] - interval to compare

   // use existing existing operators
   return (*this < TimeInt) || (*this == TimeInt);

} // end TimeInterval::operator<=

//------------------------------------------------------------------------------
// TimeInterval::operator>=
// This operator returns true if this time interval is greater than
// or equal to another. Time intervals must be of the same type. That
// is, they both must be calendar or non-calendar intervals and if
// calendar interval, they must define the interval in the same units
// (years, months, or days).

bool TimeInterval::operator>=(
    const TimeInterval &TimeInt) const { // [in] - interval to compare

   // use existing existing operators
   return (*this > TimeInt) || (*this == TimeInt);

} // end TimeInterval::operator>=

//------------------------------------------------------------------------------
// TimeInterval::operator+
// This operator adds two time intervals and returns the result.

TimeInterval TimeInterval::operator+(
    const TimeInterval &TimeInt) const { // [in] - interval to add

   // start by copying result
   TimeInterval Result = *this;

   // can only really add compatible time intervals
   // add calendar intervals
   if (Result.IsCalendar) {
      if (!(TimeInt.IsCalendar)) {
         LOG_ERROR("TimeMgr: TimeInterval::operator+ attempt to add "
                   "incompatible time intervals");
         return Result;
      }
      // calendar can only add intervals with same units
      if (this->Units == TimeInt.Units) {
         Result.CalInterval += TimeInt.CalInterval;
      } else {
         LOG_ERROR("TimeMgr: TimeInterval::operator+ attempt to add "
                   "incompatible calendar intervals");
         return Result;
      }
   } else {
      // add non-calendar using TimeFrac operator
      // first check that both operands are non-calendar intervals
      if (TimeInt.IsCalendar) {
         LOG_ERROR("TimeMgr: TimeInterval::operator+ attempt to add "
                   "incompatible time intervals");
         return Result;
      }

      // use underlying TimeFrac operator to add the TimeFrac component
      Result.Interval += TimeInt.Interval;
   }

   // return the result
   return Result;

} // end TimeInterval::operator+

//------------------------------------------------------------------------------
// TimeInterval::operator-
// This operator subtracts two time intervals and returns the result.

TimeInterval TimeInterval::operator-(
    const TimeInterval &TimeInt) const { // [in] - interval to subtract

   // start by copying result
   TimeInterval Result = *this;

   // can only really subtract compatible time intervals
   // subtract calendar intervals
   if (Result.IsCalendar) {
      if (!(TimeInt.IsCalendar)) {
         LOG_ERROR("TimeMgr: TimeInterval::operator- attempt to subtract "
                   "incompatible time intervals");
         return Result;
      }
      // calendar can only subtract intervals with same units
      if (this->Units == TimeInt.Units) {
         Result.CalInterval -= TimeInt.CalInterval;
      } else {
         LOG_ERROR("TimeMgr: TimeInterval::operator- attempt to subtract "
                   "incompatible calendar intervals");
         return Result;
      }
   } else {
      // subtract non-calendar intervals using TimeFrac operator
      // check compatibility of operands
      if (TimeInt.IsCalendar) {
         LOG_ERROR("TimeInterval::operator- attempt to subtract incompatible "
                   "time intervals");
         return Result;
      }

      // use TimeFrac operator on the TimeFrac component
      Result.Interval -= TimeInt.Interval;
   }

   // return the result
   return Result;

} // end TimeInterval::operator-

//------------------------------------------------------------------------------
// TimeInterval::operator+=
// This operator increments a time interval by adding another interval and
// storing the result in place.

TimeInterval &TimeInterval::operator+=(const TimeInterval &TimeInt) {

   // reuse (+) operator defined above
   *this = *this + TimeInt;
   return *this;

} // end TimeInterval::operator+=

//------------------------------------------------------------------------------
// TimeInterval::operator-=
// This operator decrements a time interval by subtracting another interval and
// storing the result in place.

TimeInterval &TimeInterval::operator-=(const TimeInterval &TimeInt) {

   // reuse (-) operator defined above
   *this = *this - TimeInt;
   return *this;

} // end TimeInterval::operator-=

//------------------------------------------------------------------------------
// TimeInterval::operator*
// This operator increases the size of a time interval by an integer multiple.

TimeInterval
TimeInterval::operator*(const I4 Multiplier) const { // [in] integer multiplier

   // start with a copy of the input time interval
   TimeInterval Product = *this;

   // for all interval types, just multiply values by integer

   Product.Interval *= Multiplier;
   Product.CalInterval *= Multiplier;

   // return the result
   return Product;

} // end TimeInterval::operator* for integers

//------------------------------------------------------------------------------
// TimeInterval::operator*
// This operator increases the size of a time interval by a real
// multiple. For calendar intervals, the year, month or day is
// multiplied as a real number, then converted back to an integer
// so some rounding of the result occurs to maintain the integer interval.

TimeInterval
TimeInterval::operator*(const R8 Multiplier) const { // [in] real multiplier

   // start with a copy of the input time interval
   TimeInterval Product = *this;

   // for non-calendar types, just multiply values by real variable
   // using the base time operator
   Product.Interval *= Multiplier;

   // for calendar intervals, must convert (and round) back to integer
   R8 Tmp{0.0};
   Tmp                 = static_cast<R8>(Product.CalInterval) * Multiplier;
   Product.CalInterval = static_cast<I8>(Tmp);

   // return the result
   return Product;

} // end TimeInterval::operator* for reals

//------------------------------------------------------------------------------
// TimeInterval::operator*
// This operator increases the size of a time interval by an integer
// multiple. It is the commutative complement to the member operator
// in which the integer multiplier is the first argument.

TimeInterval
operator*(const I4 &Multiplier,     // [in] integer multiplier
          const TimeInterval &TI) { // [in] TimeInterval multiplicand

   // just call the member operator in the right order
   return TI * Multiplier;

} // end TimeInterval::operator* for integers commutative version

//------------------------------------------------------------------------------
// TimeInterval::operator*
// This operator increases the size of a time interval by an real
// multiple. It is the commutative complement to the member operator
// in which the real multiplier is the first argument. As in that case,
// calendar intervals are rounded to an integer using a static cast.

TimeInterval
operator*(const R8 &Multiplier,     // [in] real (dble) multiplier
          const TimeInterval &TI) { // [in] TimeInterval multiplicand

   // just call the member operator in the right order
   return TI * Multiplier;

} // end TimeInterval::operator* for reals commutative version

//------------------------------------------------------------------------------
// TimeInterval::operator*=
// This operator increases the size of a time interval by an integer multiple.

TimeInterval &
TimeInterval::operator*=(const I4 Multiplier) { // [in] integer multiplier

   // call the member multiply operator and return
   return *this = *this * Multiplier;

} // end TimeInterval::operator*= for integer multipliers

//------------------------------------------------------------------------------
// TimeInterval::operator*=
// This operator increases the size of a time interval by an real
// multiple and stores the result in place. For calendar-dependent intervals,
// the multiplication is carried out as a real operation, but the result is
// converted back to integer form, so some rounding occurs.

TimeInterval &
TimeInterval::operator*=(const R8 Multiplier) { // [in] real multiplier

   // call the member multiply operator and return
   *this = *this * Multiplier;
   return *this;

} // end TimeInterval::operator*= for real multipliers

//------------------------------------------------------------------------------
// TimeInterval::operator/
// This routine creates a new time interval by subdividing a current
// interval by an integral number. For calendar-independent intervals
// this simply divides the fractional integer interval. For calendar-
// dependent intervals, this will divide the year, month or day
// interval using integer division. If the year, month or day does
// not divide evenly, a warning is generated but the result of the
// integer division is returned.

TimeInterval
TimeInterval::operator/(const I4 Divisor) const { //! [in] integer divisor

   // check for divide by zero and return a zero default interval
   TimeInterval Quotient;
   if (Divisor == 0) {
      LOG_ERROR("TimeMgr: TimeInterval::operator/ (int) attempt to divide by "
                "zero");
      return Quotient;
   }

   // copy the time interval to the result
   Quotient = *this;

   // for calendar-dependent intervals, divide the
   // year, month or day intervals
   if (Quotient.IsCalendar) {
      // generate warning if divisor does not divide evenly
      if (Quotient.CalInterval != 0 && Quotient.CalInterval % Divisor != 0) {
         LOG_WARN("TimeMgr: TimeInterval::operator/ (int) interval does not "
                  "divide evenly");
      }
      // return the result of the integer division
      Quotient.CalInterval /= Divisor;

   } else {
      // for calendar-independent intervals, use the TimeFrac
      // operator to divide fractional seconds
      Quotient.Interval /= Divisor;
   }

   // return the result
   return Quotient;

} // end TimeInterval::operator/ for int divisor

//------------------------------------------------------------------------------
// TimeInterval::operator/=
// This routine subdivides a current interval by an integral number
// and returns the result in place. For calendar-independent intervals
// this simply divides the fractional integer interval. For calendar-
// dependent intervals, this will divide the year, month or day
// interval using integer division. If the year, month or day does
// not divide evenly, a warning is generated but the result of the
// integer division is returned.

TimeInterval &
TimeInterval::operator/=(const I4 Divisor) { // in - integer divisor

   // simply call the related integer division function and
   // return in place
   return *this = *this / Divisor;

} // end TimeInterval::operator/= for int division

// TimeInterval methods
//------------------------------------------------------------------------------
// TimeInterval::absValue
// This function returns the absolute value of the given time interval.
// For calendar-independent intervals, this is just the positive value
// of the fractional second representation.  For calendar-dependent
// intervals it is the abs value of the year, month or day interval.

TimeInterval TimeInterval::absValue(
    const TimeInterval &TimeInt) { // [in] interval for abs value

   // initialize result to a copy of the current interval
   TimeInterval AbsValue = TimeInt;

   // take abs value of calendar-independent components
   I8 W{0};
   I8 N{0};
   I8 D{0};
   AbsValue.Interval.get(W, N, D);
   if (W < 0)
      AbsValue.Interval.setWhole(-1 * W);
   if (N < 0)
      AbsValue.Interval.setNumer(-1 * N);
   if (D < 0)
      AbsValue.Interval.setDenom(-1 * D);

   // take abs value of calendar-dependent components
   // if this is a calendar-dependent interval
   if (AbsValue.IsCalendar) {
      if (AbsValue.CalInterval < 0)
         AbsValue.CalInterval = -1 * AbsValue.CalInterval;
   }

   return AbsValue;

} // end TimeInterval::absValue

//------------------------------------------------------------------------------
// TimeInterval::negAbsValue
// This function returns the negative absolute value of the given time
// interval. For calendar-independent intervals, this is just the
// negative absolute value of the fractional second representation.
// For calendar-dependent intervals it is the negative abs value of
// the year, month or day interval.

TimeInterval TimeInterval::negAbsValue(
    const TimeInterval &TimeInt) { // [in] interval for neg abs val

   // initialize result to a copy of the current interval
   TimeInterval NegAbsValue = TimeInt;

   // take negative abs value of calendar-independent components
   I8 W{0};
   I8 N{0};
   I8 D{0};
   NegAbsValue.Interval.get(W, N, D);
   if (W > 0)
      NegAbsValue.Interval.setWhole(-1 * W);
   if (N > 0)
      NegAbsValue.Interval.setNumer(-1 * N);
   if (D < 0)
      NegAbsValue.Interval.setDenom(-1 * D); // must ensure denom positive

   // take negative abs value of calendar-dependent components
   // if this is a calendar-dependent interval
   if (NegAbsValue.IsCalendar) {
      if (NegAbsValue.CalInterval > 0)
         NegAbsValue.CalInterval = -1 * NegAbsValue.CalInterval;
   }

   return NegAbsValue;

} // end TimeInterval::negAbsValue

//------------------------------------------------------------------------------
// TimeInterval::isPositive
// Function that returns true if the time interval is positive

bool TimeInterval::isPositive(void) {

   // initialize result to false
   bool Result = false;

   // Check components of a non-calendar interval
   I8 W{0};
   I8 N{0};
   I8 D{0};
   // make sure fraction is in proper form, particularly that the
   // whole and fraction have same sign and numerator and denomintator
   // have the same sign.
   I4 Err = Interval.simplify();

   // now retrieve fractional seconds and check components
   Err = Interval.get(W, N, D);
   if (W >= 0) { // whole part of seconds non-negative
      // only need check the numerator with simplified form
      if (N > 0)
         Result = true;
   }

   // Check calendar interval quantities
   if (IsCalendar) {
      if (CalInterval > 0)
         Result = true;
   }

   return Result;

} // end TimeInterval::isPositive

//------------------------------------------------------------------------------
// TimeInstant definitions
//------------------------------------------------------------------------------

// TimeInstant accessors
// Define first so constructors can use set methods
//------------------------------------------------------------------------------
// TimeInstant::set (date and time, real seconds)
// Sets a time instant from a date and time, where the seconds is supplied as
// a real number. The calendar must already be set in order to correctly use
// this function.

I4 TimeInstant::set(const I8 Year,   //< [in] year
                    const I8 Month,  //< [in] month
                    const I8 Day,    //< [in] day
                    const I8 Hour,   //< [in] hour
                    const I8 Minute, //< [in] minute
                    const R8 RSecond //< [in] second (real)
) {

   // Initialize error code
   I4 Err{0};

   // Convert real seconds to fractional representation
   TimeFrac TmpSeconds(RSecond);
   I8 Whole{0};
   I8 Numer{0};
   I8 Denom{0};
   Err = TmpSeconds.get(Whole, Numer, Denom);

   // Call Calendar function to determine elapsed time since reference time
   ElapsedTime = Calendar::getElapsedTime(Year, Month, Day, Hour, Minute, Whole,
                                          Numer, Denom);

   return Err;

} // end TimeInstant::set (date and time, real seconds)

//------------------------------------------------------------------------------
// TimeInstant::set (date and time, fractional integer seconds)
// Sets a time instant from a date and time, where the seconds is supplied as
// a fractional integer. The calendar must already be set in order to correctly
// use this function.

int TimeInstant::set(const I8 Year,   //< [in] year
                     const I8 Month,  //< [in] month
                     const I8 Day,    //< [in] day
                     const I8 Hour,   //< [in] hour
                     const I8 Minute, //< [in] minute
                     const I8 Whole,  //< [in] second (whole integer)
                     const I8 Numer,  //< [in] second (fraction numerator)
                     const I8 Denom   //< [in] second (fraction denominator)
) {

   // Initialize error code
   I4 Err{0};

   // Call Calendar function to determine elapsed time since reference time
   ElapsedTime = Calendar::getElapsedTime(Year, Month, Day, Hour, Minute, Whole,
                                          Numer, Denom);

   return Err;

} // end TimeInstant::set (date and time, fractional integer seconds)

//------------------------------------------------------------------------------
// TimeInstant::get (date and time, real seconds)
// Retrieve the time instant in terms of a date and time with seconds returned
// as a real (double) value.

I4 TimeInstant::get(I8 &Year,   //< [out] year   of this time instant
                    I8 &Month,  //< [out] month  of this time instant
                    I8 &Day,    //< [out] day    of this time instant
                    I8 &Hour,   //< [out] hour   of this time instant
                    I8 &Minute, //< [out] minute of this time instant
                    R8 &Second  //< [out] second of this time instant
) const {

   // Initialize error code
   I4 Err{0};

   Second = 0.0;

   // Call the calendar function to convert the time to a date and time
   // appropriate for the calendar
   I8 Whole{0};
   I8 Numer{0};
   I8 Denom{0};
   Err = Calendar::getDateTime(ElapsedTime, Year, Month, Day, Hour, Minute,
                               Whole, Numer, Denom);
   if (Err != 0) {
      LOG_ERROR("TimeMgr: TimeInstant::get(date and time, real seconds) "
                "error converting base time to date/time");
   } else {
      // Convert retrieved seconds in TimeFrac and then
      // to real seconds for return
      TimeFrac TmpSeconds(Whole, Numer, Denom);
      Second = TmpSeconds.getSeconds();
   }

   // return the error code
   return Err;

} // end TimeInstant::get (date and time, real seconds)

//------------------------------------------------------------------------------
// TimeInstant::get (date and time, fractional integer seconds)
// Retrieve the time instant in terms of a date and time with seconds returned
// as a fractional integer.

I4 TimeInstant::get(I8 &Year,   //< [out] year   of this time instant
                    I8 &Month,  //< [out] month  of this time instant
                    I8 &Day,    //< [out] day    of this time instant
                    I8 &Hour,   //< [out] hour   of this time instant
                    I8 &Minute, //< [out] minute of this time instant
                    I8 &Whole,  //< [out] whole seconds of this time
                    I8 &Numer,  //< [out] frac second numerator
                    I8 &Denom   //< [out] frac second denominator
) const {

   // Initialize error code
   I4 Err{0};

   // Call the calendar function to convert the time to a date and time
   // appropriate for that calendar
   Err = Calendar::getDateTime(ElapsedTime, Year, Month, Day, Hour, Minute,
                               Whole, Numer, Denom);
   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInstant::get(date and time, frac int seconds) "
                "error converting base time to date/time");

   // return the error code
   return Err;

} // end TimeInstant::get (date and time, fractional integer seconds)

// TimeInstant constructors/destructors
//------------------------------------------------------------------------------
// TimeInstant::TimeInstant - default constructor
// This is a default constructor that creates an empty time instant with zero
// values and a null calendar.

TimeInstant::TimeInstant(void) : ElapsedTime(0, 0, 1) {
   // Everything taken care of in initialization list above
} // end TimeInstant::TimeInstant

//------------------------------------------------------------------------------
// TimeInstant::TimeInstant
// Constructs a time instant from a date, time, calendar, where the seconds
// part of time is supplied as a real number.

TimeInstant::TimeInstant(const I8 Year,   //< [in] year
                         const I8 Month,  //< [in] month
                         const I8 Day,    //< [in] day
                         const I8 Hour,   //< [in] hour
                         const I8 Minute, //< [in] minute
                         const R8 RSecond //< [in] second (real)
) {

   // Reuse the set function to define the time component
   I4 Err = this->set(Year, Month, Day, Hour, Minute, RSecond);
   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInstant::TimeInstant "
                "(date and time, real seconds) "
                "error using set function to construct TimeInstant.");

} // end TimeInstant::TimeInstant constructor from yy/mm/dd-hh:mm:ss(real)

//------------------------------------------------------------------------------
// TimeInstant::TimeInstant
// Constructs a time instant from a calendar, date and time, where the seconds
// is supplied as an integer fraction (whole + numrator/denominator)

TimeInstant::TimeInstant(const I8 Year,   //< [in] year
                         const I8 Month,  //< [in] month
                         const I8 Day,    //< [in] day
                         const I8 Hour,   //< [in] hour
                         const I8 Minute, //< [in] minute
                         const I8 Whole,  //< [in] second (whole integer)
                         const I8 Numer,  //< [in] second (fraction numerator)
                         const I8 Denom   //< [in] second (fraction denominator)
) {

   // Reuse the set function to define the time component
   I4 Err = this->set(Year, Month, Day, Hour, Minute, Whole, Numer, Denom);
   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInstant::TimeInstant "
                "(date and time, frac int seconds) "
                "error using set function to construct TimeInstant.");

} // end TimeInstant::TimeInstant constructor from yy/mm/dd-hh:mm:ss(real ss)

//------------------------------------------------------------------------------
// TimeInstant::TimeInstant
// Constructs a time instant from a string of the general form
// YYYYY-MM-DD_HH:MM:SS.SSS where the width of the YY and SS fields can
// be of arbitrary width and the separators can be any non-numeric character

TimeInstant::TimeInstant(std::string &TimeString //< [in] string with date/time
) {

   // Extract variables from string
   I8 Year    = 0;
   I8 Month   = 0;
   I8 Day     = 0;
   I8 Hour    = 0;
   I8 Minute  = 0;
   R8 RSecond = 0.;

   std::istringstream ss(TimeString);
   char discard;
   ss >> Year >> discard >> Month >> discard >> Day >> discard >> Hour >>
       discard >> Minute >> discard >> RSecond;

   // Use the set function to define the Time Instant
   I4 Err = this->set(Year, Month, Day, Hour, Minute, RSecond);
   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInstant constructor from string "
                "error using set function to construct TimeInstant.");

} // end TimeInstant::TimeInstant constructor from string

//------------------------------------------------------------------------------
// TimeInstant::~TimeInstant - destructor for time instant
// Destructor for a time instant. No allocated space so does nothing.

TimeInstant::~TimeInstant(void) { // Nothing to be done

} // end TimeInstant destructor

// TimeInstant operators
//------------------------------------------------------------------------------
// TimeInstant::operator==
// This equivalence operator returns true if two time instants are
// equivalent. If time instants are set correctly, this simply means
// checking whether all components are equivalent using the basetime and
// calendar equivalence operators.

bool TimeInstant::operator==(
    const TimeInstant &Instant) const { // [in] - instant to compare

   // check equivalence of both calendar and basetime components
   // Note that the calendars do not need to be the same calendar,
   //   simply equivalent calendars

   if (this->ElapsedTime == Instant.ElapsedTime)
      return true;
   else
      return false;

} // end TimeInstant::operator==

//------------------------------------------------------------------------------
// TimeInstant::operator!=
// This non-equivalence operator returns true if two time instants are
// not equal. If time instants are set correctly, this simply means
// checking whether all components are not equivalent.

bool TimeInstant::operator!=(
    const TimeInstant &Instant) const { // [in] - instant to compare

   // use the equivalence operator above
   return (!(*this == Instant));

} // end TimeInstant::operator!=

//------------------------------------------------------------------------------
// TimeInstant::operator<
// This operator returns true if this time instant is less than
// another. Time instant must be defined on the same calendar.

bool TimeInstant::operator<(
    const TimeInstant &Instant) const { // [in] - instant to compare

   return (this->ElapsedTime < Instant.ElapsedTime);

} // end TimeInstant::operator<

//------------------------------------------------------------------------------
// TimeInstant::operator>
// This operator returns true if this time instant is greater than
// another. Time instant must be defined on the same calendar.

bool TimeInstant::operator>(
    const TimeInstant &Instant) const { // [in] - instant to compare

   return (this->ElapsedTime > Instant.ElapsedTime);

} // end TimeInstant::operator>

//------------------------------------------------------------------------------
// TimeInstant::operator<=
// This operator returns true if this time instant is less than or equal to
// another. Time instant must be defined on the same calendar.

bool TimeInstant::operator<=(
    const TimeInstant &Instant) const { // [in] - instant to compare

   return (this->ElapsedTime <= Instant.ElapsedTime);

} // end TimeInstant::operator<=

//------------------------------------------------------------------------------
// TimeInstant::operator>=
// This operator returns true if this time instant is greater than or equal to
// another. Time instant must be defined on the same calendar.

bool TimeInstant::operator>=(
    const TimeInstant &Instant) const { // [in] - instant to compare

   return (this->ElapsedTime >= Instant.ElapsedTime);

} // end TimeInstant::operator>=

//------------------------------------------------------------------------------
// TimeInstant::operator+
// This operator is used to advance time by adding a time interval.

TimeInstant TimeInstant::operator+(
    const TimeInterval &Interval) const { // [in] - interval to add

   I4 Err{0};

   // start by copying result
   TimeInstant Result = *this;

   // for non-calendar intervals, we can simply add the ElapsedTime components
   if (!(Interval.IsCalendar))
      Result.ElapsedTime = ElapsedTime + Interval.Interval;

   else {
      // for calendar intervals, interval depends on current date
      // convert current time instant to date/time form
      I8 Year{0};
      I8 Month{0};
      I8 Day{0};
      I8 Hour{0};
      I8 Minute{0};
      I8 Whole{0};
      I8 Numer{0};
      I8 Denom{1};
      Err = Calendar::getDateTime(ElapsedTime, Year, Month, Day, Hour, Minute,
                                  Whole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: TimeInstant::operator+ error converting instant "
                   "to date/time");
      // use calendar function to increment date based on calendar kind
      Err = Calendar::incrementDate(Interval.CalInterval, Interval.Units, Year,
                                    Month, Day);
      if (Err != 0)
         LOG_ERROR("TimeMgr: TimeInstant::operator+ error advancing calendar "
                   "date/time");

      // convert new date/time back to time instant
      Err = Result.set(Year, Month, Day, Hour, Minute, Whole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: TimeInstant::operator+ error converting date/time "
                   "to instant");

   } // end calendar else clause

   // return the result
   return Result;

} // end TimeInstant::operator+

//------------------------------------------------------------------------------
// TimeInstant::operator-
// This operator is used to reverse time by subtracting a time interval.

TimeInstant TimeInstant::operator-(
    const TimeInterval &Interval) const { // [in] - interval to subtract

   I4 Err{0};

   // start by copying result
   TimeInstant Result = *this;

   // for non-calendar intervals, we can simply subtract the ElapsedTime
   // components
   if (!(Interval.IsCalendar))
      Result.ElapsedTime = ElapsedTime - Interval.Interval;

   else {
      // for calendar intervals, interval depends on current date
      // convert current time instant to date/time form
      I8 Year{0};
      I8 Month{0};
      I8 Day{0};
      I8 Hour{0};
      I8 Minute{0};
      I8 Whole{0};
      I8 Numer{0};
      I8 Denom{1};
      Err = Calendar::getDateTime(ElapsedTime, Year, Month, Day, Hour, Minute,
                                  Whole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: TimeInstant::operator- error converting instant "
                   "to date/time");
      // use calendar function to decrement date based on calendar kind
      I8 TmpInterval = -(Interval.CalInterval);
      Err = Calendar::incrementDate(TmpInterval, Interval.Units, Year, Month,
                                    Day);
      if (Err != 0)
         LOG_ERROR("TimeMgr: TimeInstant::operator- error reversing calendar "
                   "date/time");

      // convert new date/time back to time instant
      Err = Result.set(Year, Month, Day, Hour, Minute, Whole, Numer, Denom);
      if (Err != 0)
         LOG_ERROR("TimeMgr: TimeInstant::operator- cannot compare time "
                   "instants on different calendars");

   } // end calendar else clause

   // return the result
   return Result;

} // end TimeInstant::operator-

//------------------------------------------------------------------------------
// TimeInstant::operator- (time interval)
// This operator creates a time interval as a difference of two time instants

TimeInterval TimeInstant::operator-(const TimeInstant &Instant) const {

   I4 Err{0};
   I8 Whole{0};
   I8 Numer{0};
   I8 Denom{1};
   TimeInterval Result(Whole, Numer, Denom); // zero interval as default

   // subtract the ElapsedTime components of time instants
   TimeFrac Timediff = ElapsedTime - Instant.ElapsedTime;

   // extract ElapsedTime components of time difference and use them to
   // set time interval result
   Err = Timediff.get(Whole, Numer, Denom);
   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInstant::operator- (interval) error retrieving "
                "time difference components");

   Err = Result.set(Whole, Numer, Denom);
   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInstant::operator- (interval) error setting "
                "TimeFrac components of interval");

   // finally return time interval result
   return Result;

} // end TimeInstant::operator- (create time interval)

//------------------------------------------------------------------------------
// TimeInstant::operator+=
// This operator increments a time instant by adding a time interval and
// storing the result in place.

TimeInstant &TimeInstant::operator+=(const TimeInterval &Timeint) {

   // reuse (+) operator defined above
   *this = *this + Timeint;
   return *this;

} // end TimeInstant::operator+=

//------------------------------------------------------------------------------
// TimeInstant::operator-=
// This operator decrements a time instant by subtracting a time interval and
// storing the result in place.

TimeInstant &TimeInstant::operator-=(const TimeInterval &Timeint) {

   // reuse (-) operator defined above
   *this = *this - Timeint;
   return *this;

} // end TimeInstant::operator-=

// utility methods
//------------------------------------------------------------------------------
// TimeInstant::getString
// Get time instant as a string in the format
// 'YYYYYY-MM-DD_HH:MM:SS.SSSSSS' where number of digits in year, and number
// of digits after the decimal are inputs. The string separating date and time
// is also input.

std::string TimeInstant::getString(
    const I4 YearWidth,         // [in] number of digits in year
    const I4 SecondWidth,       // [in] num of decimal digits in seconds
    const std::string Separator // [in] string to separate date/time
) const {

   // First retrieve time/day
   I8 Year{0};
   I8 Month{0};
   I8 Day{0};
   I8 Hour{0};
   I8 Minute{0};
   R8 Second{0.0};
   I4 Err = this->get(Year, Month, Day, Hour, Minute, Second);
   if (Err != 0)
      LOG_ERROR("TimeMgr: TimeInstant::getString error retrieving time/date");

   // allocate character string for date/time
   I4 TimeLength = YearWidth + 3 + 3 + 1 + 3 + 3 + 3 + SecondWidth + 1;
   char *TimeStr = new char[TimeLength];

   // Use formatted print to create string

   // First define format
   std::ostringstream FmtString;
   FmtString.str(std::string());          // clear format stream
   FmtString << "%0" << YearWidth << "d"; // 0 padded year of defined width
   FmtString << "-";
   FmtString << "%02d"; // 2 digit month with leading 0
   FmtString << "-";
   FmtString << "%02d";    // 2 digit day with leading 0
   FmtString << Separator; // user-supplied separator between date/time
   FmtString << "%02d";    // 2 digit hour with leading 0
   FmtString << ":";
   FmtString << "%02d"; // 2 digit minute with leading 0
   FmtString << ":";
   if (SecondWidth == 0) {
      FmtString << "%02.0f"; // seconds with no decimal digits
   } else {
      // seconds with SecondWidth number of decimal digits
      FmtString << "%0" << SecondWidth + 3 << "." << SecondWidth << "f";
   }

   // convert format string to an actual string for formatted write
   const std::string Tmp = FmtString.str();
   const char *FormatStr = Tmp.c_str();

   // now use sprintf to print to C string
   Err = sprintf(TimeStr, FormatStr, Year, Month, Day, Hour, Minute, Second);
   if (Err < 0)
      LOG_ERROR("TimeMgr: TimeInstant::getString error in sprintf");

   // now convert C string to std::string for result
   std::string Result(TimeStr);

   // deallocate temp strings
   delete[] TimeStr;

   // All done
   return Result;

} // end TimeInstant:: getString

//------------------------------------------------------------------------------
// Alarm definitions
//------------------------------------------------------------------------------

// Alarm constructors/destructors
//------------------------------------------------------------------------------
// Alarm::Alarm - default constructor
// The default constructor creates an alarm with ring time equal to zero
Alarm::Alarm(void) {

   Name    = "";
   Ringing = false;
   Stopped = false;

   Periodic = false;

   TimeInstant AlarmTime;
   RingTime = AlarmTime;

} // end Alarm::Alarm - default

// Alarm::Alarm - construct single instance alarm
// Construct a single alarm using the input ring time.

Alarm::Alarm(
    const std::string InName,   //< [in] Name of alarm
    const TimeInstant AlarmTime //< [in] Time at/after which alarm rings
) {

   Name    = InName; // assign input name to alarm
   Ringing = false;  // initialize ringer to false
   Stopped = false;  // alarm is initially on

   // This is not an interval alarm, so set interval properties appropriately
   Periodic = false; // this is not an interval alarm
   // ringInterval - will be set using default interval constructor
   // ringTimePrev - will be set using default instant constructor

   // Now set ringing time to desired input value
   RingTime = AlarmTime;

} // end Alarm::Alarm (single instance alarm)

//------------------------------------------------------------------------------
// Alarm::Alarm - Construct a periodic/interval alarm
// Construct a periodic/interval alarm based on an input time interval and a
// start time for the interval period.

Alarm::Alarm(
    const std::string InName,         //< [in] Name of alarm
    const TimeInterval AlarmInterval, //< [in] interval at which alarm rings
    const TimeInstant IntervalStart   //< [in] start time of first interval
) {

   Name    = InName; // set the alarm name to the input value
   Ringing = false;  // initialize alarm ringer to false
   Stopped = false;  // alarm is initially on

   // this is an interval alarm so set flags and times appropriately
   Periodic     = true;
   RingInterval = AlarmInterval;

   // Set the previous ring time to the start time
   RingTimePrev = IntervalStart;
   // Next alarm time is initially set to be start + interval
   // If the start time is set in the distant past (multiple intervals from
   // present) later calls to the updateStatus method will move the previous
   // ring time forward to be closer to the actual time
   RingTime = IntervalStart + AlarmInterval;

} // end Alarm::Alarm constructor for periodic/interval alarms

//------------------------------------------------------------------------------
// Alarm::~Alarm - destructor for alarm
// Destructor for an alarm. No allocated space so does nothing.

Alarm::~Alarm(void) { // Nothing to be done

} // end Alarm destructor

// Alarm methods
//------------------------------------------------------------------------------
// Alarm::isRinging - Check whether an alarm is ringing
// This function checks to see if an alarm is ringing by simply returning
// the ringing flag that is true if alarm is ringing, false otherwise.

bool Alarm::isRinging(void) { return Ringing; }

// end Alarm::isRinging

//------------------------------------------------------------------------------
// Alarm::updateStatus - Changes the alarm status based on current time
// Checks whether the alarm should ring based on the current (or supplied)
// time instant. If the instant is equal to, or later than the ring time,
// the alarm begins ringing until stopped.  Returns an error code that is
// generally true unless you try to check a stopped alarm.

I4 Alarm::updateStatus(const TimeInstant CurrentTime // [in] current time
) {

   I4 Err{0};
   if (!Stopped) { // if stopped, then status never changes
      if (CurrentTime >= RingTime)
         Ringing = true;
   }

   return Err;

} // end Alarm::updateStatus

//------------------------------------------------------------------------------
// Alarm::reset - Stops a ringing alarm and resets to a new alarm time
// Stops a ringing alarm and sets next a new ring time. If the alarm
// is a periodic/interval alarm, the next ring time is set to be the
// next interval boundary after the input time.  If the alarm is a
// single instance, the input time is used as the next alarm time.

I4 Alarm::reset(const TimeInstant InTime // [in] new alarm time
) {

   // initialize error code and stop the ringing alarm
   I4 Err{0};
   Ringing = false;

   // if this is an interval alarm, find the next alarm time after
   // the input time
   if (Periodic) {

      // first check that the input time is valid
      if (InTime < RingTime) {
         Err = 1;
         LOG_ERROR("TimeMgr: Alarm::reset input time less than the current "
                   "ring time");
      } else {
         // now move forward in time until the next interval is greater
         // than the input time
         while (RingTime <= InTime) {
            RingTimePrev = RingTime;
            RingTime += RingInterval;
         }
      }
   } else {
      // if this is a single instance, just set the new ring time to the
      // input value
      RingTime = InTime;
   }

   return Err;

} // end Alarm::reset

//------------------------------------------------------------------------------
// Alarm::stop - Stops a ringing alarm
// This function permanently stops a ringing alarm. Use the reset function
// if the alarm needs to be set for a later time.

I4 Alarm::stop(void) {

   I4 Err{0};
   Stopped = true;
   Ringing = false;
   return Err;

} // end Alarm::stop

//------------------------------------------------------------------------------
// Alarm::rename - Renames an existing timer

I4 Alarm::rename(const std::string NewName ///< [in] new name for alarm
) {

   I4 Err{0};
   Name = NewName;
   return Err;

} // end Alarm:: rename

//------------------------------------------------------------------------------
// Alarm::getName - retrieves name of an alarm
// Returns the name of an alarm.

std::string Alarm::getName() const { return Name; }

//------------------------------------------------------------------------------
// Alarm::getInterval - retrieves the time interval for periodic alarms

const TimeInterval *Alarm::getInterval() const { return &RingInterval; }

//------------------------------------------------------------------------------
// Alarm::getRingTimePrev - get last time the alarm rang

const TimeInstant *Alarm::getRingTimePrev(void) const { return &RingTimePrev; }

//------------------------------------------------------------------------------
// Clock definitions
//------------------------------------------------------------------------------

// Clock constructors/destructors
//------------------------------------------------------------------------------
// Clock::Clock - constructs a clock for tracking time within a model
// Construct a clock using a start time and time step.

Clock::Clock(const TimeInstant InStartTime, //< [in] Start time for clock
             const TimeInterval InTimeStep  //< [in] Time step to advance clock
) {

   // Set start time and time step to input values
   StartTime = InStartTime;
   TimeStep  = InTimeStep;

   // The default for current time on creation is the start time
   CurrTime = StartTime;

   // Compute prev, next time based on time step
   PrevTime = CurrTime - TimeStep;
   NextTime = CurrTime + TimeStep;

   // Initialize number of alarms
   NumAlarms = 0;

   // create an initial vector with null pointers for later attached alarms
   Alarms.resize(MAX_ALARMS);
   for (I4 I = 0; I < MAX_ALARMS; ++I) {
      Alarms[I] = nullptr;
   }

} // end Clock::Clock

//------------------------------------------------------------------------------
// Clock::~Clock - destructor for clock
// Destructor for a clock. No allocated space and std::vector does it's own
// cleanup, so nothing to do here.

Clock::~Clock(void) { // Nothing to be done

} // end Clock destructor

// Clock accessor methods
//------------------------------------------------------------------------------
// Clock::setCurrentTime - Sets the current time
// Sets the current time to an input value. Also must reset previous time and
// next time to be consistent. Check that new time does not precede start time.

I4 Clock::setCurrentTime(
    const TimeInstant InCurrTime // new value for current time
) {

   I4 Err{0};

   // Check that the new value does not precede start time
   if (InCurrTime < StartTime) {
      LOG_ERROR("TimeMgr: Clock::setCurrentTime value for current time "
                "precedes start time");

   } else {
      // otherwise reset currTime to input value and
      // recompute prev time, next time
      CurrTime = InCurrTime;
      PrevTime = CurrTime - TimeStep;
      NextTime = CurrTime + TimeStep;
   }

   // Update status of all attached alarms based on new time
   I4 Err1{0};
   for (I4 N = 0; N < NumAlarms; ++N) {
      Err1 = Alarms[N]->updateStatus(CurrTime);
      if (Err1 != 0) {
         ++Err;
         LOG_ERROR("Clock::setCurrentTime error updating alarm # {}", N);
         break;
      }
   }
   return Err;

} // end Clock::setCurrentTime

//------------------------------------------------------------------------------
// Clock::changeTimeStep - Changes the time step for this clock
// Simply changes the time step. Does not check to see if the new time step
// is reasonable. With the new time step, must also update the next time.

I4 Clock::changeTimeStep(
    const TimeInterval InTimeStep // [in] new time step to use
) {

   I4 Err{0};

   // Assign the new time step to this clock
   TimeStep = InTimeStep;

   // Update the next time based on new time step
   NextTime = CurrTime + TimeStep;

   return Err;

} // end Clock::changeTimeStep

//------------------------------------------------------------------------------
// Clock::getCurrentTime - Retrieves current time of this clock
// Returns the value of the current time in this clock.

TimeInstant Clock::getCurrentTime(void) const { return CurrTime; }

//------------------------------------------------------------------------------
// Clock::getPreviousTime - Retrieves time at previous time step
// Returns the time associated with the previous time step.

TimeInstant Clock::getPreviousTime(void) const { return PrevTime; }

//------------------------------------------------------------------------------
// Clock::getNextTime - Retrieves time at next time step
// Returns the time associated with the next time step.

TimeInstant Clock::getNextTime(void) const { return NextTime; }

//------------------------------------------------------------------------------
// Clock::getStartTime - Retrieves start time of this clock
// Returns the start time for this clock

TimeInstant Clock::getStartTime(void) const { return StartTime; }

//------------------------------------------------------------------------------
// Clock::getTimeStep - Retrieves the time step used by this clock
// Returns the time step for this clock

TimeInterval Clock::getTimeStep(void) const { return TimeStep; }

//------------------------------------------------------------------------------
// Clock::attachAlarm - Attaches an alarm to this clock
// Attaches an alarm to this clock. The clock simply stores a pointer to this
// alarm, so user must be careful not to destroy any alarms before the clock.

I4 Clock::attachAlarm(Alarm *InAlarm // [in] pointer to alarm to attach
) {

   I4 Err{0};

   // First update the alarm count and resize array if needed
   NumAlarms += 1;

   I4 AlarmsSize = Alarms.size();
   if (NumAlarms > AlarmsSize) {
      AlarmsSize += MAX_ALARMS;
      Alarms.resize(AlarmsSize);
   }

   // Now add the pointer to the next available slot in the array
   Alarms[NumAlarms - 1] = InAlarm;

   return Err;

} // end Clock::attachAlarm

// Clock methods
//------------------------------------------------------------------------------
// Clock::advance - Advances a clock one timestep and updates alarms
// This method advances a clock one interval and updates the status of all
// attached alarms.

I4 Clock::advance(void) {

   I4 Err{0};

   // Update previous time and current time from previously computed times
   PrevTime = CurrTime;
   CurrTime = NextTime;

   // Advance next time
   NextTime += TimeStep;

   // Update status of all attached alarms
   I4 Err1{0};
   for (I4 N = 0; N < NumAlarms; ++N) {
      Err1 = Alarms[N]->updateStatus(CurrTime);
      if (Err1 != 0) {
         ++Err;
         LOG_ERROR("TimeMgr: Clock::advance error updating alarm # {}", N);
         break;
      }
   }

   return Err;

} // end Clock::advance

} // namespace OMEGA
//===-----------------------------------------------------------------------===/
