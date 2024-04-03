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

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>

// max/min macros if they don't already exist
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

namespace OMEGA {

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

} // namespace OMEGA
//===-----------------------------------------------------------------------===/
