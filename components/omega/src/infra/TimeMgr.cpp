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

//------------------------------------------------------------------------------
// Calendar definitions
//------------------------------------------------------------------------------

// initialize static calendar instance counter
I4 Calendar::NumCalendars = 0;

// Calendar accessors
//------------------------------------------------------------------------------
// Calendar::rename - Resets a calendar's name
// Resets the name of a calendar. This is the only variable that
// should be changed outside the constructor.

I4 Calendar::rename(const std::string InName) { // input - name of calendar

   // set error code
   I4 Err = 0;

   // set calendar name
   Name = InName;

   // Return error code
   return Err;

} // end Calendar::rename

//-------------------------------------------------------------------------
// Calendar::get - retrieves calendar properties
// Retrieves selected components of a Calendar

I4 Calendar::get(I4 *OutID,             // [out] assigned id
                 std::string *OutName,  // [out] calendar name
                 CalendarKind *OutKind, // [out] calendar type
                 I4 *OutDaysPerMonth,   // [out] array of days/month
                 I4 *OutMonthsPerYear,  // [out] months per year
                 I4 *OutSecondsPerDay,  // [out] seconds per day
                 I4 *OutSecondsPerYear, // [out] seconds per normal year
                 I4 *OutDaysPerYear     // [out] days per normal year
) const {

   // initial error code
   I4 Err = 0;

   // check for valid calendar
   if (this->CalKind == CalendarUnknown) {
      Err = 1;
      LOG_ERROR("TimeMgr: Calendar::get unknown calendar");
      return Err;
   }

   // if id requested, return id
   if (OutID != NULL)
      *OutID = this->ID;

   // if name requested, return calendar name
   if (OutName != NULL)
      *OutName = this->Name;

   // if calendar type requested, return calendar type
   if (OutKind != NULL)
      *OutKind = this->CalKind;

   // if days per month requested, return days per month
   if (OutDaysPerMonth != NULL) {
      for (I4 I = 0; I < this->MonthsPerYear; I++) {
         OutDaysPerMonth[I] = this->DaysPerMonth[I];
      }
   }

   // if months per year requested, return months per year
   if (OutMonthsPerYear != NULL)
      *OutMonthsPerYear = this->MonthsPerYear;

   // if seconds per day requested, return seconds per day
   if (OutSecondsPerDay != NULL)
      *OutSecondsPerDay = this->SecondsPerDay;

   // if seconds per year requested, return seconds per year
   if (OutSecondsPerYear != NULL)
      *OutSecondsPerYear = this->SecondsPerYear;

   // if days per year requested, return seconds per year
   if (OutDaysPerYear != NULL)
      *OutDaysPerYear = this->DaysPerYear;

   // all done, return
   return Err;

} // end Calendar::get

// Calendar constructors/destructors
//------------------------------------------------------------------------------
// Calendar::Calendar - default constructor creates no calendar option
// This constructor routine creates a calendar that uses no calendar,
// meaning the time variable stays in seconds with no notion of days
// or years. Construction by calendar kind is the most preferred
// option for construction, so this version should not be used often.

Calendar::Calendar(void) {

   // set default name and id
   Name = ' ';
   ID   = ++NumCalendars;

   // define calendar as using no calendar
   CalKind     = CalendarNoCalendar;
   CalKindName = CalendarKindName[CalKind - 1];

   MonthsPerYear = MONTHS_PER_YEAR;
   for (I4 I = 0; I < MonthsPerYear; I++)
      DaysPerMonth[I] = 0;
   SecondsPerDay  = 0;
   SecondsPerYear = 0;
   DaysPerYear    = 0;

} // end Calendar::Calendar default constructor

//------------------------------------------------------------------------------
// Calendar::Calendar - constructor creates standard calendar by kind
// This routine creates a standard calendar from a selection of
// standard supported calendars (e.g. Gregorian, Julian, no-leap).
// It is the preferred method for creating a calendar.

Calendar::Calendar(std::string InName,       // [in] name for calendar
                   CalendarKind InCalKind) { // [in] calendar type

   // set calendar name and id
   Name = InName;
   ID   = ++NumCalendars;

   // set calendar kind
   CalKind     = InCalKind;
   CalKindName = CalendarKindName[CalKind - 1];

   // months per year fixed
   MonthsPerYear = MONTHS_PER_YEAR;

   // set remaining variables based on type (kind) of calendar
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
      CalKindName = CalendarKindName[CalKind - 1];
      LOG_ERROR("TimeMgr: Calendar constructor unknown calendar kind");
      break;

   } // end switch on CalKind

} // end Calendar::Calendar standard constructor

//-----------------------------------------------------------------------------
// Calendar::Calendar - constructor that creates by copying another
// This calendar constructor creates a new calendar by copying from
// an existing calendar.

Calendar::Calendar(const Calendar &InCalendar) { // in - calendar to copy

   // copy from the input calendar
   *this = InCalendar;

   // but give it a new id
   ID = ++NumCalendars;

} // end Calendar::Calendar copy constructor

//------------------------------------------------------------------------------
// Calendar::Calendar - constructs a custom calendar
// In some cases, the standard calendars are inappropriate, so this
// custom constructor permits users to define their own calendar.

Calendar::Calendar(const std::string InName, // [in] name for calendar
                   int *InDaysPerMonth,      // [in] array of days/month
                   int InSecondsPerDay,      // [in] seconds per day
                   int InSecondsPerYear,     // [in] seconds per year
                   int InDaysPerYear) {      // [in] days per year

   // set name and id
   Name = InName;
   ID   = ++NumCalendars;

   // define calendar as custom
   CalKind     = CalendarCustom;
   CalKindName = CalendarKindName[CalKind - 1];

   // months per year is fixed
   MonthsPerYear = MONTHS_PER_YEAR;

   // set remaining calendar vars from input values
   for (I4 I = 0; I < MonthsPerYear; I++)
      DaysPerMonth[I] = InDaysPerMonth[I];
   SecondsPerDay  = InSecondsPerDay;
   SecondsPerYear = InSecondsPerYear;
   DaysPerYear    = InDaysPerYear;

} // end Calendar::Calendar custom constructor

//------------------------------------------------------------------------------
// Calendar::~Calendar - destructor for calendar

Calendar::~Calendar(void) {

   // decrement counter
   --NumCalendars;

} // end ~Calendar

// Calendar operators
//------------------------------------------------------------------------------
// Calendar(==) - Calendar equality comparison
// The equivalence operator to determine whether two calendars are the
// same. For most, it is sufficient to compare the calendar type, but
// for custom calendars, must check equivalence of each component.

bool Calendar::operator==(const Calendar &Cal) const {

   // if custom calendars, must check all properties for equality;
   // return false as soon as inequality is known, otherwise return true
   if (CalKind == CalendarCustom && Cal.CalKind == CalendarCustom) {
      if (MonthsPerYear != Cal.MonthsPerYear)
         return false;
      for (I4 I = 0; I < MonthsPerYear; I++) {
         if (DaysPerMonth[I] != Cal.DaysPerMonth[I])
            return false;
      }
      if (SecondsPerDay != Cal.SecondsPerDay)
         return false;
      if (SecondsPerYear != Cal.SecondsPerYear)
         return false;
      if (DaysPerYear != Cal.DaysPerYear)
         return false;
      return true; // custom calendars are equal
   } else {
      // for all other calendars, sufficient to just check kind
      return CalKind == Cal.CalKind;
   }
   return true;

} // end Calendar::operator==

//------------------------------------------------------------------------------
// Calendar(!=) - Calendar inequality comparison
// Compare two calendars for inequality. For most calendars, sufficient
// to check calendar kind, but custom calendars must compare all
// properties.

bool Calendar::operator!=(const Calendar &Cal) const {

   // if custom calendars, must check all properties for inequality
   // can return true as soon as any inequality is known
   if (CalKind == CalendarCustom && Cal.CalKind == CalendarCustom) {
      if (MonthsPerYear != Cal.MonthsPerYear)
         return true;
      if (SecondsPerDay != Cal.SecondsPerDay)
         return true;
      if (SecondsPerYear != Cal.SecondsPerYear)
         return true;
      if (DaysPerYear != Cal.DaysPerYear)
         return true;
      for (I4 I = 0; I < MonthsPerYear; I++) {
         if (DaysPerMonth[I] != Cal.DaysPerMonth[I])
            return true;
      }
      return false; // custom calendars are equal
   }

   // for all other calendars, just check calendar kind
   else
      return CalKind != Cal.CalKind;

} // end Calendar::operator!=

// Calendar methods
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

bool Calendar::isLeapYear(I8 Year,         // [in] a calendar year
                          I4 &Err) const { // [out] error return code

   // initialize error code
   Err = 0;

   // now check for leap year depending on calendar kind
   switch (CalKind) {
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
) const {

   // initialize error code, the basetime result, and common temps
   I4 Err{0};
   TimeFrac Result(0, 0, 1);
   I8 ResultWhole{0}; // whole seconds for result
   I8 DayOfYear{0};   // day since beginning of year
   I8 JD{0};          // Julian Day used for many conversions
   I8 HourTmp{0};     // For half-day JD conversions
   I4 FebDays{0};     // For tracking leap days

   // Branch on calendar type for actual calculation
   switch (CalKind) {

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
      if (Month != 2 && Day > DaysPerMonth[Month - 1]) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid day of month for "
                   "Gregorian calendar");
         break;
      }
      // if February, take leap year into account before checking
      //   day of the month
      if (Month == 2) {
         FebDays = DaysPerMonth[1];
         if (isLeapYear(Year, Err))
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
      ResultWhole = JD * SecondsPerDay + HourTmp * SECONDS_PER_HOUR +
                    Minute * SECONDS_PER_MINUTE + Whole;
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
      if (Month != 2 && Day > DaysPerMonth[Month - 1]) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid day of month for "
                   "Julian calendar");
         break;
      }
      // if February, take leap year into account before checking
      //   day of the month
      if (Month == 2) {
         FebDays = DaysPerMonth[1];
         if (isLeapYear(Year, Err))
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
      ResultWhole = JD * SecondsPerDay + HourTmp * SECONDS_PER_HOUR +
                    Minute * SECONDS_PER_MINUTE + Whole;
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
      ResultWhole = (Day - 1) * SecondsPerDay + Hour * SECONDS_PER_HOUR +
                    Minute * SECONDS_PER_MINUTE + Whole;
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
      if (Year < 0 || Month < 1 || Month > MonthsPerYear || Day < 1) {
         Err = 1;
         LOG_ERROR("TimeMgr: Calendar::getElapsedTime invalid date for "
                   "fixed-length calendars");
         break;
      }
      // check day of the month for any month except February
      if (Day > DaysPerMonth[Month - 1]) {
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
      ResultWhole = Year * SecondsPerYear; // secs at beg of year
      for (I4 Imonth = 0; Imonth < Month - 1; Imonth++)
         ResultWhole += DaysPerMonth[Imonth] * SecondsPerDay; // add months
      ResultWhole += (Day - 1) * SecondsPerDay;               // add days
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
) const {

   // initialize error code to success and make sure input elapsed
   // time is reduced to its simplest form
   I4 Error{0};
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
   switch (CalKind) {
   case CalendarGregorian: {
      // Reference time same as Julian Days
      // Convert whole seconds to Julian Day (number of days)
      // and remove that from elapsed seconds
      JD = Whole / SecondsPerDay;
      Whole -= (JD * SecondsPerDay);

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
      JD = Whole / SecondsPerDay;
      Whole -= JD * SecondsPerDay;

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
      JD = Whole / SecondsPerDay;
      Whole -= JD * SecondsPerDay;
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
      Year = Whole / SecondsPerYear;      // determine year
      Whole -= SecondsPerYear * Year;     // modify Whole to hold remainder
      DayInYear = Whole / SecondsPerDay;  // determine day in current year
      Whole -= SecondsPerDay * DayInYear; // modify Whole to hold remainder
      DayInYear += 1;                     // correct for 1-based counting
      // find month, day
      CountDays = 0;
      PrevDays  = 0;
      for (I4 I = 0; I < MonthsPerYear; I++) {
         CountDays += DaysPerMonth[I];
         if (DayInYear <= CountDays) { // day is in this month
            Month = I + 1;
            Day   = DayInYear - PrevDays;
            break;
         }
         PrevDays += DaysPerMonth[I];
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
) const {

   // initialize error code
   I4 Err{0};

   // Increment date based on supported interval types
   switch (Units) {
   // For year intervals, mostly just increment/decrement the year
   // while keeping all other units fixed, though check for leap day error
   case TimeUnits::Years: {
      switch (CalKind) {
      case CalendarGregorian:
      case CalendarNoLeap:
      case CalendarJulian:
      case Calendar360Day:
      case CalendarCustom: {
         Year += Interval;
         if (!(this->isLeapYear(Year, Err)) && Month == 2 && Day == 29) {
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
      switch (CalKind) {
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
         while (TmpMonth > MonthsPerYear) {
            Year += 1;
            TmpMonth -= MonthsPerYear;
         }
         // For negative intervals, check the opposite direction
         while (TmpMonth < 1) {
            Year -= 1;
            TmpMonth += MonthsPerYear;
         }
         // Set month and check that the day of the month is valid
         Month      = TmpMonth;
         I8 MaxDays = DaysPerMonth[Month - 1];
         if (this->isLeapYear(Year, Err) && Month == 2)
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
      switch (CalKind) {
      case CalendarGregorian:
      case CalendarNoLeap:
      case CalendarJulian:
      case Calendar360Day:
      case CalendarCustom: {
         if (Interval > 0) { // positive interval - step forward
            I4 DayMax = DaysPerMonth[Month - 1];
            if (Month == 2 && this->isLeapYear(Year, Err))
               DayMax += 1;
            for (I4 I = 1; I <= Interval; ++I) {
               Day += 1;
               if (Day > DayMax) { // rollover month boundary
                  Day = 1;
                  Month += 1;
                  if (Month > MonthsPerYear) { // rollover year boundary
                     Month = 1;
                     Year += 1;
                  }
                  DayMax = DaysPerMonth[Month - 1];
                  if (Month == 2 && this->isLeapYear(Year, Err))
                     DayMax += 1;
               }
            }
         } else { // negative interval - step backward
            for (I4 I = 1; I <= llabs(Interval); ++I) {
               Day -= 1;
               if (Day < 1) {
                  Month -= 1;
                  if (Month < 1) {
                     Month = MonthsPerYear;
                     Year -= 1;
                  }
                  Day = DaysPerMonth[Month - 1];
                  if (Month == 2 && this->isLeapYear(Year, Err))
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

} // namespace OMEGA
//===-----------------------------------------------------------------------===/
