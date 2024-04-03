#ifndef OMEGA_TIMEMGR_H
#define OMEGA_TIMEMGR_H
//===-- infra/TimeMgr.h - definitions for time manager ----------*- C++ -*-===//
//
/// \file
/// /brief This file defines classes that comprise the time manager for OMEGA
/// and declares associated methods.
///
/// The OMEGA time manager tracks simulation time during a forward-integration
/// in a manner to avoid accumulated roundoff in simulations that can extend
/// over millions of time steps, is compatible with numerous calendar systems,
/// and manages alarms for triggering events at precise times. This time manager
/// module consists of six classes:
/// 1. The TimeFrac class is the base time representation for a number of
///    seconds as an integer fraction
/// 2. The Calendar class holds important information for supported
///    calendar kinds
/// 3. The TimeInstant class represents a point in time on a given calendar
/// 4. The TimeInterval class represents the amount of time between different
///    points in time
/// 5. The Alarm class allows the user to trigger events with single-use or
///    periodic alarms
/// 6. The Clock class tracks time during an integration and helps manage
///    attached alarms
///
/// Based in part on ESMF time manager, so copyright attached here:
/// Earth System Modeling Framework
/// Copyright 2002-2024, University Corporation for Atmospheric Research,
/// Massachusetts Institute of Technology, Geophysical Fluid Dynamics
/// Laboratory, University of Michigan, National Centers for Environmental
/// Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
/// NASA Goddard Space Flight Center.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"

#include <string>

// Definitions of conversions
/// Define seconds per day
#define SECONDS_PER_DAY 86400
/// Define seconds per hour
#define SECONDS_PER_HOUR 3600
/// Define seconds per minute
#define SECONDS_PER_MINUTE 60
/// Define months per year
#define MONTHS_PER_YEAR 12

namespace OMEGA {

/// This enum class defines time units used throughout the OMEGA
/// time manager
enum class TimeUnits {
   None = 0, ///< value for undefined units
   Seconds,  ///< time units in seconds (typical)
   Minutes,  ///< time units in minutes
   Hours,    ///< time units in hours
   Days,     ///< time units in days
   Months,   ///< time units in months
   Years,    ///< time units in years
};

#define NUM_SUPPORTED_CALENDARS 9

/// This enum represents supported calendar types
enum CalendarKind {
   CalendarGregorian = 1, ///< usual Gregorian calendar
   CalendarNoLeap,        ///< Gregorian, but without leap yrs
   CalendarJulian,        ///< Julian
   CalendarJulianDay,     ///< Julian day
   CalendarModJulianDay,  ///< modified Julian day
   Calendar360Day,        ///< 12 months, 30 days each
   CalendarCustom,        ///< user defined
   CalendarNoCalendar,    ///< track elapsed time only
   CalendarUnknown        ///< uninitialized or invalid
};

/// String name associated with each supported calendar type
const std::string CalendarKindName[NUM_SUPPORTED_CALENDARS] = {
    "Gregorian", "No Leap", "Julian",      "Julian Day", "Modified Julian Day",
    "360 Day",   "Custom",  "No Calendar", "Invalid"};

/// The TimeFrac class includes the core representation, functions and
/// operators for time within OMEGA.
///
/// The TimeFrac class represents time as an integer fraction number of
/// seconds. The integer representation is to enable accumulation over
/// long times without accumulating roundoff. However, interfaces are provided
/// for accessing and creating a time using more user-friendly units of hours
/// minutes, seconds, etc.  To handle the large time range, a 64-bit integer
/// is used. Any rational fractional second is expressed using two additional
/// integers: a numerator and a denominator.  Both the whole seconds and
/// fractional numerator are signed to handle negative time intervals and
/// instants. For arithmetic consistency both must carry the same sign (both
/// positive or both negative), except, of course, for zero values. The
/// fractional seconds element (numerator) is normalized (bounded) with respect
/// to whole seconds. If the  absolute value of the numerator becomes greater
/// than or equal to the denominator, the whole seconds is incremented or
/// decremented accordingly and the numerator is reset to the remainder.
///
/// For ease in calendar conversions, a time value of zero (both whole
/// and numerator) will correspond to the Julian date of zero.
///
class TimeFrac {
   // private variables
 private:
   I8 Whole; ///< Integer (whole) seconds
   I8 Numer; ///< Integer fractional second (n/d) numerator
   I8 Denom; ///< Integer fractional second (n/d) denominator

 public:
   // Accessor methods

   /// Single call to set all native base time components
   /// \return error code
   I4 set(const I8 W, ///< [in] whole seconds
          const I8 N, ///< [in] fractional second numerator
          const I8 D  ///< [in] fractional second denominator
   );
   /// Set base time by converting from integer hours, minutes, seconds
   /// \return error code
   I4 setHMS(const I4 Hours,   ///< [in] integer hours
             const I4 Minutes, ///< [in] integer minutes
             const I4 Seconds  ///< [in] integer seconds
   );

   /// Set base time by converting from a real number of seconds
   /// \return error code
   I4 setSeconds(const R8 Seconds ///< [in] Time in real seconds
   );
   /// Set base time by converting from a real number of hours
   /// \return error code
   I4 setHours(const R8 Hours ///< [in] Time in real hours
   );
   /// Set base time by converting from a real number of minutes
   /// \return error code
   I4 setMinutes(const R8 Minutes ///< [in] Time in real minutes
   );

   /// Set whole seconds separately
   /// \return error code
   I4 setWhole(const I8 W ///< [in] Whole number of seconds
   );
   /// Set numerator of fractional seconds separately
   /// \return error code
   I4 setNumer(const I8 N ///< [in] Numerator of fractional seconds
   );
   /// Set denominator of fractional seconds separately
   /// \return error code
   I4 setDenom(const I8 D ///< [in] Denominator of fractional seconds
   );

   /// Single call to retrieve native base time components
   /// \return error code
   I4 get(I8 &W, ///< [out] whole seconds
          I8 &N, ///< [out] fractional second numerator
          I8 &D  ///< [out] fractional second denominator
   ) const;
   /// Get base time converted to integer hours, minutes, seconds
   /// \return error code
   I4 getHMS(I4 &Hours,   ///< [out] integer hours
             I4 &Minutes, ///< [out] integer minutes
             I4 &Seconds  ///< [out] integer seconds
   ) const;
   /// Get base time and convert to a real number of seconds
   /// \return Time in real seconds
   R8 getSeconds(void) const;
   /// Get base time and convert to a real number of hours
   /// \return Time in real hours
   R8 getHours(void) const;
   /// Get base time and convert to a real number of minutes
   /// \return Time in real minutes
   R8 getMinutes(void) const;
   /// Retrieve the whole seconds component of base time
   /// \return Whole number of seconds
   I8 getWhole(void) const;
   /// Retrieve the numerator component of fractional base time
   /// \return Numerator of fractional seconds
   I8 getNumer(void) const;
   /// Retrieve the denominator component of fractional base time
   /// \return Denominator of fractional seconds
   I8 getDenom(void) const;

   // constructors/destructors

   /// Default base time constructor
   TimeFrac(void);
   /// Copy constructor for base time
   TimeFrac(const TimeFrac & ///< [in] existing base time to be copied
   );
   /// Construct base time by component
   TimeFrac(const I8 Whole, ///< [in] whole seconds
            const I8 Numer, ///< [in] fractional second numerator
            const I8 Denom  ///< [in] fractional second denominator
   );
   /// Construct base time by converting from a real number of seconds
   TimeFrac(const R8 Seconds ///< [in] Time in real seconds
   );
   /// Destructor for base time
   ~TimeFrac(void);

   // operators

   /// Equivalence comparison operator for TimeFrac
   bool operator==(const TimeFrac &) const;
   /// Non-equivalence comparison operator for TimeFrac
   bool operator!=(const TimeFrac &) const;
   /// Less than comparison operator for TimeFrac
   bool operator<(const TimeFrac &) const;
   /// Greater than comparison operator for TimeFrac
   bool operator>(const TimeFrac &) const;
   /// Less than or equal comparison operator for TimeFrac
   bool operator<=(const TimeFrac &) const;
   /// Greater than or equal comparison operator for TimeFrac
   bool operator>=(const TimeFrac &) const;
   /// Addition operator for TimeFrac
   TimeFrac operator+(const TimeFrac &) const;
   /// Subtraction operator for TimeFrac
   TimeFrac operator-(const TimeFrac &) const;
   /// Increment operator for TimeFrac
   TimeFrac &operator+=(const TimeFrac &);
   /// Decrement operator for TimeFrac
   TimeFrac &operator-=(const TimeFrac &);
   /// Multiplication by integer scalar
   TimeFrac operator*(const I4 Multiplier) const;
   /// Multiplication in place by integer scalar
   TimeFrac &operator*=(const I4 Multiplier);
   /// Multiplication by real scalar
   TimeFrac operator*(const R8 Multiplier) const;
   /// Multiplication in place by real scalar
   TimeFrac &operator*=(const R8 Multiplier);
   /// Divide TimeFrac by integer scalar
   TimeFrac operator/(const I4 Divisor) const;
   /// Divide TimeFrac in place by integer scalar
   TimeFrac &operator/=(const I4 Divisor);
   /// Divide two TimeFracs and return a real result
   R8 operator/(const TimeFrac &) const;
   /// Modulus method for TimeFrac
   TimeFrac operator%(const TimeFrac &) const;
   /// Modulus method in place
   TimeFrac &operator%=(const TimeFrac &);
   /// Assignment operator for TimeFrac
   TimeFrac &operator=(const TimeFrac &);

   // Other utility methods

   /// Convert a time fraction to new denominator
   /// \return error code
   I4 convert(const I8 Denom); ///< [in] new denominator
   /// Reduce a time fraction to simplest form
   I4 simplify(void);

}; // end class TimeFrac

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_TIMEMGR_H
