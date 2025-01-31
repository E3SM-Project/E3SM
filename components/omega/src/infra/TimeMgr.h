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

#include <memory>
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

/// Utility function to extract time units from a string
TimeUnits TimeUnitsFromString(
    const std::string TimeUnitString ///< string describing time units
);

#define NUM_SUPPORTED_CALENDARS 9

/// This enum represents supported calendar types
enum CalendarKind {
   CalendarGregorian,    ///< usual Gregorian calendar
   CalendarNoLeap,       ///< Gregorian, but without leap yrs
   CalendarJulian,       ///< Julian
   CalendarJulianDay,    ///< Julian day
   CalendarModJulianDay, ///< modified Julian day
   Calendar360Day,       ///< 12 months, 30 days each
   CalendarCustom,       ///< user defined
   CalendarNoCalendar,   ///< track elapsed time only
   CalendarUnknown       ///< uninitialized or invalid
};

/// String name associated with each supported calendar type
const std::string CalendarKindName[NUM_SUPPORTED_CALENDARS] = {
    "Gregorian", "No Leap", "Julian",      "Julian Day", "Modified Julian Day",
    "360 Day",   "Custom",  "No Calendar", "Invalid"};

/// Standard CF-compliant name associated with each supported calendar type
const std::string CalendarCFName[NUM_SUPPORTED_CALENDARS] = {
    "gregorian", "noleap", "julian", "julian_day", "modified_julian_day",
    "360_day",   "custom", "none",   "invalid"};

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

/// The Calendar class encapsulates the knowledge (attributes and behavior)
/// of all supported calendar kinds:  Gregorian, Julian, Julian Day, Modified
/// Julian Day, no-leap, 360-day, custom, and no-calendar.
///
/// Time in OMEGA is measured in seconds since a reference time so conversions
/// between that time and actual calendar dates must convert correctly based on
/// calendar assumptions. For the fixed calendars like 360-day, no-leap, and
/// custom calendars, the reference time is simply year 0 at midnight (where
/// midnight typically refers to Universal Time). For the other standard
/// calendars, the reference time is part of the calendar definition.
///
/// For simplicity, Julian Day is often used during calendar conversion to
/// time. The Gregorian to/from Julian day conversion algorithm is from
/// Henry F. Fliegel and Thomas C. Van Flandern, in Communications of the ACM,
/// CACM, volume 11, number 10, October 1968, p. 657. Julian day refers to the
/// number of days since a reference day. For the algorithm used, this reference
/// day is noon November 24, -4713 in the Proleptic Gregorian calendar
/// (Gregorian calendar reversed in time to that date), which is equivalent to
/// January 1, -4712 in the Proleptic Julian calendar. When converting from a
/// Julian day to a Gregorian date, this algorithm is valid from 3/1/-4900
/// Gregorian forward.  When converting from a Gregorian date to a Julian day,
/// the algorithm is valid from 3/1/-4800 forward. In both cases, the algorithm
/// correctly takes into account leap years, those that are divisible by 4
/// and not 100, or those divisible by 400.
///
class Calendar {
   // private variables
 private:
   static std::unique_ptr<Calendar> OmegaCal; ///< single instance of calendar
   CalendarKind CalKind;                      ///< enum for calendar kind
   std::string CalKindName;                   ///< name of calendar kind

   // variables defining calendar characteristics for time
   std::vector<I4> DaysPerMonth; ///< days in each month
   I4 MonthsPerYear;             ///< num months in year
   I4 SecondsPerDay;             ///< seconds per day
   I4 SecondsPerYear;            ///< seconds per normal year
   I4 DaysPerYear;               ///< days per normal year

   // Constructors are private - use set to create the only calendar
   /// Constructor based on kind of calendar
   Calendar(CalendarKind CalKind ///< [in] choice of calendar kind
   );
   /// Constructs custom calendar based in inputs
   Calendar(std::vector<I4> &InDaysPerMonth, ///< [in] array of days per month
            I4 InSecondsPerDay,              ///< [in] seconds per day
            I4 InSecondsPerYear,             ///< [in] seconds per year
            I4 InDaysPerYear                 ///< [in] days per year (dpy)
   );

   // public methods
 public:
   /// Initializes a standard model calendar
   static void
   init(std::string CalendarKindStr ///< [in] string for type of calendar
   );

   /// Creates a custom user-defined calendar
   static void
   init(std::vector<I4> &InDaysPerMonth, ///< [in] array of days per month
        I4 InSecondsPerDay,              ///< [in] seconds per day
        I4 InSecondsPerYear,             ///< [in] seconds per year
        I4 InDaysPerYear                 ///< [in] days per year (dpy)
   );

   /// Checks whether a calendar has been defined
   static bool isDefined();

   /// Retrieve pointer to the calendar
   static Calendar *get();

   /// Retrieve Type of calendar
   static CalendarKind getKind();

   /// Retrieve days per month in calendar
   static std::vector<I4> getDaysPerMonth();

   /// Retrieve months per year in calendar
   static I4 getMonthsPerYear();

   /// Retrieve days per month in calendar
   static I4 getSecondsPerDay();

   /// Retrieve days per month in calendar
   static I4 getSecondsPerYear();

   /// Retrieve days per year in calendar
   static I4 getDaysPerYear();

   /// Disable copy constructor
   Calendar(const Calendar &) = delete;
   Calendar(Calendar &&)      = delete;

   /// Destroy the single calendar instance
   /// This should only be used during testing as it will invalidate
   /// all time instants and other behavior
   static void reset();

   /// Calendar destructor
   ~Calendar(void);

   /// Validate calendar
   I4 validate() const;

   /// Checks whether input year is a leap year
   /// \return true if year is a leap year, false otherwise
   static bool isLeapYear(I8 Year ///< [in] year to check
   );

   /// Computes the total elapsed time in seconds (in TimeFrac form)
   /// since the calendar reference time, given a calendar date, time.
   /// \return Elapsed time in TimeFrac form
   static TimeFrac
   getElapsedTime(const I8 Year,   ///< [in] calendar year
                  const I8 Month,  ///< [in] calendar month
                  const I8 Day,    ///< [in] calendar day
                  const I8 Hour,   ///< [in] time of day-hour
                  const I8 Minute, ///< [in] time of day-min
                  const I8 Whole,  ///< [in] time of day-whole seconds
                  const I8 Numer,  ///< [in] time of day-frac secs (numerator)
                  const I8 denom   ///< [in] time of day-frac secs (denom)
   );

   /// Determines the calendar date and time of day, given an
   /// elapsed time since the calendar reference time.
   /// \return error code
   static I4
   getDateTime(const TimeFrac ElapsedTime, ///< [in] time in secs from ref time
               I8 &Year,                   ///< [out] calendar year
               I8 &Month,                  ///< [out] calendar month
               I8 &Day,                    ///< [out] calendar day
               I8 &Hour,                   ///< [out] time of day-hours
               I8 &Minute,                 ///< [out] time of day-minutes
               I8 &Whole,                  ///< [out] time of day-whole seconds
               I8 &Numer, ///< [out] time of day-frac secs (numerator)
               I8 &Denom  ///< [out] time of day-frac secs (denom)
   );

   /// Increments (or decrements) a calendar date by a specified
   /// interval, supplied by an integer interval in given time units.
   /// Only calendar based intervals (years, months or days) are
   /// supported. This is primarily meant to be called by other time
   /// manager routines (eg to add/subtract time instants) for those
   /// time intervals that are dependent on date and sensitive to
   /// calendar features like leap years and varying days of the month.
   /// \return error code
   static I4
   incrementDate(const I8 Interval,     ///< [in] time interval to advance date
                 const TimeUnits Units, ///< [in] time units for interval
                 I8 &Year,  ///< [in,out] calendar year for time to be advanced
                 I8 &Month, ///< [in,out] calendar month for time to be advanced
                 I8 &Day    ///< [in,out] calendar day for time to be advanced
   );
}; // end class Calendar

/// The TimeInterval class represents an interval of time -- the amount of time
/// between two instants in time.  These can either be independent of any
/// calendar (represented as fractional seconds) or dependent on a calendar
/// (eg months or years) and thought of as a calendar interval.
///
/// Calendar independent intervals are stored as fractional integer seconds and
/// use much of the functionality from the TimeFrac class (though inheritance
/// proved awkward). Interfaces creating and  accessing calendar-independent
/// intervals also support units of hours, minutes and seconds in either real
/// or integer form.
///
/// A TimeInterval can also be defined as a Calendar interval and then becomes
/// calendar-dependent. Currently, these must be specified in integer years,
/// months or days and all operations that compare or or modify these intervals
/// must be in the same units of years, months or days. Note that this
/// restriction only applies to the operations described below between time
/// intervals. When interacting with time instants, these restrictions do
/// not apply.
///
/// TimeInterval also defines methods for multiplication and division of
/// TimeIntervals by integer and real scalars. Other methods include absolute
/// value and negative absolute value for use with both positive or negative
/// time intervals.
///
class TimeInterval {
   // private variables
 private:
   TimeFrac Interval; ///< Non-calendar interval in fractional seconds
   bool IsCalendar;   ///< True if calendar interval
   I8 CalInterval;    ///< Calendar interval length
   TimeUnits Units;   ///< Calendar interval units

 public:
   // constructors/destructors

   /// Default time interval constructor
   TimeInterval(void);

   /// Construct time interval from base time fractional integer seconds
   TimeInterval(const I8 Whole, ///< Whole integer seconds
                const I8 Numer, ///< Fractional seconds numerator
                const I8 Denom  ///< Fractional seconds denominator
   );

   /// Construct time interval from an integer length and units
   TimeInterval(const I4 InLength,      ///< length of time interval
                const TimeUnits InUnits ///< unit of time for interval
   );

   /// Construct time interval from a I8 integer length and units
   TimeInterval(const I8 Length,        ///< length of time interval
                const TimeUnits InUnits ///< unit of time for interval
   );

   /// Construct time interval from an real length and units
   TimeInterval(const R8 Length,        ///< length of time interval
                const TimeUnits InUnits ///< unit of time for interval
   );

   /// Construct time interval from a standard string in the form
   /// DDDD_HH:MM:SS.SSSS where the width of DD and SS strings can be of
   /// arbitrary width (within reason) and the separators can be any single
   /// non-numeric character
   TimeInterval(std::string &TimeString ///< [in] string form of time interval
   );

   /// Destructor for time interval
   ~TimeInterval(void);

   // Accessor methods

   /// Set a non-calendar interval in native fractional integer seconds
   /// \return error code
   I4 set(const I8 Whole, ///< Whole integer seconds
          const I8 Numer, ///< Fractional seconds numerator
          const I8 Denom  ///< Fractional seconds denominator
   );

   /// Set a time interval from I4 length and units
   /// \return error code
   I4 set(const I4 InLength,      ///< length of time interval
          const TimeUnits InUnits ///< unit of time for interval
   );

   /// Set a time interval from I8 length and units
   /// \return error code
   I4 set(const I8 Length,        ///< length of time interval
          const TimeUnits InUnits ///< unit of time for interval
   );

   /// Set a time interval from a real length and units
   /// Real length is only supported for non-calendar intervals since
   /// a non-integral length has ambiguous meaning when months, years
   /// have variable length.
   /// \return error code
   I4 set(const R8 Length,        ///< length of time interval
          const TimeUnits InUnits ///< unit of time for interval
   );

   /// Retrieve non-calendar interval in native fractional integer form
   /// \return error code
   I4 get(I8 &Whole, ///< [out] whole seconds
          I8 &Numer, ///< [out] fractional second numerator
          I8 &Denom  ///< [out] fractional second denominator
   ) const;

   /// Retrieve a time interval in integer length in specified units.
   /// To avoid roundoff issues during conversions, integer retrieval
   /// is only permitted in the same units in which the interval was
   /// defined.
   /// \return error code
   I4 get(I8 &Length, ///< [out] requested integer length of interval
          const TimeUnits ReqUnits ///< [in] unit of time for interval
   ) const;

   /// Retrieve a time interval in real length and specified units.
   /// For calendar intervals, the units must match the units in which
   /// the interval was defined. For non-calendar intervals, only conversions
   /// between hours, minutes and seconds are allowed.
   /// \return error code
   I4 get(R8 &Length,              ///< [out] Requested time interval length
          const TimeUnits ReqUnits ///< [in] unit of time for interval
   ) const;

   /// Equivalence comparison operator for TimeInterval
   bool operator==(const TimeInterval &) const;
   /// Non-equivalence comparison operator for TimeInterval
   bool operator!=(const TimeInterval &) const;
   /// Less than comparison operator for TimeInterval
   bool operator<(const TimeInterval &) const;
   /// Greater than comparison operator for TimeInterval
   bool operator>(const TimeInterval &) const;
   /// Less than or equal comparison operator for TimeInterval
   bool operator<=(const TimeInterval &) const;
   /// Greater than or equal comparison operator for TimeInterval
   bool operator>=(const TimeInterval &) const;
   /// Addition operator for TimeInterval
   TimeInterval operator+(const TimeInterval &) const;
   /// Subtraction operator for TimeInterval
   TimeInterval operator-(const TimeInterval &) const;
   /// Increment operator for TimeInterval
   TimeInterval &operator+=(const TimeInterval &);
   /// Decrement operator for TimeInterval
   TimeInterval &operator-=(const TimeInterval &);
   /// Multiplication by integer scalar
   TimeInterval operator*(const I4 Multiplier) const;
   /// Multiplication by real scalar
   /// This is primarily meant for non-calendar intervals, but
   /// for calendar intervals, it will multiply the year, month or
   /// day interval and convert to the nearest integer.
   TimeInterval operator*(const R8 Multiplier) const;
   /// Multiplication in place by integer scalar
   TimeInterval &operator*=(const I4 Multiplier);
   /// Multiplication in place by real scalar
   /// This is primarily meant for non-calendar intervals, but
   /// for calendar intervals, it will multiply the year, month or
   /// day interval and convert to the nearest integer.
   TimeInterval &operator*=(const R8 Multiplier);
   /// Divide interval by integer scalar
   TimeInterval operator/(const I4 Divisor) const;
   /// Divide interval in place by integer scalar
   TimeInterval &operator/=(const I4 Divisor);

   /// Absolute value
   static TimeInterval absValue(const TimeInterval &);
   /// Negative absolute value
   static TimeInterval negAbsValue(const TimeInterval &);

   // Other utility methods
   /// Check whether a time interval is positive
   bool isPositive(void);

   /// commutative multiplication operators need to be defined as
   /// free functions, and therefore need to be given acces to
   /// private members of TimeInterval
   friend TimeInterval operator*(const I4 &Multiplier, const TimeInterval &TI);
   friend TimeInterval operator*(const R8 &Multiplier, const TimeInterval &TI);

   /// Give TimeInstant access to TimeInterval so that
   /// incrementing/decrementing time is easier.
   friend class TimeInstant;

}; // end class TimeInterval

/// The TimeInstant class represents an instant of time. It is represented by
/// a time in fractional seconds from a reference time, typically year or day
/// zero of the relevant calendar.
///
/// This class includes a number of functions for manipulating time for the
/// purpose of tracking time within a simulation. These can include adding a
/// TimeInterval (eg a time step) to advance the time or creating a TimeInterval
/// by subtracting two instants in time. Other useful queries like time of day
/// or day of year are also provided to place the time in a practical context.
///
class TimeInstant {
   // private variables
 private:
   TimeFrac ElapsedTime; ///< Fractional seconds since reference time

 public:
   // constructors/destructors

   /// Default constructor creates empty time instant
   TimeInstant(void);

   /// Construct time instant from date, time, calendar
   /// Where seconds is supplied as real number.
   TimeInstant(const I8 Year,   ///< [in] year
               const I8 Month,  ///< [in] month
               const I8 Day,    ///< [in] day
               const I8 Hour,   ///< [in] hour
               const I8 Minute, ///< [in] minute
               const R8 RSecond ///< [in] second (real)
   );

   /// Construct time instant from date, time, calendar
   /// Where seconds is supplied in integer fractional seconds.
   TimeInstant(const I8 Year,   ///< [in] year
               const I8 Month,  ///< [in] month
               const I8 Day,    ///< [in] day
               const I8 Hour,   ///< [in] hour
               const I8 Minute, ///< [in] minute
               const I8 Whole,  ///< [in] second (whole integer)
               const I8 Numer,  ///< [in] second (fraction numerator)
               const I8 Denom   ///< [in] second (fraction denominator)
   );

   /// Construct time instant from a standard date-time string in the
   /// form YYYY-MM-DD_HH:MM:SS.SSSS where the width of the YY and SS
   /// strings can be of arbitrary width (within reason) and the
   /// separators can be any single non-numeric character
   TimeInstant(std::string &TimeString ///< [in] string containing date/time
   );

   /// Destructor for time interval
   ~TimeInstant(void);

   // Accessor methods

   /// Set time instant from date and time, where seconds is supplied
   /// as a real number.
   /// \return error code
   I4 set(const I8 Year,   ///< [in] year
          const I8 Month,  ///< [in] month
          const I8 Day,    ///< [in] day
          const I8 Hour,   ///< [in] hour
          const I8 Minute, ///< [in] minute
          const R8 RSecond ///< [in] second (real)
   );

   /// Set time instant from date and time, where seconds is supplied
   /// in integer fractional seconds.
   /// \return error code
   I4 set(const I8 Year,   ///< [in] year
          const I8 Month,  ///< [in] month
          const I8 Day,    ///< [in] day
          const I8 Hour,   ///< [in] hour
          const I8 Minute, ///< [in] minute
          const I8 Whole,  ///< [in] second (whole integer)
          const I8 Numer,  ///< [in] second (fraction numerator)
          const I8 Denom   ///< [in] second (fraction denominator)
   );

   /// Retrieve time in date, time form with real seconds.
   /// \return error code
   I4 get(I8 &Year,   ///< [out] year   of this time instant
          I8 &Month,  ///< [out] month  of this time instant
          I8 &Day,    ///< [out] day    of this time instant
          I8 &Hour,   ///< [out] hour   of this time instant
          I8 &Minute, ///< [out] minute of this time instant
          R8 &Second  ///< [out] second of this time instant
   ) const;

   /// Retrieve time in date, time form with fractional integer seconds.
   /// \return error code
   I4 get(I8 &Year,   ///< [out] year   of this time instant
          I8 &Month,  ///< [out] month  of this time instant
          I8 &Day,    ///< [out] day    of this time instant
          I8 &Hour,   ///< [out] hour   of this time instant
          I8 &Minute, ///< [out] minute of this time instant
          I8 &Whole,  ///< [out] whole seconds of this time
          I8 &Numer,  ///< [out] frac second numerator
          I8 &Denom   ///< [out] frac second denominator
   ) const;

   // Operators on time instants

   /// Equivalence comparison for TimeInstant
   bool operator==(const TimeInstant &) const;
   /// Non-equivalence comparison operator for TimeInstant
   bool operator!=(const TimeInstant &) const;
   /// Less than comparison operator for TimeInstant
   bool operator<(const TimeInstant &) const;
   /// Greater than comparison operator for TimeInstant
   bool operator>(const TimeInstant &) const;
   /// Less than or equal comparison operator for TimeInstant
   bool operator<=(const TimeInstant &) const;
   /// Greater than or equal comparison operator for TimeInstant
   bool operator>=(const TimeInstant &) const;
   /// Increment time by adding a time interval
   TimeInstant operator+(const TimeInterval &) const;
   /// Decrement time by subtracting a time interval
   TimeInstant operator-(const TimeInterval &) const;
   /// Create a time interval by subtracting two time instants
   TimeInterval operator-(const TimeInstant &) const;
   /// Increment time in place by adding time interval
   TimeInstant &operator+=(const TimeInterval &);
   /// Decrement time in place by subtracting time interval
   TimeInstant &operator-=(const TimeInterval &);

   // Other utility methods
   /// Get time as a string in the format
   /// 'YYYYYY-MM-DD{separator}HH:MM:SS.SSSSSS' where the number
   /// of digits in year, number of digits after the decimal in
   /// seconds and the character to use as separator are all input
   /// by the user.
   /// \return time string
   std::string
   getString(const I4 YearWidth,   ///< [in] number of digits in year
             const I4 SecondWidth, ///< [in] num digits after decimal in seconds
             std::string Separator ///< [in] string(char) to separate date/time
   ) const;

}; // end class TimeInstant

/// The Alarm class mimics a clock alarm by ringing on or after a specified
/// time. Alarms can be created to ring after a single time or to ring on
/// periodic intervals. Once ringing, an alarm must be manually stopped using
/// either a stop method or a reset. The reset method will stop the current
/// alarm and set the new alarm time for periodic/interval alarms.
///
class Alarm {
   // private variables
 private:
   std::string Name; ///< name for the alarm

   bool Ringing;  ///< alarm is currently ringing
   bool Periodic; ///< alarm rings periodically on interval
   bool Stopped;  ///< alarm has been stopped and not reset

   TimeInstant RingTime;      ///< time at/after which alarm rings
   TimeInterval RingInterval; ///< interval at which this alarm rings
   TimeInstant RingTimePrev;  ///< previous alarm time for interval alarms

 public:
   // constructors/destructors

   /// Default alarm constructor
   Alarm(void);

   /// Constructs a one-time alarm using the input ring time.
   Alarm(const std::string InName,   ///< [in] Name of alarm
         const TimeInstant AlarmTime ///< [in] Time at/after which alarm rings
   );

   /// Constructs a periodic/interval alarm based on an input periodic
   /// time interval and a start time for the interval period.
   Alarm(
       const std::string InName,         ///< [in] Name of alarm
       const TimeInterval AlarmInterval, ///< [in] interval at which alarm rings
       const TimeInstant IntervalStart   ///< [in] start time of first interval
   );

   /// Destructor for alarm
   ~Alarm(void);

   /// Check whether an alarm is ringing
   /// \return true if alarm is ringing, false otherwise
   bool isRinging(void);

   /// Checks whether the alarm should ring based on the current
   /// (or supplied) time instant (returns error code)
   I4 updateStatus(const TimeInstant CurrentTime ///< [in] current time
   );

   /// Stops a ringing alarm and sets next a new ring time. If the alarm
   /// is a periodic/interval alarm, the next ring time is set to be the
   /// next interval boundary after the input time.  If the alarm is a
   /// single instance, the input time is used as the next alarm time.
   /// (returns error code)
   I4 reset(const TimeInstant InTime ///< [in] new basis for alarm time
   );

   /// Stops a ringing alarm (returns error code).
   I4 stop(void);

   /// Rename an alarm (returns error code)
   I4 rename(const std::string NewName ///< [in] new name for alarm
   );

   /// Get alarm name
   std::string getName(void) const;

   /// Get alarm interval
   const TimeInterval *getInterval(void) const;

   /// Get last time the alarm rang
   const TimeInstant *getRingTimePrev(void) const;

}; // end class Alarm

// default max number of alarms to create initial space in Alarms vector
#define MAX_ALARMS 20

/// The Clock class manages the time for a forward-integrating model. Given a
/// start time and a time step, a user can advance the clock. The clock will be
/// consistent with the Calendar associated with the start time. A user can
/// attach any number of alarms to the clock and the clock will update the
/// ringing status of these attached alarms as the clock marches forward.
///
class Clock {
   // private variables
 private:
   TimeInstant StartTime; ///< initial time for this clock
   TimeInstant CurrTime;  ///< current time
   TimeInstant PrevTime;  ///< time at previous timestep
   TimeInstant NextTime;  ///< time at next timestep
   TimeInterval TimeStep; ///< interval at which this clock advances

   I4 NumAlarms; ///< current number of attached alarms

   std::vector<Alarm *>
       Alarms; ///< pointers to alarms associated with this clock

 public:
   // constructors/destructors

   /// Construct a clock from start time and time step
   Clock(const TimeInstant StartTime, ///< [in] Start time for clock
         const TimeInterval TimeStep  ///< [in] Time step to advance clock
   );

   /// Destructor for clocks
   ~Clock(void);

   // Accessor Methods
   /// Set the current time (returns error code)
   I4 setCurrentTime(
       const TimeInstant CurrTime ///< [in] new value for current time
   );

   /// Changes the time step for this clock (returns error code)
   I4 changeTimeStep(
       const TimeInterval NewTimeStep ///< [in] new value for time step
   );

   /// Retrieves current time of this clock
   TimeInstant getCurrentTime(void) const;

   /// Retrieves time at the previous time step
   TimeInstant getPreviousTime(void) const;

   /// Retrieves time at the next time step
   TimeInstant getNextTime(void) const;

   /// Retrieves start time for this clock
   TimeInstant getStartTime(void) const;

   /// Retrieves time step for clock
   TimeInterval getTimeStep(void) const;

   /// Attaches an alarm to this clock. The clock simply stores a
   /// pointer to this alarm. (returns error code)
   I4 attachAlarm(Alarm *InAlarm ///< [in] pointer to alarm to attach
   );

   /// Advance a clock one timestep and update status of any attached
   /// alarms. (returns error code)
   I4 advance(void);

}; // end class Clock

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_TIMEMGR_H
