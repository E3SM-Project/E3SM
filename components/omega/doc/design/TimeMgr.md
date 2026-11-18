(omega-design-time-manager)=
# TimeManager

## 1 Overview

When simulating the Earth System, a model must track the simulation time
as it steps forward with a given time step. Time must be kept in the
appropriate calendar system, accounting for leap years and other
adjustments as needed. Because these simulations may extend over
millions of time steps or more, time must be accumulated without roundoff
so that accumulated time does not drift and events will occur at precise
time intervals. This document describes various classes needed to track the
model's simulation time accurately and manage events in the appropriate
calendar conventions.

## 2 Requirements

### 2.1 Required: Simulation time

The time manager must be able to track the current simulation time
as the model integration steps forward.

### 2.2 Required: No accumulation roundoff

Because integrations may extend over millions of time steps or
longer, the model must track time without accumulating roundoff.

### 2.3 Required: Calendar support

The time manager must be able to track time in supported calendars,
starting with:
- Gregorian
- Gregorian with no leap year
- Idealized 360-day calendar with equal 30-day months
- No calendar (some test cases only need elapsed time in sec)
Time for standard calendars will assume Coordinated Universal
Time (UTC).

### 2.4 Required: Time intervals

The model must be able to compute time intervals, not only to
represent the model time step itself, but also to be able to
compute time between events, track periodic events or perform
time averaging. These intervals can be both absolute time (seconds)
or represented as calendar intervals (monthly, yearly, number of
days). Conversion from a calendar interval to absolute interval
(in seconds) will be needed to compute time-related quantities
like time averaging.

### 2.5 Required: Events or alarms

The time manager must be able to trigger events that occur
at specific times or at specified periodic time intervals.
This requires being able to compare time instants and trigger an
alarm if a time instant or interval has been reached. It will
also be useful to track when a previous or next event occurs in a
periodic sequence.

### 2.6 Required: Long times

To support very long time integrations or potential paleoclimate
simulations, the model must support years that extend beyond the
4-digit representation and support negative years to represent
years Before the Common Era (BCE)

### 2.7 Required: Time string representation

For model output or metadata, the model must be able to
represent a time in accepted formats. Following ISO standards
(eg ISO 8601 and prior), these will typically represent time
in decreasing intervals, eg `YYYYYY-MM-DD HH:MM:SS.SSSSS` though
variations on this format can be supported to use, for example,
different separators and different widths for years and fractional
seconds.

## 3 Algorithmic Formulation

The requirement for no roundoff accumulation implies the need to
use an integer fraction representation as in the ESMF time
manager. Similarly, calendar days are best tracked using the
Julian Day that counts days since a specified start time
(Noon UTC on 1 Jan, 4713 BCE in the Gregorian calendar -
a beautiful Monday by most accounts) with conversion to
various calendars following the algorithms in:

Fliegel, H. F., and Van Flandern, T. C., 1968, A Machine Algorithm
for Processing Calendar Dates, Communications of the Association of
Computing Machines, 11, 657.

Hatcher, D.A., 1984. Simple Formulae for Julian Day Numbers and
Calendar Dates, Quart. J. of R. Astr. Soc., 25, 53-55.

## 4 Design

The time manager design will follow an ESMF-like approach,
though simplified and not requiring linking to the full
ESMF library. Previous simplifications were performed for
WRF by John Michalakes in Fortran and are used in MPAS, E3SM.
A similar simplification was performed by Phil Jones from the
ESMF C++ version for the SciDAC CANGA project and will be used
here.

The time manager will consist of a number of classes/modules:
- `TimeFrac`: a base fraction representation
- `TimeInstant`: representation of a point in time
- `TimeInterval`: a time step or difference between time instants
- `Calendar`: support for various calendars
- `Alarm`: alarms that trigger at time instants or periodic intervals
- `Clock`: a clock that keeps track of model time as it marches forward

### 4.1 Data types and parameters

#### 4.1.1 Parameters

A number of parameters are defined among the above classes. For
all classes, it will be useful to define times and intervals with
a number and units, so we define an enum class for time units:

```c++
enum class TimeUnits{
                 None = 0, ///< value for undefined units
                 Seconds,  ///< time units in seconds (typical)
                 Minutes,  ///< time units in minutes
                 Hours,    ///< time units in hours
                 Days,     ///< time units in days
                 Months,   ///< time units in months
                 Years,    ///< time units in years
                 };
```

For calendars, there will be an enum for supported calendars as
well as a string name for each.

```c++
#define NUM_SUPPORTED_CALENDARS=9

enum CalendarKind {
   CalendarGregorian=1,   ///< usual Gregorian calendar
   CalendarNoLeap,        ///< Gregorian, but without leap yrs
   CalendarJulian,        ///< Julian
   CalendarJulianDay,     ///< Julian day
   CalendarModJulianDay,  ///< modified Julian day
   Calendar360Day,        ///< 12 months, 30 days each
   CalendarCustom,        ///< user defined
   CalendarNoCalendar,    ///< track elapsed time only
   CalendarUnknown};      ///< uninitialized or invalid

const std::string CalendarKindName[CALENDAR_KIND_COUNT] = {
     "Gregorian",
     "No Leap",
     "Julian",
     "Julian Day",
     "Modified Julian Day",
     "360 Day",
     "Custom",
     "No Calendar",
     "Invalid" };
```

#### 4.1.2 Class/structs/data types

There are six classes that will make up the time interval.

##### 4.1.2.1 TimeFrac class

There will be a TimeFrac class for the base fractional time
representation:

```c++
class TimeFrac {

   // private variables
   private:
      long long whole;  ///< whole seconds
      long long numer;  ///< fractional second (n/d) numerator
      long long denom;  ///< fractional second (n/d) denominator

   public:
      [methods described below]
};
```

##### 4.1.2.2 Calendar class

The calendar class holds useful information for the calendar to be
used:

```c++
#define MONTHS_PER_YEAR=12

class Calendar {

  // private variables
  private:

     int id;                   ///< unique id for quick checks
     static int numCalendars;  ///< number of calendars created
     std::string name;         ///< name of calendar
     CalendarKind calKind;     ///< enum for calendar kind
     std::string calKindName;  ///< name of calendar kind

     // variables defining calendar characteristics for time
     int daysPerMonth[MONTHS_PER_YEAR]; ///< days in each month
     int monthsPerYear;             ///< num months in year
     int secondsPerDay;             ///< seconds per day
     int secondsPerYear;            ///< seconds per normal year
     int daysPerYear;               ///< days per normal year

     // public methods
     public:
        [methods described below]
};
```

##### 4.1.2.3 TimeInstant class

The time instant class represents a point in time within a
given calendar:

```c++
class TimeInstant {

   // private variables
   private:
      TimeFrac elapsedTime; ///< Fractional seconds since reference time
      Calendar *calPtr;  ///< Pointer to calendar in which time is based

   public:
      [methods described below]
};
```


##### 4.1.2.4 TimeInterval class

The time interval is used for both time steps and to compute
differences between time instants. The time interval can either
be an time interval in seconds or it can be used to hold a
time interval in calendar units (eg number of days, months or
years). The latter is useful for periodic events that occur
once per year, month or ndays.

```c++
class TimeInterval {

   // private variables
   private:
      TimeFrac interval; ///< Non-calendar interval in fractional seconds
      bool isCalendar;       ///< True if calendar interval
      long long calInterval; ///< Calendar interval length
      TimeUnits units;       ///< Calendar interval units

   public:
      [methods described below]
};
```


##### 4.1.2.5 Alarm class

The alarm class allows a user to set either one-time or periodic
alarms to trigger events, like forcing updates, I/O, etc. that
occur at specific times.

```c++
class Alarm {

   // private variables
   private:
      std::string name; ///< name for the alarm

      bool ringing;   ///< alarm is currently ringing
      bool periodic;  ///< alarm rings periodically on interval
      bool stopped;   ///< alarm has been stopped and not reset

      TimeInstant  ringTime;     ///< time at/after which alarm rings
      TimeInterval ringInterval; ///< interval at which this alarm rings
      TimeInstant  ringTimePrev; ///< previous alarm time for interval alarms

   public:
      [methods described below]
};
```

##### 4.1.2.6 Clock class

The clock class is meant to track and manage the time for a model
advancing in time. Alarms can be attached to the clock so that
the ringing status can be updated as the clock marches forward.

```c++
class Clock {

   // private variables
   private:

      TimeInstant  startTime; ///< initial time for this clock
      TimeInstant  currTime;  ///< current time
      TimeInstant  prevTime;  ///< time at previous timestep
      TimeInstant  nextTime;  ///< time at next timestep
      TimeInterval timeStep;  ///< interval at which this clock advances

      int numAlarms; ///< current number of attached alarms

      std::vector<Alarm *> alarms; ///< pointers to alarms associated with this clock

   public:
      [methods described below]
};
```

### 4.2 Methods

The classes above have a number of methods associated with them for
managing time.

#### 4.2.1 TimeFrac class

The TimeFrac class is a base class that is not meant to be
accessed by users but provides the integer fraction time
capability needed for accumulating time. It contains a number
of operators to perform arithmetic on fractional integers
and a number of accessor functions to convert units into
a fractional time.

```c++
// Accessor methods

/// Single call to set all native base time components
/// \return error code
int set(const long long whole, ///< [in] whole seconds
        const long long numer, ///< [in] fractional second numerator
        const long long denom  ///< [in] fractional second denominator
        );
/// Set base time by converting from integer hours, minutes, seconds
/// \return error code
int setHMS(const int hours,   ///< [in] integer hours
           const int minutes, ///< [in] integer minutes
           const int seconds  ///< [in] integer seconds

/// Set base time by converting from a real number of seconds
/// \return error code
int setSeconds(const double seconds ///< [in] Time in real seconds
              );
/// Set base time by converting from a real number of hours
/// \return error code
int setHours(const double hours ///< [in] Time in real hours
            );
/// Set base time by converting from a real number of minutes
/// \return error code
int setMinutes(const double minutes ///< [in] Time in real minutes
              );

/// Set whole seconds separately
/// \return error code
int setWhole(
       const long long whole ///< [in] Whole number of seconds
       );
/// Set numerator of fractional seconds separately
/// \return error code
int setNumer(
       const long long numer ///< [in] Numerator of fractional seconds
       );
/// Set denominator of fractional seconds separately
/// \return error code
int setDenom(
       const long long denom ///< [in] Denominator of fractional seconds
       );

/// Single call to retrieve native base time components
/// \return error code
int get(long long &whole, ///< [out] whole seconds
        long long &numer, ///< [out] fractional second numerator
        long long &denom  ///< [out] fractional second denominator
        ) const;
/// Get base time converted to integer hours, minutes, seconds
/// \return error code
int getHMS(int &hours,   ///< [out] integer hours
           int &minutes, ///< [out] integer minutes
           int &seconds  ///< [out] integer seconds
           ) const;
/// Get base time and convert to a real number of seconds
/// \return Time in real seconds
double getSeconds(void) const;
/// Get base time and convert to a real number of hours
/// \return Time in real hours
double getHours(void) const;
/// Get base time and convert to a real number of minutes
/// \return Time in real minutes
double getMinutes(void) const;
/// Retrieve the whole seconds component of base time
/// \return Whole number of seconds
long long getWhole(void) const;
/// Retrieve the numerator component of fractional base time
/// \return Numerator of fractional seconds
long long getNumer(void) const;
/// Retrieve the denominator component of fractional base time
/// \return Denominator of fractional seconds
long long getDenom(void) const;

// constructors/destructors
/// Default base time constructor
TimeFrac(void);
/// Copy constructor for base time
TimeFrac(const TimeFrac& ///< [in] existing base time to be copied
        );
/// Construct base time by component
TimeFrac(
    const long long whole, ///< [in] whole seconds
    const long long numer, ///< [in] fractional second numerator
    const long long denom  ///< [in] fractional second denominator
    );
/// Construct base time by converting from a real number of seconds
TimeFrac(const double seconds ///< [in] Time in real seconds
        );
/// Destructor for base time
~TimeFrac(void);

// operators
/// Equivalence comparison operator for TimeFrac
bool operator==(const TimeFrac &) const;
/// Non-equivalence comparison operator for TimeFrac
bool operator!=(const TimeFrac &) const;
/// Less than comparison operator for TimeFrac
bool operator< (const TimeFrac &) const;
/// Greater than comparison operator for TimeFrac
bool operator> (const TimeFrac &) const;
/// Less than or equal comparison operator for TimeFrac
bool operator<=(const TimeFrac &) const;
/// Greater than or equal comparison operator for TimeFrac
bool operator>=(const TimeFrac &) const;
/// Addition operator for TimeFrac
TimeFrac  operator+ (const TimeFrac &) const;
/// Subtraction operator for TimeFrac
TimeFrac  operator- (const TimeFrac &) const;
/// Increment operator for TimeFrac
TimeFrac& operator+=(const TimeFrac &);
/// Decrement operator for TimeFrac
TimeFrac& operator-=(const TimeFrac &);
/// Multiplication by integer scalar
TimeFrac  operator* (const int multiplier) const;
/// Multiplication in place by integer scalar
TimeFrac& operator*=(const int multiplier);
/// Multiplication by real scalar
TimeFrac  operator* (const double multiplier) const;
/// Multiplication in place by real scalar
TimeFrac& operator*=(const double multiplier);
/// Divide TimeFrac by integer scalar
TimeFrac  operator/ (const int divisor) const;
/// Divide TimeFrac in place by integer scalar
TimeFrac& operator/=(const int divisor);
/// Divide two TimeFracs and return a real result
double   operator/ (const TimeFrac &) const;
/// Modulus method for TimeFrac
TimeFrac  operator% (const TimeFrac &) const;
/// Modulus method in place
TimeFrac& operator%=(const TimeFrac &);
/// Assignment operator for TimeFrac
TimeFrac& operator=(const TimeFrac &);

// Other utility methods
/// Convert a time fraction to new denominator
/// \return error code
int convert(const long long denom ///< [in] new denominator
           );
/// Reduce a time fraction to simplest form
int simplify(void);
```


#### 4.2.2 Calendar class

The calendar class is a mostly immutable class that
holds information about the chosen calendar for use by other
time manager classes. It mostly constructs an instance based
on the user-selected calendar. Retrieval functions and query
functions are provided. An equivalence/non-equivalence operator
is needed for the later TimeInstant equivalence. Utility functions
to convert between elapsed time and calendar dates and incrementing
calendar time are supplied for use by other time manager routines.

```c++
// accessor functions
// this is (mostly) an immutable class so use constructors
// and provide only one set accessor for renaming

// the only set function is for renaming
/// Renames a Calendar to the input string
/// \return Error code
int rename(const std::string inName ///< [in] name to use for calendar
           );

/// Retrieve any/all calendar properties
/// \return Error code
int get(int           *outId,   ///< [out] id assigned to calendar
        std::string   *outName, ///< [out] Name of calendar
        CalendarKind  *outKind, ///< [out] Kind of calendar
        int           *outDaysPerMonth,  ///< [out] Days per month
        int           *outMonthsPerYear, ///< [out] Months per year
        int           *outSecondsPerDay, ///< [out] Seconds per day
        int           *outSecondsPerYear,///< [out] Seconds per year
        int           *outDaysPerYear    ///< [out] Days per year (DPY)
        ) const;

// Might also add specific retrievals for individual
// components, eg seconds per day/year?

/// Default constructor
Calendar(void);
/// Copy constructor
Calendar(const Calendar &calendar);
/// Constructor based on kind of calendar
Calendar(std::string inName,       ///< [in] name of calendar
         CalendarKind calKind      ///< [in] choice of calendar kind
         );
/// Constructs custom calendar based in inputs
Calendar(const std::string inName, ///< [in] name of calendar
         int *inDaysPerMonth,      ///< [in] array of days per month
         int inSecondsPerDay,      ///< [in] seconds per day
         int inSecondsPerYear,     ///< [in] seconds per year
         int inDaysPerYear         ///< [in] days per year (dpy)
         );
/// Calendar destructor
~Calendar(void);

/// Calendar equivalence operator
bool operator==(const Calendar &calendar) const;
/// Calendar non-equivalence operator
bool operator!=(const Calendar &calendar) const;

/// Checks whether input year is a leap year
/// \return true if year is a leap year, false otherwise
bool isLeapYear(long long year, ///< [in]  year to check
                int &rc         ///< [out] return code to flag errors
                ) const;

/// Computes the total elapsed time in seconds (in TimeFrac form)
/// since the calendar reference time, given a calendar date, time.
/// \return Elapsed time in TimeFrac form
TimeFrac getElapsedTime(
    const long long year,   ///< [in] calendar year
    const long long month,  ///< [in] calendar month
    const long long day,    ///< [in] calendar day
    const long long hour,   ///< [in] time of day-hour
    const long long minute, ///< [in] time of day-min
    const long long whole,  ///< [in] time of day-whole seconds
    const long long numer,  ///< [in] time of day-frac secs (numerator)
    const long long denom   ///< [in] time of day-frac secs (denom)
    ) const;

/// Determines the calendar date and time of day, given an
/// elapsed time since the calendar reference time.
/// \return error code
int getDateTime(
    const TimeFrac elapsedTime, ///< [in] time in secs from ref time
    long long &year,   ///< [out] calendar year
    long long &month,  ///< [out] calendar month
    long long &day,    ///< [out] calendar day
    long long &hour,   ///< [out] time of day-hours
    long long &minute, ///< [out] time of day-minutes
    long long &whole,  ///< [out] time of day-whole seconds
    long long &numer,  ///< [out] time of day-frac secs (numerator)
    long long &denom   ///< [out] time of day-frac secs (denom)
    ) const;

/// Increments (or decrements) a calendar date by a specified
/// interval, supplied by an integer interval in given time units.
/// Only calendar based intervals (years, months or days) are
/// supported. This is primarily meant to be called by other time
/// manager routines (eg to add/subtract time instants) for those
/// time intervals that are dependent on date and sensitive to
/// calendar features like leap years and varying days of the month.
/// \return error code
int incrementDate(
    const long long interval, ///< [in] time interval to advance date
    const TimeUnits units,    ///< [in] time units for interval
    long long &year,   ///< [in,out] calendar year of time to be changed
    long long &month,  ///< [in,out] calendar month of ...
    long long &day     ///< [in,out] calendar day
    ) const;
```

#### 4.2.3 TimeInstant class

The time instant class represents a single instant in time.
Most of the methods are accessor methods for getting/setting
a time instant or its components. In addition, a number of
operators are defined for performing basic arithmetic with
time instants and time intervals (see following interval class).
Finally, a method for creating a time string for a time
instant is supplied.

```c++
// constructors/destructors
/// Default constructor creates empty time instant
TimeInstant(void);

/// Construct time instant from date, time, calendar
/// Where seconds is supplied as real number.
TimeInstant(Calendar       *Cal,    ///< [in] Calendar to use
            const long long year,   ///< [in] year
            const long long month,  ///< [in] month
            const long long day,    ///< [in] day
            const long long hour,   ///< [in] hour
            const long long minute, ///< [in] minute
            const double    rSecond ///< [in] second (real)
            );

/// Construct time instant from date, time, calendar
/// Where seconds is supplied in integer fractional seconds.
TimeInstant(      Calendar *Cal,    ///< [in] Calendar to use
            const long long year,   ///< [in] year
            const long long month,  ///< [in] month
            const long long day,    ///< [in] day
            const long long hour,   ///< [in] hour
            const long long minute, ///< [in] minute
            const long long whole,  ///< [in] second (whole integer)
            const long long numer,  ///< [in] second (fraction numerator)
            const long long denom   ///< [in] second (fraction denominator)
            );

/// Destructor for time interval
~TimeInstant(void);

// Accessor methods
/// Set time instant calendar
/// \return error code
int set(Calendar *Cal   ///< [in] Calendar to use for this time
        );

/// Set time instant from date and time, where seconds is supplied
/// as a real number.
/// \return error code
int set(const long long year,   ///< [in] year
        const long long month,  ///< [in] month
        const long long day,    ///< [in] day
        const long long hour,   ///< [in] hour
        const long long minute, ///< [in] minute
        const double    rSecond ///< [in] second (real)
        );

/// Set time instant from date and time, where seconds is supplied
/// in integer fractional seconds.
/// \return error code
int set(const long long year,   ///< [in] year
        const long long month,  ///< [in] month
        const long long day,    ///< [in] day
        const long long hour,   ///< [in] hour
        const long long minute, ///< [in] minute
        const long long whole,  ///< [in] second (whole integer)
        const long long numer,  ///< [in] second (fraction numerator)
        const long long denom   ///< [in] second (fraction denominator)
        );

/// Retrieve calendar from time instant
/// \return error code
int get(Calendar *&cal ///< [out] Calendar ptr in which instant defined
        ) const;

/// Retrieve time in date, time form with real seconds.
/// \return error code
int get(long long &year,   ///< [out] year   of this time instant
        long long &month,  ///< [out] month  of this time instant
        long long &day,    ///< [out] day    of this time instant
        long long &hour,   ///< [out] hour   of this time instant
        long long &minute, ///< [out] minute of this time instant
        double    &second  ///< [out] second of this time instant
        ) const;

/// Retrieve time in date, time form with fractional integer seconds.
/// \return error code
int get(long long &year,   ///< [out] year   of this time instant
        long long &month,  ///< [out] month  of this time instant
        long long &day,    ///< [out] day    of this time instant
        long long &hour,   ///< [out] hour   of this time instant
        long long &minute, ///< [out] minute of this time instant
        long long &whole,  ///< [out] whole seconds of this time
        long long &numer,  ///< [out] frac second numerator
        long long &denom   ///< [out] frac second denominator
        ) const;

// Operators on time instants
/// Equivalence comparison for TimeInstant
bool operator==(const TimeInstant &) const;
/// Non-equivalence comparison operator for TimeInstant
bool operator!=(const TimeInstant &) const;
/// Less than comparison operator for TimeInstant
bool operator< (const TimeInstant &) const;
/// Greater than comparison operator for TimeInstant
bool operator> (const TimeInstant &) const;
/// Less than or equal comparison operator for TimeInstant
bool operator<=(const TimeInstant &) const;
/// Greater than or equal comparison operator for TimeInstant
bool operator>=(const TimeInstant &) const;
/// Increment time by adding a time interval
TimeInstant  operator+ (const TimeInterval &) const;
/// Decrement time by subtracting a time interval
TimeInstant  operator- (const TimeInterval &) const;
/// Create a time interval by subtracting two time instants
TimeInterval operator- (const TimeInstant &) const;
/// Increment time in place by adding time interval
/// Increment time in place by adding time interval
TimeInstant& operator+=(const TimeInterval &);
/// Decrement time in place by subtracting time interval
TimeInstant& operator-=(const TimeInterval &);

// Other utility methods
/// Get time as a string in the format
/// 'YYYYYY-MM-DD{separator}HH:MM:SS.SSSSSS' where the number
/// of digits in year, number of digits after the decimal in
/// seconds and the character to use as separator are all input
/// by the user.
/// \return time string
std::string getString(
    const int yearWidth,   ///< [in] number of digits in year
    const int secondWidth, ///< [in] num digits after decimal in seconds
    std::string separator  ///< [in] string(char) to separate date/time
    ) const;
```

#### 4.2.4 TimeInterval class

The time interval class manages differences or increments in time.
The user can create or set a time interval with a value and units.
Retrieval functions for the time interval in various forms are
also supplied. A number of operators are defined for a variety
of mathematical operations on time intervals (add, subtract,
multiply, equivalence, absolute value). Finally, the time instant
class is treated as a friend class so that time instants can
be incremented by a time interval and time intervals can be created
by subtracting two time instants.

```c++
// constructors/destructors
/// Default time interval constructor
TimeInterval(void);

/// Construct time interval from base time fractional integer seconds
TimeInterval(const long long whole, ///< Whole integer seconds
             const long long numer, ///< Fractional seconds numerator
             const long long denom  ///< Fractional seconds denominator
             );

/// Construct time interval from an integer length and units
TimeInterval(const int       length, ///< length of time interval
             const TimeUnits units   ///< unit of time for interval
             );

/// Construct time interval from a long long integer length and units
TimeInterval(const long long length, ///< length of time interval
             const TimeUnits units   ///< unit of time for interval
             );

/// Construct time interval from an real length and units
TimeInterval(const double length,    ///< length of time interval
             const TimeUnits units   ///< unit of time for interval
             );

/// Destructor for time interval
~TimeInterval(void);

// Accessor methods
/// Set a non-calendar interval in native fractional integer seconds
/// \return error code
int set(const long long whole, ///< Whole integer seconds
        const long long numer, ///< Fractional seconds numerator
        const long long denom  ///< Fractional seconds denominator
        );

/// Set a time interval from integer length and units
/// \return error code
int set(const int       length, ///< length of time interval
        const TimeUnits units   ///< unit of time for interval
        );

/// Set a time interval from long long integer length and units
/// \return error code
int set(const long long length, ///< length of time interval
        const TimeUnits units   ///< unit of time for interval
        );

/// Set a time interval from an real length and units
/// Real length is only supported for non-calendar intervals since
/// a non-integral length has ambiguous meaning when months, years
/// have variable length.
/// \return error code
int set(const double length,    ///< length of time interval
        const TimeUnits units   ///< unit of time for interval
        );

/// Retrieve non-calendar interval in native fractional integer form
/// \return error code
int get(long long &whole, ///< [out] whole seconds
        long long &numer, ///< [out] fractional second numerator
        long long &denom  ///< [out] fractional second denominator
        ) const;

/// Retrieve a time interval in integer length in specified units.
/// To avoid roundoff issues during conversions, integer retrieval
/// is only permitted in the same units in which the interval was
/// defined.
/// \return error code
int get(long long &length, ///< [out] requested integer length of interval
        const TimeUnits units ///< [in] unit of time for interval
        ) const;

/// Retrieve a time interval in real length and specified units.
/// For calendar intervals, the units must match the units in which
/// the interval was defined. For non-calendar intervals, only conversions
/// between hours, minutes and seconds are allowed.
/// \return error code
int get(double &length, ///< [out] Requested time interval length
        const TimeUnits units ///< [in] unit of time for interval
       ) const;

/// Equivalence comparison operator for TimeInterval
bool operator==(const TimeInterval &) const;
/// Non-equivalence comparison operator for TimeInterval
bool operator!=(const TimeInterval &) const;
/// Less than comparison operator for TimeInterval
bool operator< (const TimeInterval &) const;
/// Greater than comparison operator for TimeInterval
bool operator> (const TimeInterval &) const;
/// Less than or equal comparison operator for TimeInterval
bool operator<=(const TimeInterval &) const;
/// Greater than or equal comparison operator for TimeInterval
bool operator>=(const TimeInterval &) const;
/// Addition operator for TimeInterval
TimeInterval  operator+ (const TimeInterval &) const;
/// Subtraction operator for TimeInterval
TimeInterval  operator- (const TimeInterval &) const;
/// Increment operator for TimeInterval
TimeInterval& operator+=(const TimeInterval &);
/// Decrement operator for TimeInterval
TimeInterval& operator-=(const TimeInterval &);
/// Multiplication by integer scalar
TimeInterval  operator* (const int multiplier) const;
/// Multiplication by real scalar
/// This is primarily meant for non-calendar intervals, but
/// for calendar intervals, it will multiply the year, month or
/// day interval and convert to the nearest integer.
TimeInterval  operator* (const double multiplier) const;
/// Multiplication in place by integer scalar
TimeInterval& operator*=(const int multiplier);
/// Multiplication in place by real scalar
/// This is primarily meant for non-calendar intervals, but
/// for calendar intervals, it will multiply the year, month or
/// day interval and convert to the nearest integer.
TimeInterval& operator*=(const double multiplier);
/// Divide interval by integer scalar
TimeInterval  operator/ (const int divisor) const;
/// Divide interval in place by integer scalar
TimeInterval& operator/=(const int divisor);

/// Absolute value
static TimeInterval absValue(const TimeInterval &);
/// Negative absolute value
static TimeInterval negAbsValue(const TimeInterval &);

// Other utility methods
/// Check whether a time interval is positive
bool isPositive(void);

/// Give the Time Instant access to Time Interval so
/// that incrementing/decrementing time is easier.
friend class TimeInstant;

}
```

#### 4.2.5 Alarm Class

The alarm class allows the creation of both one-time and
periodic alarms. There are also methods for querying the
state of an alarm and for resetting.

```c++
// constructors/destructors
/// Constructs a one-time alarm using the input ring time.
Alarm(const std::string inName,   ///< [in] Name of alarm
      const TimeInstant alarmTime ///< [in] Time at/after which alarm rings
     );

/// Constructs a periodic/interval alarm based on an input periodic
/// time interval and a start time for the interval period.
Alarm(const std::string  inName,        ///< [in] Name of alarm
      const TimeInterval alarmInterval, ///< [in] interval at which alarm rings
      const TimeInstant  intervalStart  ///< [in] start time of first interval
     );

/// Destructor for alarm
~Alarm(void);

/// Check whether an alarm is ringing
/// \return true if alarm is ringing, false otherwise
bool isRinging(void);

/// Checks whether the alarm should ring based on the current
/// (or supplied) time instant (returns error code)
int updateStatus(const TimeInstant currentTime ///< [in] current time
                );

/// Stops a ringing alarm and sets next a new ring time. If the alarm
/// is a periodic/interval alarm, the next ring time is set to be the
/// next interval boundary after the input time.  If the alarm is a
/// single instance, the input time is used as the next alarm time.
/// (returns error code)
int reset(const TimeInstant inTime ///< [in] new basis for alarm time
         );

/// Stops a ringing alarm (returns error code).
int stop(void);

/// Rename an alarm (returns error code)
int rename(const std::string newName ///< [in] new name for alarm
          );

/// Get alarm name
std::string getName(void) const;
```

#### 4.2.6 Clock Class

The Clock class defines the model clock and is created based
on a start time and time step. There are a number of retrieval
functions to get the current, past, and next times. A number
of alarms can be attached to a clock. Finally, there is
a method to advance the clock one time and update the state
of all attached alarms.

```c++
// constructors/destructors
/// Construct a clock from start time and time step
Clock(const TimeInstant  startTime, ///< [in] Start time for clock
      const TimeInterval timeStep   ///< [in] Time step to advance clock
     );

/// Destructor for clocks
~Clock(void);

// Accessor Methods
/// Set the current time (returns error code)
int setCurrentTime(
    const TimeInstant currTime ///< [in] new value for current time
    );

/// Changes the time step for this clock (returns error code)
int changeTimeStep(
    const TimeInterval newTimeStep ///< [in] new value for time step
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
int attachAlarm(Alarm* inAlarm ///< [in] pointer to alarm to attach
               );

/// Advance a clock one timestep and update status of any attached
/// alarms. (returns error code)
int advance(void);
```


## 5 Verification and Testing

### 5.1 Test TimeFrac

We will test the base TimeFrac class by creating some reference
times, one of which will be a fractional time that cannot be
exactly represented by a float (eg 1/3 seconds). Constructors and
set functions will be paired with get functions to make sure both
set/get are working. All operators (equivalence, relational,
algebraic) operators will be tested with known values. We will
try to anticipate values that could fail.
  - tests parts of requirement 2.2, 2.5, 2.6

### 5.2 Test TimeInstant

Similar to TimeFrac, create reference times in each of the
supported calendars using constructor and set functions. Test
retrievals (get) to verify they return the reference value.
Test all the supported TimeInstant operators with known values
in each of the calendars.
  - tests parts of requirement 2.1, 2.3

### 5.3 Test Calendar

For each of the supported calendars (including a custom
calendar), create an instance of the calendar and use
the retrieval (get) functions to verify that various
calendar properties are as expected (eg days per year,
days per month, leap years/days, etc.). For some calendars,
known reference dates exist from the references above and
we will verify that these result in correct elapsed time
or Julian day. Test the few operators (equivalence, copy
constructors).
  - tests most of requirement 2.3

### 5.4 Test Time Intervals

As in the others, test constructors and set functions by
pairing with retrieval (get) functions for a given
reference time interval. Test assignment constructors
by creating a time interval for a given time step in
various time units (seconds, hours, days, etc.). Test
equivalence/non-equivalence by comparing known time intervals.
Test algebraic operators (addition, subtraction, increment,
decrement, multiply, divide, etc.) using reference values.
Test creation of time interval by subtraction of time instants.
  - tests requirement 2.4, parts of 2.1

### 5.5 Test Alarms

Test both one-time and periodic alarms at various supported
intervals, including annual, monthly/nmonths, daily/ndays,
hourly/nhours, nminutes, seconds. This is accomplished by
choosing a start time in a calendar and integrating over a
time period that covers the alarm. Use the Gregorian calendar
(as most complex) and choose an interval that includes leap
years/days. Also verify that alarms are not ringing at arbitrary
selected times in between. Check the reset function and verify
alarms are reset appropriately.
  - tests requirement 2.5

### 5.6 Test Clock

This is the primary test of putting it all together and
verifying that we can advance a model clock. Construct a
clock with a given calendar and time step. Test various
retrieval functions to get known quantities - also test
the ability to change time step with set/get pair. Set a
number of alarms at various times and periodic intervals.
Integrate forward in time for N years and verify all alarms
are triggered at proper times. Also verify that the final
time is exact with no roundoff error accumulated. Create
another clock with large time steps (years) that extend over a
potential paleoclimate interval to test large number
representations are supported correctly.
  - tests requirement 2.1, 2.2, 2.3, 2.4, 2.5, 2.6
