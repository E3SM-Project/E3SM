(omega-dev-time-manager)=

# Time Manager (TimeMgr)

The Omega Time Manager module tracks simulation time during a run and manages
alarms for triggering events at specified times This is handled by six
interacting defined classes, each of which have a number of methods defined for
accessing, comparing, and manipulating the class members.

## Class Descriptions

### 1. TimeFrac

The TimeFrac class is the fundamental base time representation. To avoid
potential rounding errors that can accumulate over millions of time steps, time
is tracked in seconds as a fractional representation using three 8-byte integers
representing the whole number, numerator, and denominator. Most of the other
classes in the Time Manager are based on the TimeFrac class and generally
convert time representations from other units and store them as a TimeFrac
object.

### 2. Calendar

The Calendar class is mostly an immutable class that stores all information for
supported calendars that can be utilized in a simulation. The supported
calendars are: Gregorian, No Leap, Julian, Julian Day, Modified Julian Day,
360 Day, Custom, and No Calendar. These supported calendars are specified via an
enum defined in the Time Manager header file:
```c++
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
```
During a simulation, only one calendar can be defined and is set with the call:
```c++
Calendar::init("CalendarName");
```
that must be set early in the model initialization before other time quantities
are defined. CalendarName is a string associated with the defined calendars
above and found in the list CalendarKindName in ``TimeMgr.h``. Typically, they
are the same as the kinds above with a space for multi-word names
(eg "No Leap"). If a custom calendar is desired rather than the supported
calendars, a custom calendar can be defined with a different init
interface in which all the relevant quantities are provided:
```c++
Calendar::init(
   std::vector<I4> &InDaysPerMonth, ///< [in] array of days per month
   I4 InSecondsPerDay,              ///< [in] seconds per day
   I4 InSecondsPerYear,             ///< [in] seconds per year
   I4 InDaysPerYear                 ///< [in] days per year (dpy)
);
```
Time is generally tracked as the amount of time elapsed from a reference date
for a given calendar. A number of public functions exist to retrieve Calendar
information, convert between dates and elapsed time, check for leap years and
other checks. However, these are typically used by other TimeMgr classes and
not directly by the developer. Examples of their use can be found in the
TimeMgr unit test.

### 3. TimeInstant

The TimeInstant class represents a moment in time with the particular calendar
that has been defined for the simulation. It consists of a TimeFrac that
represents the time since a reference time associated with the defined calendar.
There are three constructors for initializing a TimeInstant. One method is with
five 8-byte integers, and an 8-byte real value:
```c++
OMEGA::TimeInstant TI1(Year, Month, Day, Hour, Minute, RealSecond);
```
Alternatively, the TimeInstant can be initialized with eight 8-byte integers,
with the seconds value represented by three 8-byte integers:
```c++
OMEGA::TimeInstant TI2(Y, M, D, H, M, Whole, Numer, Denom);
```
A final constructor creates a time instant based on a time string that
conforms roughly to the ISO standard `"YYYYYY-MM-DD_HH:MM:SS.SSSS"` though
the constructor allows for any single-character non-numeric separator between
each of the numeric fields and the width of the YY and SS fields can be up
to the 8-byte standards for integer and floats, respectively.
```c++
OMEGA::TimeInstant TI3(TimeString);
```

Among the methods defined in the TimeInstant class is `getString` which will
produce a `std::string` representation of the time conforming to ISO standards,
 e.g. `"YYYYYY-MM-DD HH:MM:SS.SSSSS"`. The `getString` method takes two 4-byte
integers which represent the number of digits to use for the year and the
number of digits for the decimal portion of seconds, as well as a `std::string`
for the seprator between the Y/M/D calendar date and the H/M/S clock time. A
string formatted as above is returned with the following call:
```c++
TI1.getString(6, 5, " ");
```
A range of accessor functions can also get/set the TimeInstant in various
forms.

### 4. TimeInterval

The TimeInterval class represents a length of time between two moments in time.
A TimeInterval can represent the time step, but is also utilized by other
classes for computing differences between TimeInstant objects, such as when
checking to trigger an Alarm. A TimeInterval can be initialized with a
fractional second representation, by supplying three 8-byte integers:
```c++
OMEGA::TimeInterval Interval1(Whole, Numerator, Denominator);
```
Alternatively, a TimeInterval can also be initialized with either an integer
or floating point value and a specified unit of time:
```c++
OMEGA::TimeInterval Interval2(Length, UnitType);
```
The supported time units are defined via an enum class in the Time Manager
header file:
```c++
enum class TimeUnits {
   None = 0, ///< value for undefined units
   Seconds,  ///< time units in seconds (typical)
   Minutes,  ///< time units in minutes
   Hours,    ///< time units in hours
   Days,     ///< time units in days
   Months,   ///< time units in months
   Years,    ///< time units in years
};
```
Finally, a time interval can be defined as the time between two time instants:
```c++
OMEGA::TimeInterval MyDeltaTime = MyTimeInstant2 - MyTimeInstant1;
```

### 5. Alarm

The Alarm class is designed to trigger events at specified times. Alarms can be
set to ring at a single time or to ring periodically. A one-time Alarm is
initialized with a `std::string` name and a TimeInstant:
```c++
OMEGA::Alarm SingleAlarm("Single Alarm", AlarmInstant);
```
A periodic Alarm is initialized with a `std::string` name, a TimeInterval, and a
TimeInstant representing the start of the first interval:
```c++
OMEGA::Alarm PeriodicAlarm("Periodic Alarm", AlarmInterval, AlarmStart);
```
The status of an Alarm can be checked using the `isRinging` method, which
returns true or false depending on the Alarm status:
```c++
SingleAlarm.isRinging();
```
The `reset` method takes a TimeInstant as input and will switch a ringing Alarm
off, as well as set a new ring time. If it is a one-time Alarm, the input time
is set as the new ring time:
```c++
SingleAlarm.reset(NewAlarmTime);
```
If the Alarm is periodic, the ringing will be switched off as above and the
new ring time will be set to the next interval boundary after the input time.
The interval boundary is an integer number of intervals after the start time
provided on Alarm creation.
```c++
PeriodicAlarm.reset(CurrentTime);
```
An Alarm can be permanently stopped using the `stop` method:
```c++
SingleAlarm.stop();
```
It is sometimes useful to retrieve other aspects of a periodic Alarm, namely
the Alarm interval and the last time the Alarm was ringing. Two retrieval
functions are provided that return const pointers to these values:
```c++
const TimeInterval *AlarmInterval = PeriodicAlarm.getInterval();
const TimeInstant *PriorRingTime = PeriodicAlarm.getRingTimePrev();
```

### 6. Clock

The Clock class is the primary object for managing time during a simulation.
A Clock is initialized with a TimeInstant representing the start time and a
TimeInterval representing the time step the Clock advances with each integration
step:
```c++
OMEGA::Clock ModelClock(StartTime, TimeStep);
```
The class contains a `std::vector` of pointers to Alarms. Any number of Alarms
can be attached to the Clock using the `attachAlarm` method:
```c++
ModelClock.attachAlarm(&SingleAlarm);
```
The Clock is advanced during each time step with the `advance` method:
```c++
ModelClock.advance();
```
The Clock will automatically trigger any attached Alarms as the Clock marches
forward. The time step for a Clock can be changed by passing a TimeInterval
to the `changeTimeStep` method:
```c++
ModelClock.changeTimeStep(NewTimeStep);
```
During a simulation, the `getCurrentTime` method will return the current time
stored in the Clock as a TimeInstant:
```c++
OMEGA::TimeInstant CurrentTime = ModelClock.getCurrentTime();
```

## Implementation

In order to utilize the Time Manager, the first step is to initialize a
Calendar object:
```c++
OMEGA::Calendar CalGreg("Gregorian", OMEGA::CalendarGregorian);
```
A start time and a time step are needed in order to initialize the Clock. The
following lines will initialize a Clock that starts Jan 1st, 2000 at midnight
with a 20 minute time step:
```c++
OMEGA::TimeInstant StartTime(&CalGreg, 2000, 1, 1, 0, 0, 0.0);
OMEGA::TimeInterval TimeStep(20, OMEGA::TimeUnits::Minutes);

OMEGA::Clock ModelClock(StartTime, TimeStep);
```
Now, Alarms can be attached to the Clock. For example, if some output
is desired on the first of each month, the following lines will attach a
monthly Alarm:
```c++
OMEGA::TimeInterval IntervalMonthly(1, OMEGA::TimeUnits::Months);
OMEGA::Alarm AlarmMonthly("Every Month", IntervalMontly, StartTime);

Err = ModelClock.attachAlarm(AlarmMonthly);
```
The Clock can be used to control the duration of a simulation by using a
TimeInstant that represents the stop time:
```c++
OMEGA::TimeInstant StopTime(&CalGreg, 2025, 1, 1, 0, 0, 0.0);
OMEGA::TimeInstant CurrentTime = StartTime;

while (CurrentTime < StopTime) {

   Err = ModelClock.advance();
   CurrentTime = ModelClock.getCurrentTime();

  ...

}
```
To utilize any attached Alarms, the ringing status needs to be checked each
time step:
```c++
if (AlarmMonthly.isRinging()) {
   AlarmMonthly.reset(CurrentTime);

   ...
}
```
