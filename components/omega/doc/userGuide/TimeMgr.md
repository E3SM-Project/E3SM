(omega-user-time-manager)=

# Time Manager (TimeMgr)

The time manager manages various time quantities during a simulation, including
time steps, time instants, the calendar, clocks and alarms. The user can
specify a few options through the ``TimeIntegration`` portion of the Omega
configuration file. In particular, the user can specify the Calendar used
for tracking simulation time. Supported calendars include:
 - Gregorian: typical western calendar with leap years/days
 - No Leap: a Gregorian calendar without leap days
 - Julian: earlier calendar still used by Eastern Orthodox traditions
 - Julian Day: a calendar that simply tracks the number of days since
   a reference day
 - Modified Julian Day: a modified version of Julian Day that changes the
   reference time
 - 360 Day: An artificial calendar with equal 30-day months
 - Custom: A calendar in which a developer can specify the number of
   days per month and days per year. This currently cannot be configured
   entirely via the configuration and would require additional code to
   specify those variables.
 - No Calendar: This is supplied for test cases in which the calendar
   is meaningless and time is simply tracked from time zero.

A number of other time quantities are specified in the ``TimeIntegration``
section (eg time step, start time) and are described in the
[TimeStepping](#omega-user-time-stepping) section.
