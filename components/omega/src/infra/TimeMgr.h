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

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_TIMEMGR_H
