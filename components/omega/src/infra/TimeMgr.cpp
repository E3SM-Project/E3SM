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

} // namespace OMEGA
//===-----------------------------------------------------------------------===/
