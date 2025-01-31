(omega-design-timers)=
# Timers

## 1 Overview

Within OMEGA, we will need a method for profiling aspects
of code performance. In the most simple form, this requires
the ability to report the time spent in various sections of
code. It may also be desirable to track other aspects of
performance or to access hardware timers, though there may be
external tools more appropriate for more detailed profiling.

## 2 Requirements

### 2.1 Requirement: Timers

We require the ability to track the time spent in a selected
code section. This requires the ability to start, stop and
accumulate time in a timer with a given name.

### 2.2 Requirement: Timing statistics

For tracking load balance, we require the timers to report
minimum, maximum and avg time across MPI tasks in a parallel
simulation.

### 2.3 Requirement: Timing context

Related to the statistics requirement, timers may be called
from a code region that is operating within a subset of the
default MPI environment so must track the context in which
they are called.

### 2.4 Requirement: Call sequence

Because some utility functions and related timers may be
called from multiple sections of code, we require the timers
to track where in the call sequence the timer occurs and to
report the time separately when called from a different
location in the code. When reporting, it is desirable to
format timer output with indentation reflecting the level of
the call tree where the timer is called.

### 2.5 Requirement: Output and E3SM compatibility

We require timer output to be compatible with existing
E3SM performance tools. Output in other common formats
for performance tools is also desired. Output in
human-readable form is also required and it is often useful
for at least a summary profile in readable form be
included in the Log file output.  For consistency with
other E3SM components, we may wish to support the General
Purpose Timing Library (GPTL) as a timer option.

### 2.6 Requirement: Full coverage

Timers must cover all code within Omega to provide a
complete profile. Granularity should be at least at the
function/routine level, though finer granularity may
be desired for larger routines. It may be desireable
to provide a preprocessor option to support very fine
granularity and more detail if desired.

### 2.7 Requirement: MPI context

Because OMEGA allows portions of code to run on a subset
of processors or alternative parallel decompositions, the
timers must be able to distinguish the context so that
statistics across MPI tasks are computed correctly.

### 2.8 Desired: Support for multiple devices

It is useful to track performance aspects of hybrid
computing nodes, including kernel launch times and
device timers. However, it is likely that that external
profiling tools are a better way to access this information.

### 2.9 Desired: Other hardware counters

It may be desirable in the future to access other hardware
counters in addition to simple timers (eg flops, memory ops,
cache misses, etc.) to provide additional detail in diagnosing
performance bottlenecks. There are interfaces like PAPI that
may provide a portable means of accessing counters when
available, but this profiling may be better performed with
other tools.

### 2.10 Desired: Thread support

We do not anticipate having large on-node CPU threaded
regions (threads are currently loop-level), so we anticipate
all timers will be called outside of threaded regions.
However, such a capability might be needed in the future.

### 2.11 Desired: Memory profiling

It would be nice to have the capability of tracking
memory use (eg high-water mark) to diagnose memory leaks
or overall memory use by the model. Other memory metrics
are also useful but may fall under 2.9 above or be tracked
via external profiling tools.

### 2.12 Desired: Interaction with other profiling tools

We will likely make use of external profiling tools for
additional diagnostics. While some tools use sampling, others
may require inserting code to assist profiling. It may be
beneficial to have generic interfaces that could support
profiling calls and be optionally activated.

## 3 Algorithmic Formulation

There are no specific algorithms associated with this capability.

## 4 Design

The initial implementation will include only a timer class.
We will support both a reference implementation using C++
stdlib timers as well as a GPTL (General Purpose Timing Library)
option that can be selected at compile time.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

At compile time, the user can optionally select the GPTL
timers as a pre-processor option (USE_GPTL). In the future,
a TIMER_LEVEL pre-processor option might be provided to
identify the granularity of timers, but this is not currently
implemented.

Currently the only timer option set in the run time configuration
is the filename (or filename template) where the timers will
be printed - a blank entry, None or Log will send the output
to the standard Logging output.  Other options may be added
in the future.
```yaml
Timers:
   TimerFile: MyFilename.$YYYY$MM$DD
```

#### 4.1.2 Class/structs/data types

The timers are implemented with a single timer class. For
the reference implementation, this will contain:
```c++
class Timer

 private:
   /// Timer name
   std::string TimerName;

   /// Context in which this timer is called
   std::shared_ptr<MachEnv> Context;

   /// Flag to determine when the timer is running
   bool IsRunning;

   /// Number of times the timer is started/called
   /// Used for timer statistics
   I8 NumCalls;

   /// Start time for the current timer interval
   R8 StartTime;

   /// Accumulated time for this timer
   R8 AccumTime;

   /// Track child and parent timers
   std::shared_ptr<Timer> ParentTimer;
   std::vector<std::shared_ptr<Timer>> ChildTimers;

   /// Current active timer
   /// Used with the above to keep track of call sequence
   static std::shared_ptr<Timer> CurTimer;

   /// Level in the timer call hierarchy
   I4 CallLevel;

   /// Timers will be output with this filename or the
   /// constructed filename based on this template
   static std::string TimerFile;

 public:
    // see methods described later
```

When using external libraries (eg GPTL), those libraries
typically manage the timers internally and the above class will
be largely empty of data but will retain the interfaces
and some of the static members that remain relevant
(eg the output filename).

### 4.2 Methods

All of the methods below are static functions.
The user will be responsible for inserting any needed or desired
synchronization like MPI barriers or Kokkos fences before
timer calls.

#### 4.2.1 Start

The start method is used to either create a new timer with a
given name or start an existing timer of that name:
```c++
Timer::start(const std::string &TimerName);
```
This assumes the timer is called from the default machine
environment. If the timer is started in a different context
than the default machine environment (eg from a defined subset
of tasks), then the relevant MachEnv name must also be passed:
```c++
Timer::start(const std::string &TimerName,
             const std::string &MachEnvName);
```

#### 4.2.2 Stop

The stop method stops a running timer and accumulates the
time in a running sum.
```c++
Timer::stop(const std::string &TimerName);
```

#### 4.2.3 Print

The timer print method stops all timers, computes timer
statistics and prints the results of all timers either to a
file or to stdout. The output will be formatted based
on the call sequence with child timers indented
from their parents. Statistics will include the number of
calls to the timer, the min, max and mean across MPI tasks,
the percentage of total run time and the percentage of the
parent timers' run time.
```c++
Timer::print();
```

#### 4.2.5 Initialize, Finalize

An initialize function will set any run-time options from
the configuration and start a top-level timer named Total.
If an external timing library is used, it will initialize
that library and its options.

A finalize method will stop all timers (including Total),
call the print method to output the timing data, and
clear all defined timers.
```c++
Timer::initialize();
Timer::finalize();
```

## 5 Verification and Testing

### 5.1 Test all internal timers

A unit test with routines and loops of varying run
times will be used to test multiple timers and timer
statistics across MPI tasks. Visual inspection of the
timer output will be used to determine compliance with
formatting.
  - tests requirements 2.1-2.5

### 5.2 Test all timers with GPTL library

The above test code will be repeated after building
and linking with the GPTL library
  - tests requirements 2.1-2.5
