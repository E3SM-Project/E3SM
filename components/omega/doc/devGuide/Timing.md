(omega-dev-timing)=

# Timing

The timing infrastructure builds upon the E3SM Pacer library for wall clock timing.
Pacer integrates, whenever possible, platform-specific marker APIs.

## Initialization

To initialize the Pacer library call
```c++
Pacer::initialize(Comm);
```
where `Comm` is an MPI communicator.

## Writing-out timing data

To write out timing data call
```c++
Pacer::print(FilePrefix, PrintAllRanks);
```
where `FilePrefix` is the prefix for all output files and an optional argument `PrintAllRanks` determines
if all MPI ranks should print their data.

## Finalization

To finalize the Pacer library call
```c++
Pacer::finalize();
```

## Basic timer use

To time a region of code enclose it with calls to `Pacer::start` and `Pacer::stop` functions, like so
```c++
Pacer::start(Name, Level);
// region of code to be timed
Pacer::stop(Name, Level);
```
These functions take a string `Name` and a non-negative integer `Level`.
The added timer will be active only if the timing level set in the config file is greater or equal to `Level`.

## Advanced timing functions

### Conditional MPI barriers

Properly timing MPI communication might require inserting MPI barriers. It might be desirable
to remove those barriers in production runs. Pacer provides a function
```c++
Pacer::timingBarrier(TimerName, Level, Comm)
```
which adds an MPI barrier and puts a timer around it using the communicator `Comm`.
Whether barriers added by this function are actually called can be controlled by the following functions
```c++
  Pacer::enableTimingBarriers();
  Pacer::disableTimingBarriers();
```

### Adding parent prefixes

It might be desirable to add a prefix to a group of timers based on their parent timer.
To enable this Pacer provides the `addParentPrefix()` and `removeParentPrefix()` functions.
For example, the following call sequence
```c++
Pacer::start("Parent", 0);
Pacer::addParentPrefix();

Pacer::start("Child", 0);
Pacer::stop("Child", 0);

Pacer::removeParentPrefix();
Pacer::start("Parent", 0);
```
results in output where the "Child" timer shows up as "Parent:Child" in the output files.
This is useful when timers are added inside general purpose routines, that are called
from many places in the code, such as halo exchange.

### Disabling timers

It might be desirable programmatically disable or enable timing. To allow that, Pacer provides
the `disableTimers()` and `enableTimers()` functions. In the following call sequence
```c++
Pacer::disableTiming();

Pacer::start("Timer1", 0);
Pacer::stop("Timer1", 0);

Pacer::enableTiming();

Pacer::start("Timer2", 0);
Pacer::stop("Timer2", 0);
```
`Timer1` is not timed while `Timer2` is. This is useful mainly when done conditionally.
For example, the first call to some function takes much longer than subsequent calls,
and having a detailed timing breakdown of the first call is not important. In that case,
it might be desirable to have a separate timer for the first call with its child timers disabled.
