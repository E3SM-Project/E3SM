(omega-design-driver)=
# Driver and Component Layer

## 1 Overview

OMEGA can be used as either a standalone ocean model or as
the ocean component of E3SM. In either case, OMEGA requires a
top-level driver and interface layer. Here we describe the
requirements and design of both the standalone driver (main)
and the component interface for coupled simulation.

## 2 Requirements

### 2.1 Requirement: Component interfaces

When running as part of a larger coupled model, OMEGA must
supply initialize, run and finalize methods. These routines
must export all variables needed by the parent model while also
importing all fields needed by OMEGA. In most cases, there will
also need to be a thin wrapper to translate data types between
the parent model and internal OMEGA data types.

### 2.2 Requirement: Standalone driver

When OMEGA is used as a standalone model, it must supply a
driver or main routine. For consistency with coupled
simulations, this driver must mimic a parent coupled model,
calling the same init, run and finalize methods and supplying
any needed data (eg surface forcing data).

### 2.3 Requirement: Encapsulation and persistent model state

Communication with the parent model or standalone driver must
be through the component interfaces as method arguments. All
other aspects of the OMEGA model and model state that need to
be retained across component calls or subsequent run intervals
must be stored internally within OMEGA as static variables.

### 2.4 Requirement: Managing environments

Initializing and exiting environments like MPI and Kokkos must
take place at the driver level (standalone driver or coupled
model driver) as these environments are shared by other
components. The finalize method described below must clean up
all of Omega memory and data types so that these environments
can be exited cleanly.

### 2.5 Requirement: Run method

The run method must advance the model one specified time
interval based on a set of inputs (the import state) and must
return a set of outputs (the export state) at the end of the
time interval. In a coupled simulation, this interval will be
the coupling interval. In standalone simulations, this interval
is typically the fastest forcing data interval.

### 2.6 Requirement: Finalize method

The finalize method must provide a graceful exit, checkpointing
as needed and cleaning up all memory. It must not, however,
exit the MPI or other shared environments (eg Kokkos) per
requirement 2.4

### 2.7 Requirement: Init method

An initialization method must initialize all model state, mesh
and other information needed by the model itself or the parent
coupled system. Variables needed by the parent model will be
returned as arguments while all other parts of the model state
will be retained in static variables for later retrieval by the
run method as described in Req 2.3. The model state on initialization
must correspond to the initial time for the simulation integration.

## 3 Algorithmic Formulation

No algorithms are introduced.

## 4 Design

The design is essentially determined by the requirements above.
We define an init, run and finalize method. There will actually
be two layers of these functions. One will be the internal
Omega inteface used by the Omega standalone driver. A second
layer will be needed for translating between internal Omega
data types and E3SM (or other parent model) data types and
ensuring the model is synchronized correctly with the parent.

Within the directory structure of OMEGA, the `src/driver`
directory will contain two subdirectories called standalone
and E3SM. The standalone subdirectory will contain the
standalone driver (obviously) and the E3SM directory will
contain the wrapper interfaces that translate between the
E3SM components and data types and the Omega methods and data
types. The CMake build system will determine which directory
will be used in the build. The actual OcnInit, OcnRun, OcnFinalize
routines described below will reside in the `src/ocean` directory.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

There are currently no additional parameters needed for this level.
Configuration is generally determined by other modules.

#### 4.1.2 Class/structs/data types

For standalone simulation, data types are determined by other modules
(eg state and mesh). No new data types are needed here.

E3SM data types are in flux (from MCT to MOAB). We will add the coupled
model data types here in the future.

### 4.2 Methods

#### 4.2.1 Init

The Init method (OcnInit) will call individual initialization
routines for every module in Omega. On input, it requires the
MPI communicator to use as the master ocean communicator
(`MPI_COMM_WORLD` for standalone, a coupler-partitioned
communicator for coupled simulations). On output, it will
return mesh information, the current model state and the
time instant associated with that state.

The precise interface awaits the design of various other
modules, but will look something like:
```c++
int OMEGA::OcnInit(
   MPIComm Comm, ///< [in] ocean MPI communicator
   OMEGA::TimeInstant &StartTime, ///< [out] sim start time
   OMEGA::TimeInterval &RunInterval, ///< [out] interval for run method
   OMEGA::State &CurrState, ///< [out] current model state
   other args as needed
);
```
The routine will return the initial state, the start time and
the run interval computed based on input options. An integer
return value will be non-zero if errors are encountered and zero
if successful.
It is also possible that a multi-stage initialization may
be needed, especially in coupled mode and would require
additional OcnInit interfaces. This will be determined during
integration with E3SM later.

#### 4.2.2 Run

The run method will advance the model one time interval,
typically the coupling time or the fastest forcing interval.
The interface will look like:
```c++
int OMEGA::OcnRun(
   OMEGA::TimeInstant &CurrTime, ///< [inout] current sim time
   OMEGA::TimeInterval &RunInterval, ///< [in] interval to advance model
   OMEGA::Alarm &EndAlarm, ///< [out] alarm to end simulation
   OMEGA::State &CurrState, ///< [inout] current model state
   other args as needed (eg forcing)
)
```
The model state, current time and the time interval will be input.
Other variables will be needed as well, like the surface forcing
fields, and will be added as needed. On return, the time instant
will contain the end time of the interval and the end alarm will
be ringing if the end of the simulation has been reached. The
CurrState will be the ocean state at the end of the run interval.
An integer error code will be zero if successful and non-zero if an
error was encountered. This interface will be modified as needed to
include other fields.

#### 4.2.3 Finalize

The finalize method will write a checkpoint/restart file
(if not already written by driver or run method on the
final timestep) and then clean up all arrays and classes by
calling the relevant routines for all Omega modules. The
interface is similar to the Init interface:

```c++
int OMEGA::OcnFinalize(
   OMEGA::TimeInstant &CurrTime, ///< [in] current sim time
   OMEGA::State &CurrState, ///< [in] current model state
   other args as needed for restart
);
```
An integer return value will be zero if successful and
non-zero if an error is encountered either writing a
restart or deallocating memory.

#### 4.2.4 Standalone driver (main)

With the interfaces above, the standalone driver should look
something like the code below (details subject to change during implementation).
```c++

int main(int argc, char **argv) {

   MPI_Init(); // initialize MPI
   Kokkos::init(); // initialize Kokkos
   {

   OMEGA::State CurrState;
   OMEGA::TimeInstant CurrTime;
   OMEGA::TimeInterval RunInterval;
   OMEGA::Alarm EndAlarm;

   int Err = OcnInit(MPI_COMM_WORLD, CurrTime, RunInterval, CurrState,
                     EndAlarm, etc);
   if (Err != 0) LOG_ERROR("Error initializing OMEGA");


   while (Err == 0 && !(EndAlarm.isRinging()) ) {

      // call routines for forcing and other inputs

      // call run method
      Err = OMEGA::OcnRun(CurrTime, RunInterval, EndAlarm,
                          CurrState, etc);
      if (Err != 0) LOG_ERROR("Error advancing Omega one interval");

      // Other tasks if needed (eg IO could occur here or within run
      // method
   }

   int Err2 = OMEGA::OcnFinalize(CurrTime, CurrState, etc);
   if (Err2 != 0) LOG_ERROR("Error finalizing Omega");

   int ErrAll = abs(Err) + abs(Err2);
   if (ErrAll == 0){
      LOG_INFO("OMEGA successfully completed");
   } else {
      LOG_ERROR("OMEGA terminating due to error");
   }

   }
   Kokkos::finalize(); // Exit Kokkos
   MPI::Finalize(); // Exit MPI

   return ErrAll;
}

```

#### 4.2.5 Coupler-component interfaces

To be added later

## 5 Verification and Testing

### 5.1 Test forward model

A forward model smoke test is included as part of the CTest
unit test suite. This test runs the standalone model in a
minimal configuration and only tests for successful completion.
Other forward model system testing (eg in Polaris) will be
inherently testing the driver layers.
  - tests requirements 2.2-2.7

### 5.2 Coupled model testing

Once Omega is integrated into E3SM, various E3SM system tests
will be run regularly and will test all coupled interfaces.
We will add an Omega developer test suite to include these tests.
  - tests requirement 2.1
