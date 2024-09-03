(omega-dev-driver)=

## Driver

The Omega Driver consists of methods to run the Omega ocean model in either
standalone mode or as a component of E3SM. In both cases, the execution of Omega
is divided into three phases, carried out by the `ocnInit`, `ocnRun`, and
`ocnFinalize` methods. An interface layer to translate data types of the parent
model to and from internal Omega data types will also be needed for running
Omega as a component. For standalone mode, a `main` function is defined to
serve as the top-level driver. The driver methods will evolve as more modules
are developed and as Omega is integrated with E3SM.

### ocnInit

The `ocnInit` method initializes everything needed to run the Omega ocean model,
and has the following interace:
```c++
int ocnInit(MPI_Comm Comm, Calendar &OmegaCal, TimeInstant &StartTime, Alarm &EndAlarm);
```
It requires an MPI communicator as input, as well as a yaml configuration file
named `omega.yml` in the directory Omega is launched from. Default configuration
files are currently kept in the `configs` directory, eventually an automated
tool will extract default configuration options from module header files and
build the default config file. The `ocnInit` method will read the configuration
options into the static `ConfigAll` object. The time management options that
determine the duration of the simulation are read from the Config object and
stored in `OmegaCal`, `StartTime`, and `EndAlarm`, which are output by the
method. The `init()` methods for each Omega module necessary to run the ocean
model are called. An integer error code is returned.

### ocnRun

The `ocnRun` method advances the model over the interval defined by time
management options, integrating from the start time until the `EndAlarm` is
ringing. The interface is:
```c++
int ocnRun(TimeInstant &CurrTime, Alarm &EndAlarm);
```
The `TimeStep` is retrieved from the default `TimeStepper` object initialized
during `ocnInit` and used with the value of `CurrTime` on input to construct
a `Clock` object which manages the time loop. The `EndAlarm` is attached to
the `Clock` to signal the end of the time loop. Eventually, there will be
periodic alarms also attached to the `Clock` to trigger forcings, writing
output, restarts, etc. As the time loop advances, the `doStep` method of the
`TimeStepper` object is called to update the ocean model until the `EndAlarm`
is triggered. The updated `CurrTime` and status of `EndAlarm` are returned,
along with an integer error code.

### ocnFinalize

The `ocnFinalize` method is needed to clean up all the objects allocated by the
model. The interface is:
```c++
int ocnFinalize(TimeInstant &CurrTime);
```
The `CurrTime` is input in case a restart file needs to be written. An integer
error code is returned.
