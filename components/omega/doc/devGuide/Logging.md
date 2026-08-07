(omega-dev-logging)=

# Developing Omega Logging

Omega's logging system is built upon the
[spdlog](https://github.com/gabime/spdlog) logging tool.

Logging macros and custom formatters for Kokkos data types have been
integrated into the Omega logging system via the Logging.h header file
located in the src/infra directory. Users who wish to utilize Omega's
logging capabilities must include this header file in their code.

The `src/infra/Logging.cpp` file encompasses various functions related
to logging. As of the current version, this file contains the
`OMEGA::initLogging` function, which initialize the logging process.
It is recommended to invoke the `OMEGA::initLogging` function at the beginning
of an Omega application to properly initialize the logging system.

## Initializing Omega Logger

The `OMEGA::initLogging` function, located within the `src/infra/Logging.cpp`
file, serves as a pivotal component for initializing the Omega logging system.

The function establishes the default logger configuration, ensuring that
logging messages are effectively saved to a designated file. The path to
this file can be passed as the second argument to `OMEGA::initLogging`. When
it is omitted, the path comes from the `OMEGA_LOG_FILE` environment variable
if that is set and non-empty, and from `OMEGA::OmegaDefaultLogfile`
(`omega.log`) otherwise. The log file is truncated when it is opened, so it
holds only the messages from the current run.

## Log files in the unit tests

Every unit test driver calls `OMEGA::initLogging(DefEnv)` without a file name,
so without any further arrangement all of the tests would write to a single
`omega.log` in the build's `test` directory, with nothing to say which test
wrote which line. To avoid this, `test/CMakeLists.txt` sets `OMEGA_LOG_FILE`
per test in `add_omega_test`, so that each test writes to
`test/logs/<TEST_NAME>.log`. Note that the name comes from the CTest test
name rather than the source file, because several tests are built from the
same source (for example `TEND_PLANE_TEST` and
`TEND_PLANE_SINGLE_PRECISION_TEST` both come from `TendencyTermsTest.cpp`).

A new test added with `add_omega_test` gets its own log file automatically. A
test driver that hard-codes a log file name would opt out of this, so pass no
file name to `initLogging` unless the test manages its own log files, as
`LoggingTest` does.

The environment variable relies on the MPI launcher passing its environment
to the ranks it starts, which is true of `srun`. If a launcher does not, the
tests fall back to writing `omega.log`.

## Creating Logging Macros

The Omega logging macros, denoted by the prefix `LOG_`, are defined within
the `src/infra/Logging.h` file. These macros constitute a main part of
the logging infrastructure, enabling users to seamlessly incorporate
logging functionality into their Omega applications.

Furthermore, the logging framework includes a distinct set of macros that
commence with the prefix `LOGGER_`. These macros offer enhanced
versatility by accommodating the utilization of the spdlog logger as
their first argument. This approach facilitats the integration of various
logger types.

## Customer formatter for Kokkos

Within the same header file, you will encounter specialized spdlog formatter
structs designed to accommodate Kokkos data types.

For further information on customizing the spdlog formatter, refer to
[Custom formatting](https://github.com/gabime/spdlog/wiki/3.-Custom-formatting).
