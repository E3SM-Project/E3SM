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
this file is determined by utilizing the `OMEGA_LOG_FILEPATH` macro, which
allows users to specify the desired file location for logging purposes.

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
