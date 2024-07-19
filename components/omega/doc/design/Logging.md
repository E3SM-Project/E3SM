(omega-design-logging)=
# Logging

## 1 Overview

The OMEGA logging system offers a standardized Application Programming
Interface (API) that facilitates the tracking of events occurring during
the runtime of all Omega sub-components. Each event comprises a descriptive
message along with an indicator denoting its severity level.

## 2 Requirements

### 2.1 Requirement: Basic logging

It is a requirement for the OMEGA logging system to furnish Omega
sub-components with a mechanism to store log messages onto a file system,
based on their corresponding severity levels.

### 2.2 Requirement: Inter-operability with E3SM

It is a requirement for the OMEGA logging system to recognize the log root
filename defined in E3SM as a part of the input configuration.

### 2.3 Requirement: Multi-processing

It is a requirement that each process should have the capability to generate
log messages independently. Additionally, there should be a mechanism in
place to control the number of processes generating log messages
simultaneously. This functionality is desired to ensure efficient management
of system logs and prevent potential system overloading due to excessive
logging activities.

### 2.4 Requirement: Multi-threading

It is a requirement for the OMEGA logging system to provide thread-safe
logging capabilities. This functionality is necessary to ensure that
concurrent logging activities by multiple threads are executed in
a synchronized manner, without causing data corruption or inconsistencies
in the generated log files.

### 2.5 Requirement: Omega data types

It is a requirement for the OMEGA logging system to support a subset of
Omega data types, in order to generate well-formatted log messages that
include these data types. This functionality is crucial to ensure that all
relevant information is captured in the logs and to facilitate effective
debugging and troubleshooting activities.

### 2.6 Requirement: Buffered file input/output

It is a requirement for the OMEGA logging system to support both buffered
and unbuffered file IO modes. Buffered mode reduces the overhead associated
with file IO operations, improving performance and efficiency. Unbuffered mode
is crucial for debugging purposes, as it allows for immediate writing of log
messages to the file, preventing potential loss of unwritten messages that
could lead to confusion. The API should provide an option for the user to
choose between buffered and unbuffered modes.

### 2.7 Desired: Formatter

It is desired for the OMEGA logging system to support various use-cases of log
files by generating log messages through one of the installed message
formatters. These formatters should be switchable at both compile-time and
runtime to accommodate changing requirements. The default formatter should
generate an unstructured text output, while other possible formatters should
be able to generate structured output such as JSON. This functionality is
essential to ensure flexibility in the logging system, allowing users to
tailor the output format to their specific needs.

### 2.8 Desired: Output sinks

It is desired for the OMEGA logging system to be capable of transporting log
messages to various output sinks such as emails or databases. Thisfunctionality
is crucial for ensuring the accessibility and availability of log data, enabling
users to view and analyze the logs through multiple channels.

### 2.9 Desired: Local logging

It is desired for the OMEGA logging system to provide each sub-component with
its own local logger, which can be used to set parameters such as "severity,"
"message format," and "sink type" independently from other sub-components. This
functionality is necessary to ensure that each sub-component has full control
over its own logging activities, allowing for more granular and targeted logging
output.

## 3 Algorithmic Formulation

No algorithm is necessary for the logging system.

## 4 Design

The spdlog logging library has been chosen as the core logging capability for
the OMEGA logging system. spdlog is a feature-rich, header-only C++ logging
library that is widely used and highly regarded in the industry. It offers
a comprehensive set of features that satisfy most of the requirements, but
additional E3SM and Omega specific requirements must be addressed.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

The software design will implement a global switch for setting the minimum log
severity and message format using predefined global C++ defines. The severity
level can be set at compile-time using `-D OMEGA_LOG_LEVEL=<level>`, where `<level>`
can be one of the predefined levels: `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`,
or `CRITICAL`, where the severity increases from `TRACE` to `CRITICAL`. Log
messages that have severity equal to or higher than the specified level will be
stored in the log destination.

Furthermore, the message format can be set using `-D OMEGA_LOG_PATTERN=<pattern>` at
compile time, where `<pattern>` represents a format in the form of `%flag`,
similar to the strftime function. This pattern defines the layout of the log
message.

Users can control which MPI ranks generate log files using
`-D OMEGA_LOG_TASKS=<tasks-pattern>` at compile time. The `<tasks-pattern>` is
either all for all tasks to generate log files, or comma-separated MPI rank
numbers, or a range of MPI ranks with a dash. For example,
`-D OMEGA_LOG_TASKS=0,2-3` indicates that MPI ranks 0, 2, and 3 generate log
files.

#### 4.1.2 Class/structs/data types

No public data type is necessary for the logging system.

### 4.2 Methods

The OMEGA logging system includes multiple macros that allow for logging at
various severity levels. These macros guarantee that log messages with severity
levels lower than the specified `OMEGA_LOG_LEVEL` will not be compiled.

#### 4.2.1 Log write

The OMEGA logging system includes macros that allow the user to write log
messages with different severity levels. These macros follow the naming
convention `LOG_<severity level>`. For instance, the following macro can be
used to write an informational log message:

```c++
LOG_INFO("log message");
```
The OMEGA logging system provides the following macros for writing log messages
with different severity levels:

1. `LOG_CRITICAL`: Used to indicate a critical error or issue that requires
   immediate attention and could potentially cause data loss or application
   failure.
2. `LOG_ERROR`: Used to indicate an error or issue that could cause
   application instability or incorrect behavior.
3. `LOG_WARN`: Used to indicate a warning or potential issue that should be
   addressed but does not necessarily require immediate attention.
4. `LOG_INFO`: Used to provide general informational messages about the
   application's operation or progress.
5. `LOG_DEBUG`: Used to provide detailed debugging information for
   developers to track down issues or analyze behavior.
6. `LOG_TRACE`: Used to provide the most detailed level of information,
   including function calls and variable values, and is typically used for deep
   analysis and debugging.

In order to enhance the logging capabilities, OMEGA logging system provides
each macro with an optional argument of an Omega data type. If the data type
is recognized, it will be formatted automatically by the system and appended in
the appropriate format to the log message. As an example, the following macro
demonstrates how an argument of type `Array1DR8` can be appended to a log
message:

```c++
Array1DR8 var;
...
LOG_INFO("log message {}", var);
```

## 5 Verification and Testing

### 5.1 Test OMEGA logging capabilities

A unit test should be performed to verify all requirements, except for 2.2
(E3SM interoperability). The unit test should create multiple MPI ranks with
multiple threads to ensure that requirement coverage is complete. The unit
test should also include tests for both buffered and unbuffered logging
capability, with various log message sizes and frequencies, to ensure that both
modes work correctly and without loss of data.

### 5.2 Test OMEGA data type support

A unit test should be performed to verify Requirement 2.5(Omega data types).

### 5.3 Test E3SM inter-operability

Requirement 2.2 will be tested as part of larger system tests and will not have
a dedicated unit test.
