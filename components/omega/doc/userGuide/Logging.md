(omega-user-logging)=

# Omega Logging

Omega's logging system is built upon the [spdlog](https://github.com/gabime/spdlog)
logging tool. As of this version, It only supports Standalone build mode.

## Logging in Standalone Build mode

To use Omega's logging system in Standalone build mode, you need to include
the Logging.h header file in your code.

```c
#include "Logging.h"
```

You can then use the following macros to write logs of different severity
levels:

```c
LOG_TRACE   (logmsg);
LOG_DEBUG   (logmsg);
LOG_INFO    (logmsg);
LOG_WARN    (logmsg);
LOG_ERROR   (logmsg);
LOG_CRITICAL(logmsg);
```

The "logmsg" argument in these macros can be any string that follows the syntax
of C++'s std::format function. For example, the following code writes a
"Hello, World!" message as a TRACE level log:

```c
LOG_TRACE("Hello, {}!", "World");
```
For more information about the syntax of logmsg, please see the spdlog documentation.

Here is a table that summarizes the different log severities and their meanings:

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


You can replace the string "World" with any variable whose type is
a C++ basic data type, such as an int, a float, or a string. You can
also use a limited set of Kokkos variables. For example, the following code
writes the value of the variable `MyInt` as a TRACE level log:

```c
int MyInt = 10;
LOG_TRACE("The value of MyInt is: {}", MyInt);
```

By default, the logfile will be created in the build directory.

## E3SM Component Build

T.B.D.
