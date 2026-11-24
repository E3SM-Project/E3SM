(omega-design-error)=
# ErrorHandler (Error)

## 1 Overview

Any simulation code requires a means to check for errors,
report errors and abort if the error is critical. Omega
includes a logging facility that can perform some of these
functions. However, we desire a facility that can check
for errors in a manner that improves code readability and
that also provides a stack trace to better understand the
path taken that resulted in the error.

## 2 Requirements

### 2.1 Requirement: Logging errors

We require the ability to log an error with an error message
that includes the routine and line number. Note that the
current Omega Logging facility already supports this and
should be leveraged for this purpose.

### 2.2 Requirement: Levels of error severity

Omega must support several levels of error severity. The
Logging facility already supports a common hierarchy with
debug, info, warn, error and critical levels. Because the
debug and info simply write information, they will continue
to be managed by the Logging interfaces.

The warn level may be shared by both the Logging facility and
error handler. For example, if only a warning message is required,
the Logging facility should manage the message. If we wish to
inform the calling routine that an error at the warning level
has occurred, the Error handler should manage the warning.

The error handler should be the primary interface for the more
severe errors at the error and critical levels. In a canonical
hierarchy, critical errors are those that cause the application
to fail, while the error level refers to errors that are severe
enough to be addressed but still allow the application to proceed
without crashing. In Omega, these latter errors typically cause
incorrect or unintended simulation results and are considered
an application failure. We prefer to use critical errors that
abort the simulation for most such errors. The error level
should be reserved for cases in which an error is encountered
but we must pass error information back to the calling routine
(eg through a return code) to judge the response to the error.

Note that at the initial creation of this document, many of
the current errors reported by Omega are treated at the error
level and should be promoted to critical errors following the
description above.

### 2.3 Requirement: Abort on error

When a critical error is encountered, we must have the ability
to abort the simulation.

### 2.4 Requirement: Generate stack trace

The error log must include a stack trace to determine and
output the calling sequence that resulted in the error being
reported.

### 2.5 Requirement: Identify task id

We require a means of determining which parallel task is
encountering this error.

### 2.6 Requirement: Output location

We require errors and the related stack trace to be logged
in an error log file. A separate error log file for each MPI
task may be desirable both for ease of implementation and to
better satisfy requirements 2.4 and 2.5 above. An abort error
message should also be sent to the default log file when
possible so that users know to look for an error log for
details.

### 2.7 Requirement: Success/fail error codes

The error handler should include a public success error code
(typically zero) to standardize a successful return code for
Omega routines. Any return code not equal to success should
be interpreted as a failure. However, a standard Fail code
can also be provided as a default failure code.

### 2.8 Requirement: Error messages and formatting

We require the ability to define an error message that will
be output to an error log for any error encountered. We also
require the ability to include some variable or run-time
information as part of the the logged message to better
describe the cause of the error, similar to the existing
Logging capability (eg "Value of MyVar is {} but expected
{}"). Ideally, the formatting would be identical to that
provided by the Logging facility and would leverage that
facility for the logging of error messages.

### 2.9 Requirement: Accumulated error

In some cases, we wish to accumulate (sum) non-critical errors
encountered within a routine rather than returning or aborting
immediately. For example, a unit test might prefer to report
all errors rather than aborting on the first error encountered.
Among other things, this is a requirment that error codes must
all have the same sign to avoid summing to zero (the nominal
success code).

### 2.10 Desired: Compact error checking

We desire error checking to be very concise - ideally a
single line of code - so that the error checking does
not impair code readability.

### 2.11 Desired: Registering error codes and messages

To facilitate 2.10 or to provide unique error codes, it may
be desirable to pre-register an error code and related error
message during initialization of a module so that it can
be used later and referred to simply by the assigned code.

## 3 Algorithmic Formulation

No algorithms are needed.

## 4 Design

We will build the error handler on the existing Logging
functionality, the spdlog library and the cpptrace library
that provides a stack trace capability. Similar to the
Logging facility, we will make use of defined cpp macros
for many interfaces. These macros allow us to both condense
multi-line patterns into a single-line interface, but also
allow us to incorporate the LINE and FILE macros at the
actual line where the error is encountered.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

There will be two public integer parameters (constexpr)
``Err::Success`` and ``Err::Fail``. Success will be the only
identifier for success and will be equal to zero. Any error code
not equal to success will denote failure. This allows multiple
failure error codes and accumulation of error. The ``Err::Fail``
code is defined to provide a default failure code if one is
not defined or if a specific value is not required.

#### 4.1.2 Class/structs/data types

There will be an Error class that primarily holds defined
error codes and associated messages.
```c++
class Error {
   public:

      static constexpr int Success = 0;
      static constexpr int Fail    = 1;

      // public methods (see below)

   private:

      // container for pre-defined error codes/messages
      static std::map<int, std::string> AllCodes;

      // private methods (see below)
};
```

### 4.2 Methods

#### 4.2.1 Macros

We provide a number of macros for error checking. As described
in section 2.2, we expect the most common interface should be
the critical error macro:
```c++
ERROR_CRITICAL(ErrCode, ErrMsg, additional args for ErrMsg)
```
The error message (ErrMsg) can be either a simple string or
can follow the format of the Logging facility in which ``{}``
placeholders can be used with additional arguments to
include variable information
(see the [Logging documentation](#omega-user-logging)).
If the input ErrCode is ``Err::Success`` nothing happens and
the code continues. If the ErrCode matches a code predefined
using the ``defineCode`` function described below, the ErrMsg
from that predefined code is used with any of the additional
arguments to fill in any ``{}`` placeholders. Finally, if the
ErrCode is not equal to Success and is not a predefined code,
the ErrMsg is used, again with the subsequent arguments to print
an error message. For the failure modes, the macro also prints
a stack trace using cpptrace calls and then aborts the
simulation using the ``Err::abort`` function described below.
Typical use cases might be:
```c++
// Checking an error condition
if (MyFailCondition) ERROR_CRITICAL(Err::Fail, "MySimpleMessage");

// or
if (MyFailCondition)
   ERROR_CRITICAL(Err::Fail, "Variable {} has value {}",
                  MyVarName, Value);

// or if the code is predefined
if (MyFailCondition) ERROR_CRITICAL(DefinedCode, " ");

// Checking a return code from another routine
Err = MyFunction(...);
ERROR_CRITICAL(Err, "Error encountered in MyFunction");

```

Two interfaces are provided in which a condition can be checked,
rather than the error code:
```c++
ERROR_REQUIRE(condition, ErrMsg, additional args for msg);
ERROR_ASSERT(condition, ErrMsg, additional args for msg);
```
where condition is any conditional expression that evaluates to
true/false. These have the same behavior as the critical errors
above with a false condition treated as an ``Err::Fail``. The
REQUIRE interface is always evaluated; the ASSERT interface behaves
similarly to the C++ assert function and is only checked for DEBUG
builds.

For less severe errors, we provide several different macros
to match the behavior desired when continuing the simulation.
The first interface is similar to the critical interface above:
```c++
ERROR_CHECK(ErrCode, ErrMsg, addition args for message);
ERROR_WARN(ErrCode, ErrMsg, addition args for message);
```
which simply prints the message when the ErrCode is not
equal to Success. The ``ERROR_WARN`` interface is used when
only a warning message is desired. Predefined
error codes can be used here as well. Unlike the
critical interface, this does not print a stack trace or abort
the simulation and the behavior is more similar to the
``LOG_ERROR`` or ``LOG_WARN`` macros.

If a routine has a return code and a return-on-error is
desired, we provide an interface:
```c++
ERROR_RETURN(ReturnCode, ErrCode, ErrMsg, added vars if needed);
ERROR_RETURN_WARN(ReturnCode, ErrCode, ErrMsg, added vars);
```
This behaves similar to the ``ERROR_CHECK`` except that after
the message is printed, the error code is added to the return
code and a return statement is executed to exit the routine.
A use case might be:
```c++
if (MyFailCondition) ERROR_RETURN(ReturnCode, Err::Fail, "MyErrMsg");

// which is equivalent to:

if (MyFailCondition) {
   LOG_ERROR("MyErrMsg");
   ReturnCode += Err::Fail;
   return ReturnCode;
}
```

Finally, a macro is provided to accumulate errors:
```c++
ERROR_SUM(ErrSum, ErrCode, ErrMsg, added vars if needed);
```
which after printing the messages, adds the most recent
error code to the running sum stored in ``ErrSum``.

Note that in the latter two cases, care must be taken so that
the accumulated error or return codes do not match a
predefined code and cause confustion in a later use of the
code in an ``ERROR_CRITICAL`` or ``ERROR_CHECK`` call.

#### 4.2.2 Other Error methods

Some additional functions are needed to complement the macros
described above.

To allow more concise error checking, developers can pre-define
and error code and associated error message using:
```c++
Error::defineCode(ErrCode, ErrMsg);
```
This routine addes the ErrMsg to a map of predefined messages
and returns the new ErrCode assigned to that message for later
use in the macros described above. This function should be
called at some point during model (or module) initialization
and the ErrCode saved for use during the run sequence.

The abort routine:
```c++
Error::abort();
```
is provided to abort a simulation. A minimal version would
simply call the MPI abort function, but it is possible to
add additional cleanup tasks.

A finalize routine:
```c++
Error::finalize(ErrCode);
```
will clean up all of the predefined error codes. If the
ErrCode is not equal to success, it will set the error code
to the standard abnormal termination exit code (1) and avoid other
standard c++ exit codes for other run-time errors.

## 5 Verification and Testing

### 5.1 Test all failure modes

A unit test will test all non-critical error interfaces first
and then end with a critical failure. The ``WILL_FAIL`` property
in ctest will ensure the expected unit test failure will pass.
An inspection of the log and error files will be used to verify
messages are printed as expected.
  - tests all requirement
