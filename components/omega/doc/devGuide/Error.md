(omega-dev-error)=

# Error Handler

The error handling facility in Omega provides a number of macros and an
error class to determine behavior when encountering errors and to provide a
standard error return code that routines can use. Macros are used to provide
a concise single-line interface when checking errors and to ensure the error
message reports the correct line number.  It builds on the
[Logging](#omega-dev-logging) facility for writing the error messages and
for some interfaces will generate a stack trace using the
[cpptrace](https://github.com/jeremy-rifkin/cpptrace) library. All interfaces
can be used after including the `Error.h` and after initializing the Logging
facility (typically done as early as possible in Omega initialization).

The most common use case is to check for critical errors using the
`ABORT_ERROR` macro. For example:
```c++
  if (condition) ABORT_ERROR("Bad value for variable: {}", Value);
```
The macro takes a string argument for the error message and allows values
to be inserted where `{}` placeholders exist, similar to the Logging macros.
This macro writes the error message and aborts the simulation. If Omega has
been built with debug information available (either a debug build or the
-g option has been passed to the compiler), the macro also creates and writes
a stack trace to provide the calling sequence that generated the error.

Closely related to this critical error macro are the assert and require
macros:
```c++
OMEGA_ASSERT(Condition, Message, addition args for message);
OMEGA_REQUIRE(Condition, Message, addition args for message);
```
These emulate the behavior of the C++ assert function in which the condition
is checked and a critical error is thrown if the condition is not met (false).
The `ASSERT` macro is only evaluated in debug builds while the `REQUIRE`
macro is always evaluated. These have identical behavior to the critical error
macro, but take the additional boolean argument of the condition to be checked.

The remaining macros operate on or with an Error class that contains an error
code and an error message. These can be used as a return code, as in:
```c++
Error MyError = callMyFunction();
```
The error code within the class is an ErrorCode enum that defines the severity
of the error:
```c++
enum class ErrorCode {
   Success = 0, /// Standard success code
   Warn,        /// Standard warning code
   Fail,        /// Standard failure code
   Critical     /// Standard critical error code
}
```

Errors can be accumulated using either the `+` or `+=` operators. This is
useful especially in test drivers where we may not want to exit when an
error is encountered or if the error is only a warning.
```c++
Error MyError; // defaults to success
MyError += myRoutine1();
MyError += myRoutine2();
```
When adding or accumulating errors, the resulting error code is the max (most
severe) error code and the error messages are concatenated (separated by a
newline character) so that all messages are retained.

For routines that return an error, we provide an `RETURN_ERROR` macro. For
example,
```c++
Error ReturnVal; // initializes to success
[do stuff]
if (Condition) RETURN_ERROR(ReturnVal, ErrorCode::Fail, ErrMsg);
```
This macro accumulates the new error (Code, ErrMsg) into the return code
and returns control and the error to the calling routine. The error message
is _not_ printed by this macro; the calling routine must decide what action
to take with the returned error code and message.

To check these returned error codes or any other accumulated errors in a code,
the following macros can be used:
```c++
CHECK_ERROR(ErrorToBeChecked, NewMessage, MsgArgs...);
CHECK_ERROR_WARN(ErrorToBeChecked, NewMessage, MsgArgs...);
CHECK_ERROR_ABORT(ErrorToBeChecked, NewMessage, MsgArgs...);
```
These macros check the error and if it is not Success, prints the new message
to the log file along with any other accumulated messages associated with the
error. In the first two cases, the error is then reset (back to success
with an empty message) and the simulation continues. In the third case,
the new message is treated as a critical error, so after printing all the
messages, it prints a stack trace and aborts.

While the macros above cover most use cases, it is possible to use the
Error class functions directly to customize the treatment of particular errors.
For checking error codes, there are two useful member functions:
```c++
MyError.isSuccess();
MyError.isFail();
```
that check whether the code is `ErrorCode::Success` or `ErrorCode::Fail`.
These functions check for equality so it is sometimes more useful to check
their negation, as in
```c++
if (!MyError.isSuccess()) { do something; }
```
This example will pick up all non-successful conditions (warn, fail, critical).

The default constructor is often used (eg in the examples above) to create
and initialize an Error instance for use as a return code.  Three other
constructors in addition to this default constructor are used by the
macros and examples above:
```c++
Error MyError; // default constructor with success code and empty message
Error MyError(ErrorCode); // constructor with error code only and empty msg
Error MyError(ErrorCode, ErrorMessage); // constructor with both code and msg
Error MyError(ErrorCode, Line, File, ErrorMessage); // as above but adds prefix
```
In the last two cases, the constructor initializes the error message either
with or without a prefix. In cases where the message is expected to be
accumulated without an immediate print (eg the `RETURN_ERROR` macro) the prefix
case is used and the prefix includes the message severity and the routine name
and line number where the error occurs (the cpp `__LINE__` and `__FILE__` can
be used). In the other constructor with a message, the prefix is left off,
typically for cases where the relevant LOG function will be called underneath
which adds the approporiate prefix as part of Logging.

In the `CHECK_ERROR` macros, the Error is reset to Success and a blank message
using the function
```c++
MyError.reset();
```

Finally, the Error class includes the abort function `Error::abort()` which
aborts the simulation. Currently, the function simply prints a generic aborting
message and calls the `MPI_abort()` function to kill the simulation, but
additional functionality may be added later to clean up before shutting down.
The abort function assumes it is being called from one of the Error functions
or macros so does not print any previous accumulated errors and does not
print a stack trace. It is best to use one of the ABORT macros above rather
than calling the abort member function directly.
