(omega-user-error)=

# Error handling

Omega includes an error handling facility that checks for potential errors
and responds to those errors when encountered.
It uses the [Logging](#omega-dev-logging) facility for writing the error
messages so its behavior follows that of Logging. In particular, the error
messages will be output to a single log file (omega.log by default, see the
`OMEGA_LOG_FILE` environment variable) unless logging on more than one MPI
task is requested with `OMEGA_LOG_TASKS`, in which case each logging task
writes its own file. Each message
will have the format
```
[*level] [file:line number] error message
```
where level is the severity (warn, error, critical) and file and line number
refer to the source code line where the error occurs. The error facility
sometimes accumulates multiple error messages before aborting. In
these cases, the most recent error in these accumulated errors will be printed
first and the other error messages in the chain will be printed in reverse
order and indented.  For critical errors, a full stack trace is generated
if debug information is available, either with the `-g` compiler
option or a DEBUG build using `-DOMEGA_BUILD_TYPE=Debug` during the Cmake
configuration.
