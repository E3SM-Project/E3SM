(omega-user-broadcast)=

# Omega Broadcast

## Broadcasting Functions

Omega's broadcasting functions offer developers a streamlined method for
disseminating scalar and array values among MPI tasks. These functions are
categorized into two types: Blocking and Non-blocking.

### Blocking Broadcasting Functions

When it comes to MPI broadcast operations, a blocking broadcast would imply
that the function call to broadcast data to all processes in a group does
not return until the data has been sent and correctly received by all the
target processes.

### Including the Broadcasting Functions in Your Code

To leverage Omega's broadcasting functions, you must include the Broadcast.h
header file in your source code:

```c
#include "Broadcast.h"
```

### Broadcasting Scalar or Array Values

For broadcasting either scalar or array values, utilize the Broadcast
function as demonstrated below:

```c
OMEGA::I4 MyVal = 1;

// Broadcasting from the master task
OMEGA::Broadcast(MyVal);
```
The above example illustrates the broadcasting of an OMEGA::I4 type scalar
value from the master task of the default Omega environment.

### Specifying broadcasting task and/or environment in Omega

For developers seeking to modify the sending task or the Omega environment
during broadcasting, Omega offers flexible syntax options. These additional
arguments to the Broadcast function enable customization, while maintaining
the simplicity of the basic broadcasting operation.

```c
OMEGA::I4 MyVal      = 1;
const int RootTask   = 1;

// Get a specific Omega Machine Environment
OMEGA::MachEnv *SubsetEnv = OMEGA::MachEnv::get("Subset");

// broadcast from the master task of SubsetEnv environment
OMEGA::Broadcast(MyVal, SubsetEnv);

// broadcast from task 1 of SubsetEnv environment
OMEGA::Broadcast(MyVal, SubsetEnv, RootTask);

// broadcast from task 1 of the default environment
OMEGA::Broadcast(MyVal, RootTask);
```

### Broadcasting Array Values

The syntax for broadcasting arrays is analogous to that used
for scalar values:

```c
std::vector<OMEGA::R8> MyVector;

for (int i = 1; i <= 5; i++) {
    MyVector.push_back(1.0);
}

OMEGA::Broadcast(MyVector, RootTask);
```

This example shows how to broadcast an array of OMEGA::R8 type values,
specifying the root task for broadcasting.

### Supported Data Types for Broadcasting

The Broadcast functions are compatible with the following data types:

* `OMEGA::I4`
* `OMEGA::I8`
* `OMEGA::R4`
* `OMEGA::R8`
* `bool`
* `std::string` (NOTE: array of strings are not supported)

## Non-blocking Broadcasting Functions

This feature is currently under development and will offer asynchronous
broadcasting capabilities.
