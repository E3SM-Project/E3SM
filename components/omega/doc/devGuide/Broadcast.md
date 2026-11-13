(omega-dev-broadcast)=

## Broadcast

In a parallel environment, it is sometimes necessary to broadcast values
from a single task to all remaining tasks. Omega provides broadcast
functions for scalars and vectors of all supported data types with the
lone exception of a vector of strings. Scalar strings are supported and
can be used to broadcast a string vector element-wise. All functions are
simple wrappers around the relevant `MPI_Bcast` function.

To broadcast a scalar of any type, there are several supported interfaces:
```c++
Broadcast(MyScalar);
Broadcast(MyScalar, FromTask);
Broadcast(MyScalar, AltMachEnv);
Broadcast(MyScalar, AltMachEnv, FromTask);
```
The first broadcasts a scalar from the master task to all other tasks in
the default machine environment (set of MPI tasks). The second allows the
user to broadcast from a task other than the master task, but still in the
default environment.  The final two interfaces perform the same functions
but operate only in the set of tasks belonging to a different MachEnv than
the default. There are type-specific implementations for each, but all are
aliased to these same generic interfaces and the combinations above are
implemented using a combination of optional arguments and specific interfaces.

Omega also provides for the broadcasts of a `std::vector` of any supported
data type except strings. The interfaces are identical to those for scalars:
```c++
Broadcast(MyVector);
Broadcast(MyVector, FromTask);
Broadcast(MyVector, AltMachEnv);
Broadcast(MyVector, AltMachEnv, FromTask);
```
The implementation of a boolean vector is slightly different from the others
because a `std::vector<bool>` uses an optimal packing that does not allow
certain assumptions of storage for the MPI call. Instead, the boolean vector
is converted to/from an integer vector within the broadcast routine and
uses an MPI integer broadcast. This is hidden from the user. There is no
interface for a vector of strings due to the complexity of extracting
variable-length strings into a character vector for the MPI call. In the
unlikely event such a broadcast is needed, users should simply loop through
the vector and broadcast each scalar string individually.

All of the interfaces above use the blocking MPI broadcast and include an
implicit synchronization. For most use cases, this is adequate. However,
non-blocking broadcasts may be added in the future should the need arise.
