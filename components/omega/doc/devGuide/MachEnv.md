(omega-dev-mach-env)=

## Machine Environment (MachEnv)

The MachEnv class maintains a number of parameters associated with
the parallel machine environment. These include message-passing parameters
like MPI communicators and task information, the number of threads
if threaded, vector length for CPUs, and GPU and node information as
needed. Multiple environments can be defined to support running portions
of the model on subsets of tasks.

On model initiation, the initialization routine `OMEGA::MachEnv::init()`
is called to set up the default MachEnv, which can be retrieved at
any time using:
```c++
OMEGA::MachEnv DefEnv = OMEGA::MachEnv::getDefault();
```
that returns a pointer to the default environment.
Once an environment has been retrieved, the individual data members
are retrieved using various get functions:
```c++
  MPI_Comm Comm   = DefEnv->getComm();
  int MyTask      = DefEnv->getMyTask();
  int NumTasks    = DefEnv->getNumTasks();
  int MasterTask  = DefEnv->getMasterTask();
```
There is also a logical function `isMasterTask` that can be used
for work that should only be done on the master task. By default,
the master task is defined as task 0 in the environment, but if
the master task is overloaded, there is a `setMasterTask` that can
redefine any other task in the group as the master.

If OMEGA has been built with OpenMP threading, a `getNumThreads`
function is available; it returns 1 if threading is not on.
The MachEnv also has a public parameter `OMEGA::VecLength` that can
be used to tune the vector length for CPU architectures. For
GPU builds, this VecLength is set to 1.

As noted previously, additional environments can be defined for
subsets of a parent environment. There are three constructor
interfaces for creating an environment:
```c++
  OMEGA::MachEnv MyNewEnv(Name, ParentEnv, NewSize);
  OMEGA::MachEnv MyNewEnv(Name, ParentEnv, NewSize, Begin, Stride);
  OMEGA::MachEnv MyNewEnv(Name, ParentEnv, NewSize, Tasks);
```
An optional additional argument exists for all three that can supply an
alternative task to use as the Master Task for message passing. If not
provided, the master task defaults to 0. The first interface above
creates a new environment from the first `NewSize` contiguous
tasks in the parent machine environment. The second creates a new
environment from a strided set of tasks starting at the `Begin` Task
(eg all odd tasks would have Begin=1 and Stride=2). The final
interface creates a subset containing the selected tasks stored in
a vector `Tasks`. Each new environment is given a `Name` and can be
retrieved by name using `OMEGA::MachEnv::getEnv(Name)`. Once retrieved,
the other get functions listed above are used to get each data member.
Because the new environments are subsets of the default environment, one
additional function `OMEGA::MachEnv::isMember()` is provided so that
non-member tasks can be excluded from calculations. The retrieval functions
call from non-member tasks will return invalid values. Finally, there is
a `removeEnv(Name)` function that can delete any individual defined
environment by name and a `removeAll()` function that cleans up by
removing all defined environments.

As a class that is basically a container for the environment parameters,
the implementation is a simple class with several scalar data members and
the retrieval/creation functions noted above. To track all defined
environments, a c++ map container is used to pair a name with each
defined environment.
