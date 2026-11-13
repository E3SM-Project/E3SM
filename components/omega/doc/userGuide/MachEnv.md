(omega-user-mach-env)=

## Machine Environment (MachEnv)

Within Omega, many aspects of the machine environment are stored
in a class called MachEnv. These include message-passing parameters
like MPI communicators and task information, number of threads
if threaded, vector length for CPUs, GPU and node information as
needed. Multiple environments are supported in case portions of the
code need to run on subsets of tasks or in different contexts. A
default environment is defined early in the initialization of the
model and can be retrieved as described in the
[Developer's Guide](#omega-dev-mach-env).

The user is not expected to set any parameters in MachEnv. All
quantities are derived from either the job launch command
(eg mpirun or srun) that defines the number of MPI tasks and tasks
layouts or from machine parameters enforced during the build based
on supported machine xml configurations.
The latter include the pre-processing parameters
`-DOMEGA_VECTOR_LENGTH=xx` and `-DOMEGA_THREADED` that define an
optimal vector length for CPU code and turn on OpenMP threading
if desired.
