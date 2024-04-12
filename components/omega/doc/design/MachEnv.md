(omega-design-machine-env)=
# MachineEnv

## 1 Overview

On startup, OMEGA will need to initialize the message-passing and other
environments and set up parameters related to machine layout, messaging,
and hardware for use throughout OMEGA.


## 2 Requirements

### 2.1 Requirement: Initialize MPI in Standalone

In standalone mode, OMEGA will need to initialize the MPI environment
and define a default communicator.

### 2.2 Requirement: Create MPI communicator in coupled mode

When running as part of a coupled system, OMEGA will need to define a default
communicator based on a parent communicator sent by the calling routine
(coupler or parent model).

### 2.3 Requirement: Define MPI layouts

Each MPI rank will need to know its own rank id, number of ranks and
define a master rank.

### 2.4 Requirement: Kokkos initialization

Since we are using Kokkos for kernel launching and array types, we
will need to initialize Kokkos in standalone mode. It may also be
needed for coupled simulations.

### 2.5 Requirement: Vector blocking size defined at compile time

To achieve the best CPU performance, especially when using
GPU-friendly loop forms, it is useful to explicitly size inner
loops based on a compile-time length (chunk size) that is a multiple
of the vector length.

### 2.6 Desired: Set alternative master task

While rank 0 is typically used as a master task, it is sometimes desirable
to assign a different rank as master to avoid overloading rank 0, especially
when running in coupled mode with other components also using the same rank.

### 2.7 Desired: Multiple environments

In OMEGA, we may wish to run sub-components on different partitions. For
example, we might want to rearrange the communication-dominated
barotropic solve to run on fewer nodes or within a single node. We will
need to be able to create new environments based on a subset of an
existing environment. In setting up multiple environments, it will
be desireable to have some awareness of network topology for optimal
task placement.

### 2.8 Desired: Define threading parameters

If OpenMP threading is enabled for CPU, it may be useful to also define
similar master threads and thread layouts to enable some task parallelism
within an MPI rank.

### 2.9 Desired: Define other machine parameters

In the future, it may be useful to define other machine parameters, like
the specific configuration of cores, accelerators and other devices, to
manage task assignments on hybrid nodes. It should be easy to modify
this class to add additional information.

## 3 Algorithmic Formulation

No algorithms are needed beyond what is provided by standard MPI, OpenMP
libraries.

## 4 Design

In general, this is a simple class to hold information for later
retrieval.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

For improved vector performance on CPU and perhaps match thread block
sizes on GPU, we wish to block the inner loops with a compile-time
parameter. We define a CPP parameter: `OMEGA_VECTOR_SIZE`. This is
dependent on the machine and on whether GPU acceleration is enabled.
It will typically take on values like 16, 32, 64 for CPU-only builds
and 1 for GPU builds to maximize parallelism.

#### 4.1.2 Class/structs/data types

There will be a simple class MachEnv. We use a class here rather than
a struct so that we can make members private to prevent overwriting
these variables.

```c++
class MachEnv {

   private:
      int mComm;      ///< MPI communicator for this environment
      int mMyRank;    ///< rank ID for local MPI rank
      int mNumRanks;  ///< total number of MPI ranks
      int mMasterRank;///< rank ID for master rank
      bool mIsMaster; ///< true if the local rank is the master

   public:

      // Methods
      [define methods here - see below]
}
```

#### 4.1.3 Default environment

We will keep a default environment `OMEGA::defaultEnv` as a public
static instantiation that will be used by most of the OMEGA infrastructure.
If other environments are created, they must be maintained by the
defining routines or sub-components.

### 4.2 Methods

#### 4.2.1 Initialization

There will be two forms of the initialization routines. One of
these two must be called as early as possible in OMEGA initialization
(typically the first call). Both forms will return an integer
error code and will define the default environment `OMEGA::defaultEnv`.

```c++
// Initialization - standalone
int MachEnvInit();

// Initialization - coupled
int MachEnvInit(int inCommunicator, ///< [in] parent MPI communicator
                );
```

#### 4.2.2 Constructors

We provide several constructors for creating instances of the environment.
These will be used primarily by the above initialization routine, though
can also be used to create additional environments per requirement 2.7.

```c++
// Generic constructor that uses `MPI_COMM_WORLD`
MachEnv();

// Constructor that uses an assigned communicator (eg from coupler)
MachEnv(const int inCommunicator ///< [in] parent MPI communicator
        );

// Create a new environment from a contiguous subset of an
// existing environment
MachEnv(const int inCommunicator,///< [in] parent MPI communicator
        const int newSize        ///< [in] use first newSize ranks
        );

// Create a new environment from a strided subset of an
// existing environment
MachEnv(const int inCommunicator,///< [in] parent MPI communicator
        const int newSize,       ///< [in] num ranks in new env
        const int begin,         ///< [in] starting parent rank
        const int stride         ///< [in] stride for ranks to incl
        );

// Create a new environment from a custome subset of an
// existing environment, supplying list of parent ranks to include
MachEnv(const int inCommunicator ///< [in] parent MPI communicator
        const int newSize,       ///< [in] num ranks in new env
        const int ranks[]        ///< [in] vector of parent ranks to incl
        );
```

In standalone mode, the simple constructor would be used:

```c++
MachEnv omegaEnv(); // Initialize MPI and machine environment
```

and would define `MPI_COMM_WORLD` as the communicator as well
as initialize various other environments as needed.  In coupled mode,
the second form would be used with the coupler passing the
appropriate communicator to be used.

To satisfy requirement 2.7, we provide three forms to create a new
environment based on subsets of the old. The first simply creates
a subset from the first newSize ranks of the parent. The second
creates a strided subset (every "n" ranks starting from a specified
beginning rank). A third is the most general and creates a subset from
a vector of specific ranks of the parent to include in the new environment.

#### 4.2.2 Get/Retrieval

We provide specific retrieval functions for each class member to mimic
using this environment as a struct:

```c++
int Comm() const;      ///< returns MPI communicator for this environment
int MyRank() const;    ///< returns local MPI rank
int NumRanks() const;  ///< total number of MPI ranks
```

In typical use, these would look like:

```c++
myRank = OMEGA::defaultEnv.MyRank(); // retrieve local rank id
if (OMEGA::defaultEnv.IsMaster()){
   // do stuff on master rank
}
```

#### 4.2.3 Change master task

While most of the members should not be modified (read only using the get above), we will need a single set function to satisfy requirement 2.6.

```c++
int SetMaster(const int newMasterRank);
```

Note that this should occur as soon as possible after the environment is
created. Resetting the master after other activities have already assumed
the default master rank could leave variables defined on the wrong rank
and undefined on the new master rank. The integer return value is a
success/fail return code.


## 5 Verification and Testing

We will test this with a simple test driver in an 8-rank
MPI configuration. Requirement 2.4 (Kokkos initialization) will
need to be verified by visual inspection but later tests using
Kokkos functions will determine whether this has been successful.

### 5.1 Test standalone initialization

The test driver will call the standalone initialization, then
verify by retrieving the members and comparing to the equivalent
native MPI calls on `MPI_COMM_WORLD`
  * tests requirement 2.1, 2.3 and retrieval functions

### 5.2 Test changing master task

After the above test, use the set function to change the master
task to 1 and then verify using the retrieval function
  * tests requirement 2.6

### 5.3 Test multiple environments

Create three new environments using each of the three subset
constructors with the first creating a new environment using
the first 4 ranks, the second using every other rank (stride 2),
and the third using ranks 1,2,5,7. Verify that the members are
as expected when compared to the parent environment.
  * tests requirement 2.7

### 5.4 Initialization in coupled mode

While we can't test this directly in a standalone test driver,
we can test the underlying constructor by passing one of the
subset communicators as the input and verifying the resulting
environment is the same as the subset environment.
  * tests requirement 2.2 (mostly)

### 5.5 Test vector blocking factor

Within the test driver, set an int variable to `OMEGA_VECTOR_SIZE`,
and build the test with `-D OMEGA_VECTOR_SIZE=16`. Verify the
internal variable is also 16 to test that the preprocessor properly
propagates the value internally.
  * tests requirement 2.5
