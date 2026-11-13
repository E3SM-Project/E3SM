(omega-design-broadcast)=
# Broadcast

## 1 Overview

In a parallel environment, it is often useful to broadcast values from
one MPI Rank to the others. In particular, the current design for
OMEGA configuration requires us to broadcast many aspects of that
configuration from the master rank to all others.

## 2 Requirements

### 2.1 Requirement: Broadcast all supported data types

We must be able to broadcast variables of I4, I8, R4, R8, Real,
boolean, std::string.

### 2.2 Requirement: Broadcast scalars and vectors

Generally we will be broadcasting scalar values, but will require
vectors of values in some cases.

### 2.3 Requirement: Broadcast from any rank

The most common use case will broadcast from the master rank but
we must be able to support broadcasting from any rank.

### 2.4 Desired: Broadcast within alternative environments

OMEGA plans to support subset partitions of the domain for sub-models,
so we must be able to broadcast only within these sub-groups using
the appropriate communicator.

### 2.5 Desired: Non-blocking broadcast

It is desireable with any communication function to be able to
overlap communication and computation using non-blocking forms.

## 3 Algorithmic Formulation

No algorithms needed beyond the MPI broadcast internal algorithm.

## 4 Design

The OMEGA broadcast is mostly just wrappers around the MPI broadcast
function with simplified arguments.

### 4.1 Data types and parameters

For non-blocking broadcasts, we will alias the `MPI_request` type
to a Broadcast ID:

```c++
using BroadcastID = MPI_Request;
```

### 4.2 Methods

We will have four forms of the Broadcast function, each returning an
integer error code.  Note that in the interfaces below, we could
template based on data type, but prefer to overload to keep the
interfaces cleaner.

#### 4.2.1 Broadcast from master task in default environment

The first form is the simplest for a broadcast within the default
environment and from the master task:

```c++
int Broadcast([data type] Value);
```

where `[data type]` is one of the supported types (I4, I8 R4, R8, Real,
boolean, std::string). In actual use, this would look like:

```c++
int Err = Broadcast(MyIntValue);
int Err = Broadcast(MyRealValue);
[etc for all data types]
```

On the master task, the value would be broadcast and remain unchanged.
The remaining tasks would receive the broadcast and store the value.

#### 4.2.2 Broadcast from another rank within default environment

This is similar to the above, but adds the additional argument
for the source rank to broadcast from.

```c++
int Broadcast([data type] Value,  ///< [in] value to be broadcast
              const int SrcRank   ///< [in] rank to broadcast from
              );      //
```

#### 4.2.3 Broadcast from master rank within a different environment

As in 4.2.1, but adds the machine environment as an argument:

```c++
int Broadcast([data type] Value,     ///< [in] value to be broadcast
              const MachEnv *SubEnv, ///< [in] defined OMEGA environment
              );
```

#### 4.2.4 Broadcast from another rank within a different environment

As in 4.2.2, but adds the machine environment as an argument:

```c++
int Broadcast([data type] Value,     ///< [in] value to be broadcast
              const MachEnv *SubEnv, ///< [in] defined OMEGA environment
              const int SrcRank      ///< [in] rank to broadcast from
              );
```

#### 4.2.5 Broadcast of vector variables

There will be interfaces identical to the above with the value argument
replaced by `std::vector<type> value` to broadcast a vector of values.

#### 4.2.6 Non-blocking broadcasts

Non-blocking forms of the above for all options will exist that
return a request id. Following the MPI standard, this will all
be named IBroadcast. In addition, an IBroadcastWait will be
included to wait for the request to complete. A non-blocking
sequence would look like:

```c++
BroadcastID myReqID = IBroadcast(MyVar);
[ perform other tasks/computation ]
int Err = IBroadcastWait(MyReqID);
```

## 5 Verification and Testing

### 5.1 Test scalars of all data types

A multi-processor test driver will initialize variables of all supported
types to zero or empty values. The master rank will then set a non-zero
value and broadcast to the remaining ranks. Each rank will verify the
new value is the expected one.
  - tests requirement 2.1

### 5.2 Test vectors of all types

Repeat the above with vectors of all types using a set of non-zero
values.
  - tests requirement 2.2

### 5.3 Test broadcast from other ranks

Repeat the above two tests except broadcasting from a rank other
than the master rank.
  - tests requirement 2.3

### 5.4 Test non-blocking broadcasts

Initialize non-blocking broadcasts in tandem with all the
above. After verification for the blocking forms, issue the
waits to complete the non-blocking forms and then perform the
verification on the non-blocking results.

### 5.5 Test within a subset environment

Create a new MachEnv with a subset of ranks and repeat the above tests
with the extra environment argument. Check also that ranks not in
the environment are not corrupted/overwritten with the broadcast value.
