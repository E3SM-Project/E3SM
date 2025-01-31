<!--- OMEGA global reductions requirements and design ------------------------->

(omega-design-global-reductions)=


# Global Reductions

## 1 Overview

In a parallel application like OMEGA, it is often necessary to perform
averages, sums or other global collective operations on the distributed
data. These reductions are described here.

## 2 Requirements

### 2.1 Requirement: Global sums

The ability to perform the sum of all elements of a distributed array
is required

### 2.2 Requirement: Global min/max

The ability to find the minimum or maximum value in a distributed
array is required.

### 2.3 Requirement: Data types

There must be reductions for all supported arithmetic types
(integers, floating point) and all array dimensions from
0-d (scalars) to 5-d.

### 2.4 Requirement: Sums with product

In many cases, a global sum is often performed on a product of two
arrays. For example, a sum can be performed with a multiplicative
mask or a sum may be a dot product of two vectors within a solver.
A global sum interface that multiplies two arrays before the sum
is required.

### 2.5 Requirement: Restricted address space

In order to compute reductions over only owned and active
cells/edges or for regional diagnostics, there must be an option
to restrict the address space over which the local portion of sums
are computed.

### 2.6 Requirement: Multiple fields/arrays

For performance reasons, it is often beneficial to compute reductions
for multiple arrays at the same time to reduce the number of
messages and associated message latency. An interface for computing
independent reductions for multiple arrays at once is needed. Note
that the size of the arrays and any associated options must be the
same for all arrays in a multi-field reduction. It would be possible
to relax this restriction in the future.

### 2.7 Requirement: Other environments/decompositions

Because OMEGA supports multiple machine environments (communicators
and processor sets) and decompositions, the reductions must be
callable from environments other than the default.

### 2.8 Requirement: Reproducibility

The results (esp sums) must be bit-for-bit reproducible on a given
machine with respect to processor and thread count

### 2.9 Requirement: Accelerators and location

When running on accelerated architectures, the ability to specify
whether the reduction occurs on host or device is needed.

### 2.10 Requirement: CPU threading

For reductions performed on the CPU, the operations must support
and MPI+OpenMP model while still meeting reproducibility requirements.

### 2.11 Desired: Other reduction operations

It will be desireable in the future to add other reduction operations
like minloc/maxloc.

## 3 Algorithmic Formulation

For integers, the sum is straightforward. For single-precision floating
point, we will perform sums in double precision and convert back to
single before returning. In each of the above cases, the local Kokkos
versions of sum, minval, maxval will be used for the local sum before
accumulating the MPI sum across ranks.

For double-precision sums, we will begin with
a version of the double-double (DDPDD) algorithm of Donald Knuth,
further improved by David A Bailey and outlined in He and Ding (2001
J. of Supercomputing, 18, 259). In the future, we will also support
the more robust version by Pat Worley in the E3SM reprosum routines.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

In the future, if we support multiple reproducible algorithms,
we will need an enum to define choice of algorithm.

#### 4.1.2 Class/structs/data types

There are no new classes or data types associated with this
functionality.

### 4.2 Methods

#### 4.2.1 Global sum

For each supported array type ArrayTTDD (where TT is the data
type I4, I8, R4, R8, Real and DD is the dimension 1D thru 5D)
or supported scalar data type, there will be a sum function
aliased to the globalSum interface:

```c++
    sum = globalSum(const ArrayTTDD array,
                    const I4 mpiComm,
                    const std::vector<I4> indxRange);
```

where mpiComm is the MPI communicator extracted from the
relevant MachEnv (eg defaultEnv.comm) and indxRange is a
`std::vector` of length 2x the number of array dimensions.
The entries of this vector will be the min, max index of
each array dimension. So, for example an array dimensioned
`(nCellsAll+1, nVertLevels+1)` might need to be summed over
the `indxRange{0, nCellsOwned-1, 0, nVertLevels-1}`. Note
that the indxRange supports a constant scalar values only so
does not support a variable index range like
minLevelCell/maxLevelCell. That is best managed either
by masking with the sum-product interface below or
ensuring the array has been set to zero in non-active
entries.

For scalar sums, the interface is the same with the
indxRange argument absent.

We will assume that Kokkos host arrays will be summed on the
host and device arrays will be summed on the device.

#### 4.2.2 Global sum with product

There will be a similar interface for a sum with product:

```c++
    sum = globalSum(const ArrayTTDD array1,
                    const ArrayTTDD array2,
                    const I4 mpiComm,
                    const std::vector<I4> indxRange);
```

that returns the sum of `array1*array2`. The two arrays
must be of the same data type. It is not required that
the two arrays be of exactly the same size but both arrays
must be of appropriate size to accomodate the indxRange
and the indices must align appropriately. This typically
means that if the sizes differ, any extra entries occur
at the end of the respective dimension.

#### 4.2.3 Global sum multi-field

To sum multiple fields at once, a `std::vector` of
arrays must be constructed and the result will be
stored in a `std::vector` of appropriate type:

```c++
    std::vector<type> sum =
         globalSum(const std::vector<ArrayTTDD> arrays,
                   const I4 mpiComm,
                   const std::vector<I4> indxRange)
```

As is the case with sum-product, each array must be
of the same data type but can have slightly different
sizes as long as the dimension lengths can accomodate
the indxRange and that the indices align appropriately
with the indxRange requested.

#### 4.2.4 Global sum multi-field with product

The multi-field sum with product is still to be determined, but
should have an interface like:

```c++
    std::vector<type> sum =
         globalSum(const std::vector<ArrayTTDD> arrays1,
                   const std::vector<ArrayTTDD> arrays2,
                   const I4 mpiComm,
                   const std::vector<I4> indxRange)
```

The implementation detail that needs to be determined is how
to manage two use cases. A typical use case is for the second
array in the product to be the same for all fields (eg a mask or
an area array), but we may wish to support a case where each
array product has a different array for both operands.
Because Kokkos arrays contain metadata and a pointer to
the data, it may be possible to create a std::vector of
the same array, allowing both the case of a fixed array
to be used for all or for all the array products to have
different arrays. This requires some testing of approaches.

#### 4.2.5 Global minval

The global minval will support interfaces similar to the
global sums, but will return the minimum value instead.

```c++
    varMin = globalMinval(const ArrayTTDD array,
                          const I4 mpiComm,
                          const std::vector<I4> indxRange);

    varMin = globalMinval(const ArrayTTDD array1,
                          const ArrayTTDD array2,
                          const I4 mpiComm,
                          const std::vector<I4> indxRange);

    std::vector<type> varMin =
         globalMinval(const std::vector<ArrayTTDD> arrays,
                      const I4 mpiComm,
                      const std::vector<I4> indxRange)
```

Note that the minval routine will not be cognizant of any
special values so any inactive entries should either be
masked or set to an appropriately high value.

#### 4.2.6 Global maxval

The maxval will have the identical form of minval but will
return the maximum value. The same restrictions and comments
on special values apply.

## 5 Verification and Testing

### 5.1 Test basics

For each data type, a set of reference arrays will be created
with random numbers that extend across the entire range of that
type (for later reproducibility tests). A reference sum will
be computed serially using the same reproducible algorithm
as the global implementation. The basic test will compare the
global sum with the reference serial sum. Because the data types
include both host and device arrays, this will also test reductions
on either location.
  - tests requirement 2.1, 2.3, part of 2.8, 2.9

### 5.2 Reproducibility

A MachEnv and decomposition for a subset of MPI ranks will be
created and the reference arrays above will be distributed across
the new decomposition. Global sums in each case will be compared
to the serial reference value and the full MPI rank case
for bit reproducibility. Similarly, a test with CPU threading
on will test reproducibility under threading.
  - tests requirement 2.7, 2.8, 2.10

### 5.3 Restricted index range

The sums will be tested with a reduced index range on the same
reference arrays from above and compared with the
similarly-restricted serial reference sums.
  - tests requirement 2.5

### 5.4 Sum with product

A mask array will be defined for each type and used with the
sum-with-product interface and compared against a serial
version of the same.
  - tests requirement 2.4

### 5.5 Multi-field sums

Additional arrays similar to the reference arrays will be created
with reference serial sums. Then the multi-field interface will
be tested against these reference sums.
  - tests requirement 2.6

### 5.6 Min/Max

A similar approach to the sums above will be used to test the
min and max functions in all combinations of interfaces.
  - tests requirement 2.2
