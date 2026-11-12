(omega-dev-reductions)=

# Global Reductions

Global reductions are global collective operations to reduce distributed
data to a scalar value using either sum, min or max operators.

## Global sum

For scalars, `globalSum` computes a reproducible sum of all scalars
on MPI processors in a given MPI communicator:
```c++
TT Result = globalSum(const TT Val, const MPI_Comm Comm);
```
where `TT` is the data type I4, I8, R4 or R8. The computed sum across
MPI tasks is returned as the result.

For arrays, `globalSum` first computes a local sum of the array before
the sum across tasks of these local sums. The sum can be computed on a subset
of the array indices as defined by an optional vector argument `indexRange`:
```c++
TT Result = globalSum(const ArrayTTDD array,
                      const MPI_Comm Comm,
                      const std::vector<I4> *indexRange = nullptr);
```
`ArrayTTDD` is either a host or device Omega array of type I4, I8, R4 or R8
and dimension DD ranging from 1D to 5D. The `indexRange` vector is of length
2x the number of array dimensions: e.g. for a 2D array the min and
max indexes for a local sum might be
`indexRange{0, nCellsOwned-1, 0, nVertLayers-1}`. The index range vector
is an optional parameter and, when absent, local sum defaults to
min,max indexes of the array. The computed global sum is returned as
the `Result`.

## Global sum with product

This function signature
```c++
TT Result = globalSum(const ArrayTTDD Array1,
                      const ArrayTTDD Array2,
                      const MPI_Comm Comm,
                      const std::vector<I4> *indexRange = nullptr);
```
performs a sum of element-by-element product `array1*array2`. The second
array must be the same size and dimension of the first array; if `indexRange`
is provided, the indexes must be valid for both arrays. This signature can be
used to mask the local sum of `array1` with values of masking `array2`.

## Global sum multi-field

For scalars, multiple fields can be summed with a single call. A multi-field
interface for arrays is not yet supported.
This function signature accepts a vector of scalars:
```c++
std::vector<TT> Result = globalSum(const std::vector<TT> scalars,
                                   const MPI_Comm Comm);
```
The array interface (not yet supported) will be:
```c++
std::vector<TT> Result = globalSum(const std::vector<ArrayTTDD> arrays,
                                   const MPI_Comm Comm,
                                   const std::vector<I4> *indexRange = nullptr);
```
A sum-with-product for multiple fields is also planned but not yet supported.
Its signature will be:
```c++
std::vector<TT> Result = globalSum(const std::vector<ArrayTTDD> arrays1,
                                   const std::vector<ArrayTTDD> arrays2,
                                   const MPI_Comm Comm,
                                   const std::vector<I4> *indexRange = nullptr);
```

## Global minval and maxval

Functions `globalMinVal` and `globalMaxVal` provide interfaces analogous
to `globalSum` but will find the minimum or maximum value, respectively, of a
scalar or array across the global domain. The MinVal interfaces are:
```c++
TT Result = globalMinVal(const TT Val, const MPI_Comm Comm); // for a scalar

std::vector<TT> Result = globalMinVal(const std::vector<TT> Vals,
                                      const MPI_Comm Comm); // multiple scalars

TT Result = globalMinVal(const ArrayTTDD array,
                         const MPI_Comm Comm,
                         const std::vector<I4> *indexRange = nullptr);

TT Result = globalMinVal(const ArrayTTDD array1,
                         const ArrayTTDD array2,
                         const MPI_Comm Comm,
                         const std::vector<I4> *indexRange = nullptr);

TT Result = globalMinVal(const std::vector<ArrayTTDD> arrays,
                         const MPI_Comm Comm,
                         const std::vector<I4> *indexRange = nullptr)
```
with identical MaxVal interfaces to find the maximum. WARNING - the use of a
multiplicative mask to zero array entries works fine for global sums, but does
not always lead to desireable results in the MinVal and MaxVal interfaces.
Multiplying by a mask of 1s and 0s will result in 0 values in all masked
entries and the MinVal and MaxVal will view these as actual entries unless
limited by the indexRange argument. Care must be taken in the MinVal, MaxVal
case to mask appropriately.
