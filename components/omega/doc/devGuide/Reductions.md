(omega-dev-reductions)=

# Global Reductions

Global reductions are global collective operations to reduce distributed
data to a scalar value using either sum, min or max operators.


## Global sum

For scalars, `globalSum` computes a reproducible sum of all scalars
on MPI processors in a given MPI communicator:
```c++
int globalSum(const TT *Val, const MPI_Comm Comm, TT *Result)
```
where `TT` is the data type I4, I8, R4 or R8. Computed sum is stored
in `Result`. Any errors from the MPI collective call are passed back
in the function `int` return code.

For arrays, `globalSum` first computes a local sum of the array before
the global sum. The local sum is either between 0 and max indexes of
the array or between indexes defined in the specified `indexRange`:
```c++
int globalSum(const ArrayTTDD array,
              const MPI_Comm Comm,
              TT *Result,
              const std::vector<I4> *indexRange = nullptr)
```
`ArrayTTDD` is either a host or device OMEGA array of type I4, I8, R4 or R8
and dimension DD ranging from 1D to 5D. The `indexRange` vector is of length
2x the number of array dimensions: e.g. for a 2D array the min and
max indexes for a local sum might be
`indexRange{0, nCellsOwned-1, 0, nVertLevels-1}`. The index range vector
is an optional parameter and, when absent, local sum defaults to
min,max indexes of the array. Computed global sum is stored
in `Result`. Any errors from the MPI collective call are passed back
in the function `int` return code.


## Global sum with product

This function signature
```c++
int globalSum(const ArrayTTDD array1,
              const ArrayTTDD array2,
              const MPI_Comm Comm,
              TT *Result,
              const std::vector<I4> *indexRange = nullptr)
```
performs a sum of element-by-element product `array1*array2`. The second
array must be as big as the first array; if `indexRange` is provided, the
indexes must be valid for both arrays. This signature can be used to
mask the local sum of `array1` with values of masking `array2`.


## Global sum multi-field

To sum multiple fields at once, this function signature accepts a vector
of scalars:
```c++
int globalSum(const std::vector<TT> scalars,
              const MPI_Comm Comm,
              std::vector<TT> *Result)
```
or arrays:
```c++
int globalSum(const std::vector<ArrayTTDD> arrays,
              const MPI_Comm Comm,
              std::vector<TT> *Result,
              const std::vector<I4> *indexRange = nullptr)
```


## Global sum multi-field with product

This signature provides sums with masks for multiple fields:
```c++
int globalSum(const std::vector<ArrayTTDD> arrays1,
              const std::vector<ArrayTTDD> arrays2,
              const MPI_Comm Comm,
              std::vector<TT> *Result,,
              const std::vector<I4> *indexRange = nullptr)
```


## Global minval and maxval

Functions `globalMinVal` and `globalMaxVal` provide interfaces similar
to `globalSum` but reduce data with MIN and MAX operators:
```c++
int globalMinVal(const ArrayTTDD array,
                 const MPI_Comm Comm,
                 TT *Result,
                 const std::vector<I4> *indexRange = nullptr)

int globalMinVal(const ArrayTTDD array1,
                 const ArrayTTDD array2,
                 const MPI_Comm Comm,
                 TT *Result,
                 const std::vector<I4> *indexRange = nullptr)

int globalMinVal(const std::vector<ArrayTTDD> arrays,
                 const MPI_Comm Comm,
                 std::vector<TT> *Result,
                 const std::vector<I4> *indexRange = nullptr)
```


## Utility functions globalMin globalMax

These functions are utility functions to compute global MIN and MAX
across all MPI processors for a given scalar or array:
```c++
int globalMin(const TT *scalar,
              TT *Result,
              const MPI_Comm Comm)

int globalMin(const ArrayTTDD array,
              ArrayTTDD Result,
              const MPI_Comm Comm)

```
