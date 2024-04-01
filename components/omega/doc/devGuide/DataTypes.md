(omega-dev-data-types)=

## Data Types and Precision

OMEGA supports all standard data types and uses some additional defined
types to guarantee a specific level of precision. In particular, we
define I4, I8, R4 and R8 types for 4-byte (32-bit) and 8-byte (64-bit)
integer and floating point variables. Note that these exist in the
OMEGA namespace so use the scoped form OMEGA::I4, etc.
In addition, we define a Real data type that is, by default,
double precision (8 bytes/64-bit) but if the code is built with a
`-DSINGLE_PRECISION` (see insert link to build system) preprocessor flag,
the default Real becomes single precision (4-byte/32-bit). For floating
point variables, developers should use this Real type instead of the
specific R4 or R8 forms unless the specific form is required
(eg converting to single precision before output or if a particular
algorithm is known to require double precision). This allows us to
easily convert all reals to single precision to explore performance or
accuracy in single precision mode. In some cases creating literal values
of type Real is necessary to avoid unwanted promotions. For that purpose
a user-defined literal `_Real` is provided. As an example, we can compute
the inverse area of a cell using only the Real type as follows:
```c++
    Real InvAreaCell = 1._Real / AreaCell(ICell);
```

## Arrays and Kokkos

The C++ language does not have native support for multi-dimensional
arrays as part of the language standard, though there are a number
of implementations as part of the Standard Template Library and
elsewhere. OMEGA uses the [Kokkos](https://github.com/kokkos)
framework for defining and allocating arrays on both the CPU host and
any accelerator device that may be present. Because the syntax for
defining such arrays is somewhat long, we instead define a number of
alias array types of the form `ArrayNDTT` where N is the dimension of
the array and TT is the data type (I4, I8, R4, R8 or Real) corresponding
to the types described above. The dimension refers to the number of
ranks in the array and not the physical dimension. Although Kokkos
supports Fortran ordering, we will use C ordering for array indices.
Within OMEGA the default location for an Array should be on the device
with a similar type HostArrayNDTT defined for arrays needed on the host.
As an example, we can define and allocate a device and host array using:
```c++
   Array3dReal Temperature("Temperature",nTimeLevels, nCells, nVertLevels);
   HostArray3dReal TemperatureHost("Temperature",nTimeLevels, nCells, nVertLevels);
```
Alternatively, you can use the copy functions to create a host copy
from the device or vice versa.
```c++
   auto TemperatureHost = OMEGA::createHostCopy(Temperature);
```
Finally, the arrays can be deallocated explicity using the class
deallocate method, eg `Temperature.deallocate();` or if they are local
to a routine, they will be automatically deallocated when they fall out
of scope on exit. More details on Kokkos arrays are available in the Kokkos
documentation.
