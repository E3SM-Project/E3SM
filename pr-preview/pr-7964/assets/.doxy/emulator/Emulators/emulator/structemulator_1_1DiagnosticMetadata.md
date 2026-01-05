

# Struct emulator::DiagnosticMetadata



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**DiagnosticMetadata**](structemulator_1_1DiagnosticMetadata.md)



_Metadata for creating diagnostics._ 

* `#include <diagnostic_factory.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; double &gt; | [**area\_weights**](#variable-area_weights)  <br>_Area weights for horiz averaging._  |
|  [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) | [**comm**](#variable-comm)   = `[**MPI\_COMM\_WORLD**](inference__backend_8hpp.md#define-mpi_comm_world)`<br>_MPI communicator._  |
|  int | [**nlevs**](#variable-nlevs)   = `1`<br>_Default number of levels._  |












































## Public Attributes Documentation




### variable area\_weights 

_Area weights for horiz averaging._ 
```C++
std::vector<double> emulator::DiagnosticMetadata::area_weights;
```




<hr>



### variable comm 

_MPI communicator._ 
```C++
MPI_Comm emulator::DiagnosticMetadata::comm;
```




<hr>



### variable nlevs 

_Default number of levels._ 
```C++
int emulator::DiagnosticMetadata::nlevs;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/diagnostics/diagnostic_factory.hpp`

