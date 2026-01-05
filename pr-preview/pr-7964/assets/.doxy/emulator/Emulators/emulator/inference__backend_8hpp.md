

# File inference\_backend.hpp



[**FileList**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md) **>** [**inference\_backend.hpp**](inference__backend_8hpp.md)

[Go to the source code of this file](inference__backend_8hpp_source.md)

_Abstract interface for neural network inference backends._ [More...](#detailed-description)

* `#include <memory>`
* `#include <string>`
* `#include <vector>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**emulator**](namespaceemulator.md) <br> |
| namespace | [**inference**](namespaceemulator_1_1inference.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md) <br>_Abstract interface for inference backends._  |
| struct | [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) <br>_Configuration options for inference backends._  |
| struct | [**ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md) <br>_Result of configuration validation._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef int | [**MPI\_Comm**](#typedef-mpi_comm)  <br> |















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**MPI\_COMM\_NULL**](inference__backend_8hpp.md#define-mpi_comm_null)  `0`<br> |
| define  | [**MPI\_COMM\_WORLD**](inference__backend_8hpp.md#define-mpi_comm_world)  `0`<br> |

## Detailed Description


Defines the common interface for all inference backends, allowing pluggable AI/ML frameworks for model execution. 


    
## Public Types Documentation




### typedef MPI\_Comm 

```C++
typedef int MPI_Comm;
```




<hr>
## Macro Definition Documentation





### define MPI\_COMM\_NULL 

```C++
#define MPI_COMM_NULL `0`
```




<hr>



### define MPI\_COMM\_WORLD 

```C++
#define MPI_COMM_WORLD `0`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/inference_backend.hpp`

