

# File data\_view.hpp



[**FileList**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**data\_view.hpp**](data__view_8hpp.md)

[Go to the source code of this file](data__view_8hpp_source.md)

_Non-owning view abstraction for zero-copy data access._ [More...](#detailed-description)

* `#include <cstddef>`
* `#include <vector>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**emulator**](namespaceemulator.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**DataView**](classemulator_1_1DataView.md) &lt;typename T&gt;<br>_Non-owning view over contiguous data._  |


















































## Detailed Description


Provides a lightweight, non-owning view over contiguous data that can wrap raw pointers, std::vector, or (in the future) Kokkos views without copying data. 


    

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/data_view.hpp`

