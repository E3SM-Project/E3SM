

# File libtorch\_backend.cpp



[**FileList**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md) **>** [**libtorch\_backend.cpp**](libtorch__backend_8cpp.md)

[Go to the source code of this file](libtorch__backend_8cpp_source.md)

_LibTorch inference backend implementation._ [More...](#detailed-description)

* `#include "libtorch_backend.hpp"`
* `#include <iostream>`
* `#include <torch/script.h>`
* `#include <torch/torch.h>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**emulator**](namespaceemulator.md) <br> |
| namespace | [**inference**](namespaceemulator_1_1inference.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Impl**](structemulator_1_1inference_1_1LibTorchBackend_1_1Impl.md) <br>_Private implementation details for_ [_**LibTorchBackend**_](classemulator_1_1inference_1_1LibTorchBackend.md) _._ |


















































## Detailed Description


Provides native C++ neural network inference using LibTorch (PyTorch C++ API). This backend loads TorchScript models and executes inference without Python.




**Note:**

Models must be exported to TorchScript format (.pt) using torch.jit.trace() or torch.jit.script() before use with this backend.




**See also:** LibTorchBackend 



    

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/libtorch_backend.cpp`

