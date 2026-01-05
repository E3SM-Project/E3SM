

# File emulator\_output\_stream.hpp



[**FileList**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_output\_stream.hpp**](emulator__output__stream_8hpp.md)

[Go to the source code of this file](emulator__output__stream_8hpp_source.md)

_Single output stream handler for diagnostic output._ [More...](#detailed-description)

* `#include "emulator_diagnostics.hpp"`
* `#include "emulator_logger.hpp"`
* `#include <functional>`
* `#include <map>`
* `#include <mpi.h>`
* `#include <string>`
* `#include <vector>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**emulator**](namespaceemulator.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**EmulatorOutputStream**](classemulator_1_1EmulatorOutputStream.md) <br>_Manages a single diagnostic output stream._  |
| class | [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) <br>_Interface for providing field data to output streams._  |
| struct | [**OutputControl**](structemulator_1_1OutputControl.md) <br>_Controls output timing and tracks averaging state._  |


















































## Detailed Description


Manages a single output stream including file creation, averaging, and snapshot writing at configured intervals. 


    

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_output_stream.hpp`

