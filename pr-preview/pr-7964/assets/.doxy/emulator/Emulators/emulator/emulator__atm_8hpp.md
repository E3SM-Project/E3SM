

# File emulator\_atm.hpp



[**FileList**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**emulator\_atm.hpp**](emulator__atm_8hpp.md)

[Go to the source code of this file](emulator__atm_8hpp_source.md)

_Atmosphere emulator component declaration._ [More...](#detailed-description)

* `#include "../../common/src/coupling_fields.hpp"`
* `#include "../../common/src/emulator_comp.hpp"`
* `#include "../../common/src/emulator_config.hpp"`
* `#include "../../common/src/emulator_output_manager.hpp"`
* `#include "../../common/src/inference/inference_backend.hpp"`
* `#include "impl/atm_coupling.hpp"`
* `#include "impl/atm_field_data_provider.hpp"`
* `#include "impl/atm_field_manager.hpp"`
* `#include <fstream>`
* `#include <memory>`
* `#include <string>`
* `#include <vector>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**emulator**](namespaceemulator.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**EmulatorAtm**](classemulator_1_1EmulatorAtm.md) <br>_Atmosphere emulator component._  |


















































## Detailed Description


Defines the EmulatorAtm class, which implements an AI-based atmosphere component for E3SM. Supports multiple inference backends (STUB, LibTorch) for executing neural network models like ACE2. 


    

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/eatm/src/emulator_atm.hpp`

