
# File List

Here is a list of all files with brief descriptions:


* **dir** [**components**](dir_409f97388efe006bc3438b95e9edef48.md)     
    * **dir** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md)     
        * **dir** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md)     
            * **dir** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md)     
                * **dir** [**diagnostics**](dir_8c59288332f7499e80b506ac6e5a3903.md)     
                    * **file** [**derived\_diagnostic.hpp**](derived__diagnostic_8hpp.md) _Base class for derived diagnostics._     
                    * **file** [**diagnostic\_factory.cpp**](diagnostic__factory_8cpp.md) _Implementation of diagnostic factory functions._     
                    * **file** [**diagnostic\_factory.hpp**](diagnostic__factory_8hpp.md) _Factory functions for creating diagnostics._     
                    * **file** [**horiz\_avg\_diagnostic.cpp**](horiz__avg__diagnostic_8cpp.md) _Implementation of horizontal averaging diagnostic._     
                    * **file** [**horiz\_avg\_diagnostic.hpp**](horiz__avg__diagnostic_8hpp.md) _Horizontal averaging diagnostic._     
                    * **file** [**vert\_slice\_diagnostic.cpp**](vert__slice__diagnostic_8cpp.md) _Implementation of vertical slicing diagnostic._     
                    * **file** [**vert\_slice\_diagnostic.hpp**](vert__slice__diagnostic_8hpp.md) _Vertical slicing diagnostic._     
                * **dir** [**inference**](dir_93d8073dc3188aa38589a4c97b36101e.md)     
                    * **file** [**inference\_backend.hpp**](inference__backend_8hpp.md) _Abstract interface for neural network inference backends._     
                    * **file** [**inference\_factory.cpp**](inference__factory_8cpp.md) _Factory functions for creating inference backends._     
                    * **file** [**libtorch\_backend.cpp**](libtorch__backend_8cpp.md) _LibTorch inference backend implementation._     
                    * **file** [**libtorch\_backend.hpp**](libtorch__backend_8hpp.md) _LibTorch inference backend for native C++ PyTorch inference._     
                    * **file** [**stub\_backend.cpp**](stub__backend_8cpp.md) _Stub inference backend implementation._     
                    * **file** [**stub\_backend.hpp**](stub__backend_8hpp.md) _Stub inference backend for testing without ML dependencies._     
                * **file** [**coupling\_fields.cpp**](coupling__fields_8cpp.md) _Implementation file for coupling field utilities._     
                * **file** [**coupling\_fields.hpp**](coupling__fields_8hpp.md) _Base class for managing coupling field indices between components._     
                * **file** [**data\_view.hpp**](data__view_8hpp.md) _Non-owning view abstraction for zero-copy data access._     
                * **file** [**emulator\_comp.cpp**](emulator__comp_8cpp.md) _Implementation of the EmulatorComp base class._     
                * **file** [**emulator\_comp.hpp**](emulator__comp_8hpp.md) _Abstract base class for all emulated E3SM components._     
                * **file** [**emulator\_config.cpp**](emulator__config_8cpp.md) _Implementation of configuration parsing functions._     
                * **file** [**emulator\_config.hpp**](emulator__config_8hpp.md) _Configuration structures for emulator components._     
                * **file** [**emulator\_context.hpp**](emulator__context_8hpp.md) _Singleton context for managing emulator component instances._     
                * **file** [**emulator\_diagnostics.cpp**](emulator__diagnostics_8cpp.md) _Implementation of diagnostic configuration utilities._     
                * **file** [**emulator\_diagnostics.hpp**](emulator__diagnostics_8hpp.md) _Configuration structures for diagnostic output._     
                * **file** [**emulator\_io.cpp**](emulator__io_8cpp.md) _Implementation of the EmulatorIO class for parallel I/O._     
                * **file** [**emulator\_io.hpp**](emulator__io_8hpp.md) _I/O wrapper using SCORPIO/PIO C interface._     
                * **file** [**emulator\_logger.cpp**](emulator__logger_8cpp.md) _Implementation of the Logger class._     
                * **file** [**emulator\_logger.hpp**](emulator__logger_8hpp.md) _Simple logging utility for emulator components._     
                * **file** [**emulator\_output\_manager.cpp**](emulator__output__manager_8cpp.md) _Implementation of EmulatorOutputManager._     
                * **file** [**emulator\_output\_manager.hpp**](emulator__output__manager_8hpp.md) _Manages all diagnostic output streams._     
                * **file** [**emulator\_output\_stream.cpp**](emulator__output__stream_8cpp.md) _Implementation of EmulatorOutputStream._     
                * **file** [**emulator\_output\_stream.hpp**](emulator__output__stream_8hpp.md) _Single output stream handler for diagnostic output._     
        * **dir** [**eatm**](dir_54689134e1a693092e83f56806593839.md)     
            * **dir** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md)     
                * **dir** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md)     
                    * **file** [**atm\_coupling.cpp**](atm__coupling_8cpp.md) _Implementation of atmosphere coupling functions._     
                    * **file** [**atm\_coupling.hpp**](atm__coupling_8hpp.md) _Atmosphere coupling field indices and transfer functions._     
                    * **file** [**atm\_field\_data\_provider.cpp**](atm__field__data__provider_8cpp.md) _Implementation of AtmFieldDataProvider._     
                    * **file** [**atm\_field\_data\_provider.hpp**](atm__field__data__provider_8hpp.md) _FieldDataProvider adapter for AtmFieldManager._     
                    * **file** [**atm\_field\_manager.cpp**](atm__field__manager_8cpp.md) _Implementation of the AtmFieldManager class._     
                    * **file** [**atm\_field\_manager.hpp**](atm__field__manager_8hpp.md) _Field storage container for atmosphere emulator._     
                    * **file** [**atm\_io.cpp**](atm__io_8cpp.md) _Implementation of atmosphere I/O functions._     
                    * **file** [**atm\_io.hpp**](atm__io_8hpp.md) _I/O functions for atmosphere emulator initial conditions._     
                * **file** [**emulator\_atm.cpp**](emulator__atm_8cpp.md) _Atmosphere emulator component implementation._     
                * **file** [**emulator\_atm.hpp**](emulator__atm_8hpp.md) _Atmosphere emulator component declaration._     
                * **file** [**emulator\_atm\_interface.cpp**](emulator__atm__interface_8cpp.md) _C interface for the atmosphere emulator (Fortran-callable)._     

