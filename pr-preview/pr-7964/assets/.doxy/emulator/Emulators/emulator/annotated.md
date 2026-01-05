
# Class List


Here are the classes, structs, unions and interfaces with brief descriptions:

* **namespace** [**anonymous namespace{components/emulator\_comps/eatm/src/emulator\_atm\_interface.cpp}**](namespace_0d070153036360144114261021302354156044166014006223.md) 
* **namespace** [**emulator**](namespaceemulator.md)     
    * **struct** [**BuildConfig**](structemulator_1_1BuildConfig.md) _Build-time configuration (set by buildnml, fixed for a case)._     
    * **struct** [**CouplingConfig**](structemulator_1_1CouplingConfig.md) _Coupling configuration options._     
    * **class** [**CouplingFieldsBase**](classemulator_1_1CouplingFieldsBase.md) _Base utility class for coupling field index management._     
    * **class** [**DataView**](classemulator_1_1DataView.md) _Non-owning view over contiguous data._     
    * **class** [**DerivedDiagnostic**](classemulator_1_1DerivedDiagnostic.md) _Base class for derived diagnostics._     
    * **struct** [**DiagnosticConfig**](structemulator_1_1DiagnosticConfig.md) _Complete diagnostic output configuration._     
    * **struct** [**DiagnosticMetadata**](structemulator_1_1DiagnosticMetadata.md) _Metadata for creating diagnostics._     
    * **class** [**EmulatorAtm**](classemulator_1_1EmulatorAtm.md) _Atmosphere emulator component._     
    * **class** [**EmulatorComp**](classemulator_1_1EmulatorComp.md) _Abstract base class for all emulated E3SM components._     
    * **struct** [**EmulatorConfig**](structemulator_1_1EmulatorConfig.md) _Base configuration structure for all emulator components._     
    * **class** [**EmulatorContext**](classemulator_1_1EmulatorContext.md) _Singleton context for managing emulator component instances._     
    * **class** [**EmulatorIO**](classemulator_1_1EmulatorIO.md) _Static class providing parallel I/O using SCORPIO/PIO._     
    * **class** [**EmulatorOutputManager**](classemulator_1_1EmulatorOutputManager.md) _Manages all diagnostic output for an emulator component._     
    * **class** [**EmulatorOutputStream**](classemulator_1_1EmulatorOutputStream.md) _Manages a single diagnostic output stream._     
    * **class** [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) _Interface for providing field data to output streams._     
    * **struct** [**HistoryRestartConfig**](structemulator_1_1HistoryRestartConfig.md) _History restart configuration for averaging buffers._     
    * **class** [**HorizAvgDiagnostic**](classemulator_1_1HorizAvgDiagnostic.md) _Computes area-weighted horizontal average of a field._     
    * **class** [**Logger**](classemulator_1_1Logger.md) _Simple logger with optional file output._     
    * **struct** [**ModelIOConfig**](structemulator_1_1ModelIOConfig.md) _Model I/O variable configuration._     
    * **struct** [**OutputControl**](structemulator_1_1OutputControl.md) _Controls output timing and tracks averaging state._     
    * **struct** [**OutputStreamConfig**](structemulator_1_1OutputStreamConfig.md) _Configuration for a single output stream._     
    * **struct** [**RestartConfig**](structemulator_1_1RestartConfig.md) _Restart output configuration._     
    * **struct** [**RuntimeConfig**](structemulator_1_1RuntimeConfig.md) _Runtime configuration (can change per run via atm\_in YAML)._     
    * **class** [**VertSliceDiagnostic**](classemulator_1_1VertSliceDiagnostic.md) _Extracts a single vertical level from a 3D field._     
    * **namespace** [**impl**](namespaceemulator_1_1impl.md)     
        * **struct** [**AtmCouplingIndices**](structemulator_1_1impl_1_1AtmCouplingIndices.md) _Coupling field indices for atmosphere component._     
        * **class** [**AtmFieldDataProvider**](classemulator_1_1impl_1_1AtmFieldDataProvider.md) _Adapter implementing_ [_**FieldDataProvider**_](classemulator_1_1FieldDataProvider.md) _for_[_**AtmFieldManager**_](classemulator_1_1impl_1_1AtmFieldManager.md) _._    
        * **class** [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) _Field storage container for atmosphere emulator._     
    * **namespace** [**inference**](namespaceemulator_1_1inference.md)     
        * **class** [**InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md) _Abstract interface for inference backends._     
        * **struct** [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) _Configuration options for inference backends._     
        * **class** [**LibTorchBackend**](classemulator_1_1inference_1_1LibTorchBackend.md) _LibTorch backend for native C++ PyTorch inference._     
            * **struct** [**Impl**](structemulator_1_1inference_1_1LibTorchBackend_1_1Impl.md) _Private implementation details for_ [_**LibTorchBackend**_](classemulator_1_1inference_1_1LibTorchBackend.md) _._    
        * **class** [**StubBackend**](classemulator_1_1inference_1_1StubBackend.md) _Stub backend for testing without actual inference._     
        * **struct** [**ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md) _Result of configuration validation._     
* **namespace** [**emulator**](namespaceemulator_1_1_0d123014276007031377224173102036175132222221352317.md) 
* **namespace** [**emulator**](namespaceemulator_1_1_0d351360357202143324273130142075163307310363145270.md) 

