
# Class Hierarchy

This inheritance list is sorted roughly, but not completely, alphabetically:


* **class** [**emulator::CouplingFieldsBase**](classemulator_1_1CouplingFieldsBase.md) _Base utility class for coupling field index management._ 
* **class** [**emulator::DataView**](classemulator_1_1DataView.md) _Non-owning view over contiguous data._ 
* **class** [**emulator::DerivedDiagnostic**](classemulator_1_1DerivedDiagnostic.md) _Base class for derived diagnostics._     
    * **class** [**emulator::HorizAvgDiagnostic**](classemulator_1_1HorizAvgDiagnostic.md) _Computes area-weighted horizontal average of a field._ 
    * **class** [**emulator::VertSliceDiagnostic**](classemulator_1_1VertSliceDiagnostic.md) _Extracts a single vertical level from a 3D field._ 
* **class** [**emulator::EmulatorComp**](classemulator_1_1EmulatorComp.md) _Abstract base class for all emulated E3SM components._     
    * **class** [**emulator::EmulatorAtm**](classemulator_1_1EmulatorAtm.md) _Atmosphere emulator component._ 
* **class** [**emulator::EmulatorContext**](classemulator_1_1EmulatorContext.md) _Singleton context for managing emulator component instances._ 
* **class** [**emulator::EmulatorIO**](classemulator_1_1EmulatorIO.md) _Static class providing parallel I/O using SCORPIO/PIO._ 
* **class** [**emulator::EmulatorOutputManager**](classemulator_1_1EmulatorOutputManager.md) _Manages all diagnostic output for an emulator component._ 
* **class** [**emulator::EmulatorOutputStream**](classemulator_1_1EmulatorOutputStream.md) _Manages a single diagnostic output stream._ 
* **class** [**emulator::FieldDataProvider**](classemulator_1_1FieldDataProvider.md) _Interface for providing field data to output streams._     
    * **class** [**emulator::impl::AtmFieldDataProvider**](classemulator_1_1impl_1_1AtmFieldDataProvider.md) _Adapter implementing_ [_**FieldDataProvider**_](classemulator_1_1FieldDataProvider.md) _for_[_**AtmFieldManager**_](classemulator_1_1impl_1_1AtmFieldManager.md) _._
* **class** [**emulator::Logger**](classemulator_1_1Logger.md) _Simple logger with optional file output._ 
* **class** [**emulator::impl::AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) _Field storage container for atmosphere emulator._ 
* **class** [**emulator::inference::InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md) _Abstract interface for inference backends._     
    * **class** [**emulator::inference::LibTorchBackend**](classemulator_1_1inference_1_1LibTorchBackend.md) _LibTorch backend for native C++ PyTorch inference._ 
    * **class** [**emulator::inference::StubBackend**](classemulator_1_1inference_1_1StubBackend.md) _Stub backend for testing without actual inference._ 
* **struct** [**emulator::BuildConfig**](structemulator_1_1BuildConfig.md) _Build-time configuration (set by buildnml, fixed for a case)._ 
* **struct** [**emulator::CouplingConfig**](structemulator_1_1CouplingConfig.md) _Coupling configuration options._ 
* **struct** [**emulator::DiagnosticConfig**](structemulator_1_1DiagnosticConfig.md) _Complete diagnostic output configuration._ 
* **struct** [**emulator::DiagnosticMetadata**](structemulator_1_1DiagnosticMetadata.md) _Metadata for creating diagnostics._ 
* **struct** [**emulator::EmulatorConfig**](structemulator_1_1EmulatorConfig.md) _Base configuration structure for all emulator components._ 
* **struct** [**emulator::HistoryRestartConfig**](structemulator_1_1HistoryRestartConfig.md) _History restart configuration for averaging buffers._ 
* **struct** [**emulator::ModelIOConfig**](structemulator_1_1ModelIOConfig.md) _Model I/O variable configuration._ 
* **struct** [**emulator::OutputControl**](structemulator_1_1OutputControl.md) _Controls output timing and tracks averaging state._ 
* **struct** [**emulator::OutputStreamConfig**](structemulator_1_1OutputStreamConfig.md) _Configuration for a single output stream._ 
* **struct** [**emulator::RestartConfig**](structemulator_1_1RestartConfig.md) _Restart output configuration._ 
* **struct** [**emulator::RuntimeConfig**](structemulator_1_1RuntimeConfig.md) _Runtime configuration (can change per run via atm\_in YAML)._ 
* **struct** [**emulator::impl::AtmCouplingIndices**](structemulator_1_1impl_1_1AtmCouplingIndices.md) _Coupling field indices for atmosphere component._ 
* **struct** [**emulator::inference::InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) _Configuration options for inference backends._ 
* **struct** [**emulator::inference::LibTorchBackend::Impl**](structemulator_1_1inference_1_1LibTorchBackend_1_1Impl.md) _Private implementation details for_ [_**LibTorchBackend**_](classemulator_1_1inference_1_1LibTorchBackend.md) _._
* **struct** [**emulator::inference::ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md) _Result of configuration validation._ 

