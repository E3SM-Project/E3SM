

# Class emulator::inference::InferenceBackend



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**inference**](namespaceemulator_1_1inference.md) **>** [**InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md)



_Abstract interface for inference backends._ [More...](#detailed-description)

* `#include <inference_backend.hpp>`





Inherited by the following classes: [emulator::inference::LibTorchBackend](classemulator_1_1inference_1_1LibTorchBackend.md),  [emulator::inference::StubBackend](classemulator_1_1inference_1_1StubBackend.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**finalize**](#function-finalize) () = 0<br>_Release resources and finalize the backend._  |
| virtual bool | [**infer**](#function-infer) (const double \* inputs, double \* outputs, int batch\_size) = 0<br>_Run inference on input data._  |
| virtual bool | [**initialize**](#function-initialize) (const [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) & config) = 0<br>_Initialize the backend._  |
| virtual bool | [**is\_initialized**](#function-is_initialized) () const = 0<br>_Check if the backend is ready for inference._  |
| virtual std::string | [**name**](#function-name) () const = 0<br>_Get the human-readable name of this backend._  |
| virtual [**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) | [**type**](#function-type) () const = 0<br>_Get the backend type enumeration._  |
| virtual [**ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md) | [**validate**](#function-validate) () const<br>_Validate configuration before running._  |
| virtual  | [**~InferenceBackend**](#function-inferencebackend) () = default<br> |




























## Detailed Description


## Public Functions Documentation




### function finalize 

_Release resources and finalize the backend._ 
```C++
virtual void emulator::inference::InferenceBackend::finalize () = 0
```



After calling this, the backend is no longer usable until [**initialize()**](classemulator_1_1inference_1_1InferenceBackend.md#function-initialize) is called again. 


        

<hr>



### function infer 

_Run inference on input data._ 
```C++
virtual bool emulator::inference::InferenceBackend::infer (
    const double * inputs,
    double * outputs,
    int batch_size
) = 0
```



Executes the model on the provided input batch and writes results to the output buffer.




**Parameters:**


* `inputs` Input data array, size = batch\_size \* input\_channels 
* `outputs` Output data array, size = batch\_size \* output\_channels 
* `batch_size` Number of samples in the batch 



**Returns:**

true if inference succeeded, false on error




**Precondition:**

[**initialize()**](classemulator_1_1inference_1_1InferenceBackend.md#function-initialize) must have been called successfully 




**Precondition:**

outputs must be pre-allocated with sufficient size 





        

<hr>



### function initialize 

_Initialize the backend._ 
```C++
virtual bool emulator::inference::InferenceBackend::initialize (
    const InferenceConfig & config
) = 0
```



Loads the model, allocates resources, and prepares for inference. Must be called before [**infer()**](classemulator_1_1inference_1_1InferenceBackend.md#function-infer).




**Parameters:**


* `config` Configuration options 



**Returns:**

true if initialization succeeded, false on error 





        

<hr>



### function is\_initialized 

_Check if the backend is ready for inference._ 
```C++
virtual bool emulator::inference::InferenceBackend::is_initialized () const = 0
```





**Returns:**

true if initialized and ready 





        

<hr>



### function name 

_Get the human-readable name of this backend._ 
```C++
virtual std::string emulator::inference::InferenceBackend::name () const = 0
```





**Returns:**

Backend name (e.g., "LibTorch", "Stub") 





        

<hr>



### function type 

_Get the backend type enumeration._ 
```C++
virtual BackendType emulator::inference::InferenceBackend::type () const = 0
```





**Returns:**

[**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) value 





        

<hr>



### function validate 

_Validate configuration before running._ 
```C++
inline virtual ValidationResult emulator::inference::InferenceBackend::validate () const
```



Checks that the model file exists, dimensions match, device is available, etc. Call this after [**initialize()**](classemulator_1_1inference_1_1InferenceBackend.md#function-initialize) to detect configuration errors early.




**Returns:**

[**ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md) with errors/warnings if any 





        

<hr>



### function ~InferenceBackend 

```C++
virtual emulator::inference::InferenceBackend::~InferenceBackend () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/inference_backend.hpp`

