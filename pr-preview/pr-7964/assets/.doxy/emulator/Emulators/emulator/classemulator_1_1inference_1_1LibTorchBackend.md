

# Class emulator::inference::LibTorchBackend



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**inference**](namespaceemulator_1_1inference.md) **>** [**LibTorchBackend**](classemulator_1_1inference_1_1LibTorchBackend.md)



_LibTorch backend for native C++ PyTorch inference._ [More...](#detailed-description)

* `#include <libtorch_backend.hpp>`



Inherits the following classes: [emulator::inference::InferenceBackend](classemulator_1_1inference_1_1InferenceBackend.md)












## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Impl**](structemulator_1_1inference_1_1LibTorchBackend_1_1Impl.md) <br>_Private implementation details for_ [_**LibTorchBackend**_](classemulator_1_1inference_1_1LibTorchBackend.md) _._ |










































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LibTorchBackend**](#function-libtorchbackend) () <br> |
| virtual void | [**finalize**](#function-finalize) () override<br>_Release resources and finalize the backend._  |
|  size\_t | [**get\_memory\_usage\_bytes**](#function-get_memory_usage_bytes) () const<br>_Get approximate memory usage in bytes._  |
| virtual bool | [**infer**](#function-infer) (const double \* inputs, double \* outputs, int batch\_size) override<br>_Run inference on input data._  |
| virtual bool | [**initialize**](#function-initialize) (const [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) & config) override<br>_Initialize the backend._  |
| virtual bool | [**is\_initialized**](#function-is_initialized) () override const<br>_Check if the backend is ready for inference._  |
| virtual std::string | [**name**](#function-name) () override const<br>_Get the human-readable name of this backend._  |
| virtual [**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) | [**type**](#function-type) () override const<br>_Get the backend type enumeration._  |
|   | [**~LibTorchBackend**](#function-libtorchbackend) () override<br> |


## Public Functions inherited from emulator::inference::InferenceBackend

See [emulator::inference::InferenceBackend](classemulator_1_1inference_1_1InferenceBackend.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**finalize**](classemulator_1_1inference_1_1InferenceBackend.md#function-finalize) () = 0<br>_Release resources and finalize the backend._  |
| virtual bool | [**infer**](classemulator_1_1inference_1_1InferenceBackend.md#function-infer) (const double \* inputs, double \* outputs, int batch\_size) = 0<br>_Run inference on input data._  |
| virtual bool | [**initialize**](classemulator_1_1inference_1_1InferenceBackend.md#function-initialize) (const [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) & config) = 0<br>_Initialize the backend._  |
| virtual bool | [**is\_initialized**](classemulator_1_1inference_1_1InferenceBackend.md#function-is_initialized) () const = 0<br>_Check if the backend is ready for inference._  |
| virtual std::string | [**name**](classemulator_1_1inference_1_1InferenceBackend.md#function-name) () const = 0<br>_Get the human-readable name of this backend._  |
| virtual [**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) | [**type**](classemulator_1_1inference_1_1InferenceBackend.md#function-type) () const = 0<br>_Get the backend type enumeration._  |
| virtual [**ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md) | [**validate**](classemulator_1_1inference_1_1InferenceBackend.md#function-validate) () const<br>_Validate configuration before running._  |
| virtual  | [**~InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md#function-inferencebackend) () = default<br> |






















































## Detailed Description


## Public Functions Documentation




### function LibTorchBackend 

```C++
emulator::inference::LibTorchBackend::LibTorchBackend () 
```




<hr>



### function finalize 

_Release resources and finalize the backend._ 
```C++
virtual void emulator::inference::LibTorchBackend::finalize () override
```



After calling this, the backend is no longer usable until [**initialize()**](classemulator_1_1inference_1_1LibTorchBackend.md#function-initialize) is called again.    


        
Implements [*emulator::inference::InferenceBackend::finalize*](classemulator_1_1inference_1_1InferenceBackend.md#function-finalize)


<hr>



### function get\_memory\_usage\_bytes 

_Get approximate memory usage in bytes._ 
```C++
size_t emulator::inference::LibTorchBackend::get_memory_usage_bytes () const
```





**Returns:**

Estimated memory used by model and buffers 





        

<hr>



### function infer 

_Run inference on input data._ 
```C++
virtual bool emulator::inference::LibTorchBackend::infer (
    const double * inputs,
    double * outputs,
    int batch_size
) override
```



Run inference on input data.


Executes the model on the provided input batch and writes results to the output buffer.




**Parameters:**


* `inputs` Input data array, size = batch\_size \* input\_channels 
* `outputs` Output data array, size = batch\_size \* output\_channels 
* `batch_size` Number of samples in the batch 



**Returns:**

true if inference succeeded, false on error




**Precondition:**

[**initialize()**](classemulator_1_1inference_1_1LibTorchBackend.md#function-initialize) must have been called successfully 




**Precondition:**

outputs must be pre-allocated with sufficient size   


Executes the loaded model on the provided input data. The backend expects input in [batch\_size, input\_channels] format. For CNN models requiring [N, C, H, W], the caller ([**EmulatorComp**](classemulator_1_1EmulatorComp.md)) must reshape first.




**Parameters:**


* `inputs` Input data array of size [batch\_size \* input\_channels] 
* `outputs` Output data array of size [batch\_size \* output\_channels] 
* `batch_size` Number of samples in the batch 



**Returns:**

true if inference succeeded, false on error 





        
Implements [*emulator::inference::InferenceBackend::infer*](classemulator_1_1inference_1_1InferenceBackend.md#function-infer)


<hr>



### function initialize 

_Initialize the backend._ 
```C++
virtual bool emulator::inference::LibTorchBackend::initialize (
    const InferenceConfig & config
) override
```



Initialize the LibTorch backend.


Loads the model, allocates resources, and prepares for inference. Must be called before [**infer()**](classemulator_1_1inference_1_1LibTorchBackend.md#function-infer).




**Parameters:**


* `config` Configuration options 



**Returns:**

true if initialization succeeded, false on error   


Loads the TorchScript model from the configured path and sets up the execution device and precision. 


        
Implements [*emulator::inference::InferenceBackend::initialize*](classemulator_1_1inference_1_1InferenceBackend.md#function-initialize)


<hr>



### function is\_initialized 

_Check if the backend is ready for inference._ 
```C++
inline virtual bool emulator::inference::LibTorchBackend::is_initialized () override const
```





**Returns:**

true if initialized and ready    





        
Implements [*emulator::inference::InferenceBackend::is\_initialized*](classemulator_1_1inference_1_1InferenceBackend.md#function-is_initialized)


<hr>



### function name 

_Get the human-readable name of this backend._ 
```C++
inline virtual std::string emulator::inference::LibTorchBackend::name () override const
```





**Returns:**

Backend name (e.g., "LibTorch", "Stub")    





        
Implements [*emulator::inference::InferenceBackend::name*](classemulator_1_1inference_1_1InferenceBackend.md#function-name)


<hr>



### function type 

_Get the backend type enumeration._ 
```C++
inline virtual BackendType emulator::inference::LibTorchBackend::type () override const
```





**Returns:**

[**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) value    





        
Implements [*emulator::inference::InferenceBackend::type*](classemulator_1_1inference_1_1InferenceBackend.md#function-type)


<hr>



### function ~LibTorchBackend 

```C++
emulator::inference::LibTorchBackend::~LibTorchBackend () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/libtorch_backend.hpp`

