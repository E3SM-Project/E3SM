

# Namespace emulator::inference



[**Namespace List**](namespaces.md) **>** [**emulator**](namespaceemulator.md) **>** [**inference**](namespaceemulator_1_1inference.md)




















## Classes

| Type | Name |
| ---: | :--- |
| class | [**InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md) <br>_Abstract interface for inference backends._  |
| struct | [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) <br>_Configuration options for inference backends._  |
| class | [**LibTorchBackend**](classemulator_1_1inference_1_1LibTorchBackend.md) <br>_LibTorch backend for native C++ PyTorch inference._  |
| class | [**StubBackend**](classemulator_1_1inference_1_1StubBackend.md) <br>_Stub backend for testing without actual inference._  |
| struct | [**ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md) <br>_Result of configuration validation._  |


## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**BackendType**](#enum-backendtype)  <br>_Enumeration of available inference backend types._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  std::string | [**backend\_type\_to\_string**](#function-backend_type_to_string) ([**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) type) <br>_Convert a_ [_**BackendType**_](namespaceemulator_1_1inference.md#enum-backendtype) _to its string representation._ |
|  std::unique\_ptr&lt; [**InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md) &gt; | [**create\_backend**](#function-create_backend) ([**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) type) <br>_Factory function to create an inference backend by type._  |
|  std::unique\_ptr&lt; [**InferenceBackend**](classemulator_1_1inference_1_1InferenceBackend.md) &gt; | [**create\_backend**](#function-create_backend) (const [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) & config) <br>_Factory function to create and initialize a backend._  |
|  [**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) | [**parse\_backend\_type**](#function-parse_backend_type) (const std::string & str) <br>_Parse a backend type from a string._  |




























## Public Types Documentation




### enum BackendType 

_Enumeration of available inference backend types._ 
```C++
enum emulator::inference::BackendType {
    STUB,
    LIBTORCH
};
```




<hr>
## Public Functions Documentation




### function backend\_type\_to\_string 

_Convert a_ [_**BackendType**_](namespaceemulator_1_1inference.md#enum-backendtype) _to its string representation._
```C++
inline std::string emulator::inference::backend_type_to_string (
    BackendType type
) 
```





**Parameters:**


* `type` The backend type 



**Returns:**

String name of the backend 





        

<hr>



### function create\_backend 

_Factory function to create an inference backend by type._ 
```C++
std::unique_ptr< InferenceBackend > emulator::inference::create_backend (
    BackendType type
) 
```



Create an inference backend by type.


Creates an uninitialized backend instance. Call initialize() on the returned object before use.




**Parameters:**


* `type` Backend type to create 



**Returns:**

Unique pointer to new backend instance, or nullptr on error


Creates an uninitialized backend instance. The caller must call initialize() on the returned object before using it.




**Parameters:**


* `type` Backend type to create 



**Returns:**

Unique pointer to new backend, or [**StubBackend**](classemulator_1_1inference_1_1StubBackend.md) if type unavailable




**Note:**

If a requested backend is not available at compile time (e.g., LIBTORCH without EMULATOR\_HAS\_LIBTORCH defined), falls back to Stub. 





        

<hr>



### function create\_backend 

_Factory function to create and initialize a backend._ 
```C++
std::unique_ptr< InferenceBackend > emulator::inference::create_backend (
    const InferenceConfig & config
) 
```



Create and initialize a backend from configuration.


Convenience function that creates a backend and calls initialize().




**Parameters:**


* `config` Configuration including backend type and init parameters 



**Returns:**

Initialized backend, or nullptr if creation/init failed


Convenience function that creates a backend and calls initialize(). If initialization fails, returns nullptr.




**Parameters:**


* `config` Configuration with backend type and parameters 



**Returns:**

Initialized backend, or nullptr on failure 





        

<hr>



### function parse\_backend\_type 

_Parse a backend type from a string._ 
```C++
inline BackendType emulator::inference::parse_backend_type (
    const std::string & str
) 
```



Supports case-insensitive matching and common aliases.




**Parameters:**


* `str` String to parse ("libtorch", "torch", "stub", etc.) 



**Returns:**

Corresponding [**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) (defaults to STUB if unrecognized) 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/inference_backend.hpp`

