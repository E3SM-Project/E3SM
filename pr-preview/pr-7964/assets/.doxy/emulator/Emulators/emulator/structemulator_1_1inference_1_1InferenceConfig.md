

# Struct emulator::inference::InferenceConfig



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**inference**](namespaceemulator_1_1inference.md) **>** [**InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md)



_Configuration options for inference backends._ [More...](#detailed-description)

* `#include <inference_backend.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**BackendType**](namespaceemulator_1_1inference.md#enum-backendtype) | [**backend**](#variable-backend)   = `BackendType::STUB`<br>_Backend implementation to use._  |
|  int | [**device\_id**](#variable-device_id)   = `-1`<br>_GPU device ID (-1 for CPU)_  |
|  bool | [**dry\_run**](#variable-dry_run)   = `false`<br>_If true, validate config and exit without running._  |
|  std::vector&lt; std::string &gt; | [**expected\_input\_vars**](#variable-expected_input_vars)  <br>_Expected input variable names._  |
|  std::vector&lt; std::string &gt; | [**expected\_output\_vars**](#variable-expected_output_vars)  <br>_Expected output variable names._  |
|  int | [**grid\_height**](#variable-grid_height)   = `0`<br>_Height dimension (for spatial\_mode)_  |
|  int | [**grid\_width**](#variable-grid_width)   = `0`<br>_Width dimension (for spatial\_mode)_  |
|  int | [**input\_channels**](#variable-input_channels)   = `44`<br>_Number of input features per sample._  |
|  std::string | [**model\_path**](#variable-model_path)  <br>_Path to model file (TorchScript .pt for LibTorch)_  |
|  int | [**output\_channels**](#variable-output_channels)   = `50`<br>_Number of output features per sample._  |
|  bool | [**spatial\_mode**](#variable-spatial_mode)   = `/* multi line expression */`<br>_If true, reshape to [N, C, H, W] for CNN models._  |
|  bool | [**use\_fp16**](#variable-use_fp16)   = `false`<br>_Use half precision (requires CUDA)_  |
|  bool | [**verbose**](#variable-verbose)   = `false`<br>_Enable verbose output (for debugging)_  |












































## Detailed Description


Contains all parameters needed to initialize an inference backend, including model path, device selection, and tensor dimensions. 


    
## Public Attributes Documentation




### variable backend 

_Backend implementation to use._ 
```C++
BackendType emulator::inference::InferenceConfig::backend;
```




<hr>



### variable device\_id 

_GPU device ID (-1 for CPU)_ 
```C++
int emulator::inference::InferenceConfig::device_id;
```




<hr>



### variable dry\_run 

_If true, validate config and exit without running._ 
```C++
bool emulator::inference::InferenceConfig::dry_run;
```




<hr>



### variable expected\_input\_vars 

_Expected input variable names._ 
```C++
std::vector<std::string> emulator::inference::InferenceConfig::expected_input_vars;
```




<hr>



### variable expected\_output\_vars 

_Expected output variable names._ 
```C++
std::vector<std::string> emulator::inference::InferenceConfig::expected_output_vars;
```




<hr>



### variable grid\_height 

_Height dimension (for spatial\_mode)_ 
```C++
int emulator::inference::InferenceConfig::grid_height;
```




<hr>



### variable grid\_width 

_Width dimension (for spatial\_mode)_ 
```C++
int emulator::inference::InferenceConfig::grid_width;
```




<hr>



### variable input\_channels 

_Number of input features per sample._ 
```C++
int emulator::inference::InferenceConfig::input_channels;
```




<hr>



### variable model\_path 

_Path to model file (TorchScript .pt for LibTorch)_ 
```C++
std::string emulator::inference::InferenceConfig::model_path;
```




<hr>



### variable output\_channels 

_Number of output features per sample._ 
```C++
int emulator::inference::InferenceConfig::output_channels;
```




<hr>



### variable spatial\_mode 

_If true, reshape to [N, C, H, W] for CNN models._ 
```C++
bool emulator::inference::InferenceConfig::spatial_mode;
```




<hr>



### variable use\_fp16 

_Use half precision (requires CUDA)_ 
```C++
bool emulator::inference::InferenceConfig::use_fp16;
```




<hr>



### variable verbose 

_Enable verbose output (for debugging)_ 
```C++
bool emulator::inference::InferenceConfig::verbose;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/inference_backend.hpp`

