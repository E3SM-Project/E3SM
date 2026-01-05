

# Struct emulator::RuntimeConfig



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**RuntimeConfig**](structemulator_1_1RuntimeConfig.md)



_Runtime configuration (can change per run via atm\_in YAML)._ [More...](#detailed-description)

* `#include <emulator_config.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  bool | [**enabled**](#variable-enabled)   = `true`<br>_Enable inference (false = pass-through mode)_  |
|  std::string | [**ic\_file**](#variable-ic_file)  <br>_Initial conditions file path._  |
|  std::string | [**model\_path**](#variable-model_path)  <br>_Path to AI model file (TorchScript .pt)_  |












































## Detailed Description


These values can be modified between runs without rebuilding the case. 


    
## Public Attributes Documentation




### variable enabled 

_Enable inference (false = pass-through mode)_ 
```C++
bool emulator::RuntimeConfig::enabled;
```




<hr>



### variable ic\_file 

_Initial conditions file path._ 
```C++
std::string emulator::RuntimeConfig::ic_file;
```




<hr>



### variable model\_path 

_Path to AI model file (TorchScript .pt)_ 
```C++
std::string emulator::RuntimeConfig::model_path;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_config.hpp`

