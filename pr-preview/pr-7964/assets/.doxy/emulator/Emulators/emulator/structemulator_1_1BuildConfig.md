

# Struct emulator::BuildConfig



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**BuildConfig**](structemulator_1_1BuildConfig.md)



_Build-time configuration (set by buildnml, fixed for a case)._ [More...](#detailed-description)

* `#include <emulator_config.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::string | [**grid\_name**](#variable-grid_name)  <br>_Grid name (e.g., "gauss180x360")_  |
|  std::string | [**inference\_backend**](#variable-inference_backend)   = `"stub"`<br>_Backend: "stub" or "libtorch"._  |












































## Detailed Description


These values are determined at case creation time and cannot change between runs without rebuilding the case. 


    
## Public Attributes Documentation




### variable grid\_name 

_Grid name (e.g., "gauss180x360")_ 
```C++
std::string emulator::BuildConfig::grid_name;
```




<hr>



### variable inference\_backend 

_Backend: "stub" or "libtorch"._ 
```C++
std::string emulator::BuildConfig::inference_backend;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_config.hpp`

