

# Struct emulator::EmulatorConfig



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**EmulatorConfig**](structemulator_1_1EmulatorConfig.md)



_Base configuration structure for all emulator components._ [More...](#detailed-description)

* `#include <emulator_config.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**BuildConfig**](structemulator_1_1BuildConfig.md) | [**build**](#variable-build)  <br>_Build-time configuration._  |
|  [**CouplingConfig**](structemulator_1_1CouplingConfig.md) | [**coupling**](#variable-coupling)  <br>_Coupling configuration._  |
|  [**DiagnosticConfig**](structemulator_1_1DiagnosticConfig.md) | [**diagnostics**](#variable-diagnostics)  <br>_Diagnostic output configuration._  |
|  std::string | [**grid\_file**](#variable-grid_file)  <br> |
|  [**ModelIOConfig**](structemulator_1_1ModelIOConfig.md) | [**model\_io**](#variable-model_io)  <br>_Model I/O configuration._  |
|  [**RuntimeConfig**](structemulator_1_1RuntimeConfig.md) | [**runtime**](#variable-runtime)  <br>_Runtime configuration._  |












































## Detailed Description


## Public Attributes Documentation




### variable build 

_Build-time configuration._ 
```C++
BuildConfig emulator::EmulatorConfig::build;
```




<hr>



### variable coupling 

_Coupling configuration._ 
```C++
CouplingConfig emulator::EmulatorConfig::coupling;
```




<hr>



### variable diagnostics 

_Diagnostic output configuration._ 
```C++
DiagnosticConfig emulator::EmulatorConfig::diagnostics;
```




<hr>



### variable grid\_file 

```C++
std::string emulator::EmulatorConfig::grid_file;
```





**Deprecated**

Use domain info instead 




        

<hr>



### variable model\_io 

_Model I/O configuration._ 
```C++
ModelIOConfig emulator::EmulatorConfig::model_io;
```




<hr>



### variable runtime 

_Runtime configuration._ 
```C++
RuntimeConfig emulator::EmulatorConfig::runtime;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_config.hpp`

