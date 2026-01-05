

# Struct emulator::ModelIOConfig



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**ModelIOConfig**](structemulator_1_1ModelIOConfig.md)



_Model I/O variable configuration._ [More...](#detailed-description)

* `#include <emulator_config.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; std::string &gt; | [**input\_variables**](#variable-input_variables)  <br>_Fields to pack as model input._  |
|  std::vector&lt; std::string &gt; | [**output\_variables**](#variable-output_variables)  <br>_Fields to extract from output._  |
|  bool | [**spatial\_mode**](#variable-spatial_mode)   = `true`<br>_Enable spatial mode for CNN-based models (e.g., ACE2)._  |












































## Detailed Description


Specifies which fields are used as inputs and outputs for the AI model, and controls the data layout during inference. 


    
## Public Attributes Documentation




### variable input\_variables 

_Fields to pack as model input._ 
```C++
std::vector<std::string> emulator::ModelIOConfig::input_variables;
```




<hr>



### variable output\_variables 

_Fields to extract from output._ 
```C++
std::vector<std::string> emulator::ModelIOConfig::output_variables;
```




<hr>



### variable spatial\_mode 

_Enable spatial mode for CNN-based models (e.g., ACE2)._ 
```C++
bool emulator::ModelIOConfig::spatial_mode;
```



When true (CNN mode):



* Input data is reshaped from [H\*W, C] to [1, C, H, W] before inference
* Output data is reshaped from [1, C, H, W] to [H\*W, C] after inference
* The entire grid is processed as a single sample (batch\_size=1)




When false (pointwise/MLP mode):



* Data remains in [batch\_size, channels] format
* Each grid point is processed independently (batch\_size=H\*W)






**Note:**

Default: true (assumes CNN model like ACE2) 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_config.hpp`

