

# Struct emulator::inference::ValidationResult



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**inference**](namespaceemulator_1_1inference.md) **>** [**ValidationResult**](structemulator_1_1inference_1_1ValidationResult.md)



_Result of configuration validation._ [More...](#detailed-description)

* `#include <inference_backend.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; std::string &gt; | [**errors**](#variable-errors)  <br>_List of error messages._  |
|  bool | [**valid**](#variable-valid)   = `true`<br>_True if all checks passed._  |
|  std::vector&lt; std::string &gt; | [**warnings**](#variable-warnings)  <br>_List of warning messages._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**add\_error**](#function-add_error) (const std::string & msg) <br>_Add an error and mark as invalid._  |
|  void | [**add\_warning**](#function-add_warning) (const std::string & msg) <br>_Add a warning (does not affect validity)_  |
|  bool | [**has\_warnings**](#function-has_warnings) () const<br>_Check if there are any warnings._  |




























## Detailed Description


Returned by [**InferenceBackend::validate()**](classemulator_1_1inference_1_1InferenceBackend.md#function-validate) to report validation status and any errors found. 


    
## Public Attributes Documentation




### variable errors 

_List of error messages._ 
```C++
std::vector<std::string> emulator::inference::ValidationResult::errors;
```




<hr>



### variable valid 

_True if all checks passed._ 
```C++
bool emulator::inference::ValidationResult::valid;
```




<hr>



### variable warnings 

_List of warning messages._ 
```C++
std::vector<std::string> emulator::inference::ValidationResult::warnings;
```




<hr>
## Public Functions Documentation




### function add\_error 

_Add an error and mark as invalid._ 
```C++
inline void emulator::inference::ValidationResult::add_error (
    const std::string & msg
) 
```




<hr>



### function add\_warning 

_Add a warning (does not affect validity)_ 
```C++
inline void emulator::inference::ValidationResult::add_warning (
    const std::string & msg
) 
```




<hr>



### function has\_warnings 

_Check if there are any warnings._ 
```C++
inline bool emulator::inference::ValidationResult::has_warnings () const
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/inference_backend.hpp`

