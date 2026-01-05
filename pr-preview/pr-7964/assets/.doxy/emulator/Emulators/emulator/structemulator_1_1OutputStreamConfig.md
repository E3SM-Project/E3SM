

# Struct emulator::OutputStreamConfig



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**OutputStreamConfig**](structemulator_1_1OutputStreamConfig.md)



_Configuration for a single output stream._ [More...](#detailed-description)

* `#include <emulator_diagnostics.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**OutputAvgType**](namespaceemulator.md#enum-outputavgtype) | [**avg\_type**](#variable-avg_type)   = `OutputAvgType::INSTANT`<br> |
|  std::vector&lt; std::string &gt; | [**fields**](#variable-fields)  <br>_Fields to output._  |
|  std::string | [**filename\_prefix**](#variable-filename_prefix)   = `"emulator"`<br>_Output filename prefix._  |
|  int | [**frequency**](#variable-frequency)   = `1`<br>_Output every N units._  |
|  [**FrequencyUnit**](namespaceemulator.md#enum-frequencyunit) | [**frequency\_unit**](#variable-frequency_unit)   = `FrequencyUnit::NDAYS`<br> |
|  int | [**max\_snapshots\_per\_file**](#variable-max_snapshots_per_file)   = `1`<br>_Snapshots before new file._  |
|  [**OutputPrecision**](namespaceemulator.md#enum-outputprecision) | [**precision**](#variable-precision)   = `OutputPrecision::FLOAT32`<br> |
|  std::string | [**stream\_name**](#variable-stream_name)   = `"h0"`<br>_Stream identifier._  |












































## Detailed Description


Each stream writes to its own set of NetCDF files with configurable output frequency, averaging, and field selection. 


    
## Public Attributes Documentation




### variable avg\_type 

```C++
OutputAvgType emulator::OutputStreamConfig::avg_type;
```




<hr>



### variable fields 

_Fields to output._ 
```C++
std::vector<std::string> emulator::OutputStreamConfig::fields;
```




<hr>



### variable filename\_prefix 

_Output filename prefix._ 
```C++
std::string emulator::OutputStreamConfig::filename_prefix;
```




<hr>



### variable frequency 

_Output every N units._ 
```C++
int emulator::OutputStreamConfig::frequency;
```




<hr>



### variable frequency\_unit 

```C++
FrequencyUnit emulator::OutputStreamConfig::frequency_unit;
```




<hr>



### variable max\_snapshots\_per\_file 

_Snapshots before new file._ 
```C++
int emulator::OutputStreamConfig::max_snapshots_per_file;
```




<hr>



### variable precision 

```C++
OutputPrecision emulator::OutputStreamConfig::precision;
```




<hr>



### variable stream\_name 

_Stream identifier._ 
```C++
std::string emulator::OutputStreamConfig::stream_name;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_diagnostics.hpp`

