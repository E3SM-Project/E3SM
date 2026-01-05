

# Namespace emulator



[**Namespace List**](namespaces.md) **>** [**emulator**](namespaceemulator.md)


















## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**impl**](namespaceemulator_1_1impl.md) <br> |
| namespace | [**inference**](namespaceemulator_1_1inference.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BuildConfig**](structemulator_1_1BuildConfig.md) <br>_Build-time configuration (set by buildnml, fixed for a case)._  |
| struct | [**CouplingConfig**](structemulator_1_1CouplingConfig.md) <br>_Coupling configuration options._  |
| class | [**CouplingFieldsBase**](classemulator_1_1CouplingFieldsBase.md) <br>_Base utility class for coupling field index management._  |
| class | [**DataView**](classemulator_1_1DataView.md) &lt;typename T&gt;<br>_Non-owning view over contiguous data._  |
| class | [**DerivedDiagnostic**](classemulator_1_1DerivedDiagnostic.md) <br>_Base class for derived diagnostics._  |
| struct | [**DiagnosticConfig**](structemulator_1_1DiagnosticConfig.md) <br>_Complete diagnostic output configuration._  |
| struct | [**DiagnosticMetadata**](structemulator_1_1DiagnosticMetadata.md) <br>_Metadata for creating diagnostics._  |
| class | [**EmulatorAtm**](classemulator_1_1EmulatorAtm.md) <br>_Atmosphere emulator component._  |
| class | [**EmulatorComp**](classemulator_1_1EmulatorComp.md) <br>_Abstract base class for all emulated E3SM components._  |
| struct | [**EmulatorConfig**](structemulator_1_1EmulatorConfig.md) <br>_Base configuration structure for all emulator components._  |
| class | [**EmulatorContext**](classemulator_1_1EmulatorContext.md) <br>_Singleton context for managing emulator component instances._  |
| class | [**EmulatorIO**](classemulator_1_1EmulatorIO.md) <br>_Static class providing parallel I/O using SCORPIO/PIO._  |
| class | [**EmulatorOutputManager**](classemulator_1_1EmulatorOutputManager.md) <br>_Manages all diagnostic output for an emulator component._  |
| class | [**EmulatorOutputStream**](classemulator_1_1EmulatorOutputStream.md) <br>_Manages a single diagnostic output stream._  |
| class | [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) <br>_Interface for providing field data to output streams._  |
| struct | [**HistoryRestartConfig**](structemulator_1_1HistoryRestartConfig.md) <br>_History restart configuration for averaging buffers._  |
| class | [**HorizAvgDiagnostic**](classemulator_1_1HorizAvgDiagnostic.md) <br>_Computes area-weighted horizontal average of a field._  |
| class | [**Logger**](classemulator_1_1Logger.md) <br>_Simple logger with optional file output._  |
| struct | [**ModelIOConfig**](structemulator_1_1ModelIOConfig.md) <br>_Model I/O variable configuration._  |
| struct | [**OutputControl**](structemulator_1_1OutputControl.md) <br>_Controls output timing and tracks averaging state._  |
| struct | [**OutputStreamConfig**](structemulator_1_1OutputStreamConfig.md) <br>_Configuration for a single output stream._  |
| struct | [**RestartConfig**](structemulator_1_1RestartConfig.md) <br>_Restart output configuration._  |
| struct | [**RuntimeConfig**](structemulator_1_1RuntimeConfig.md) <br>_Runtime configuration (can change per run via atm\_in YAML)._  |
| class | [**VertSliceDiagnostic**](classemulator_1_1VertSliceDiagnostic.md) <br>_Extracts a single vertical level from a 3D field._  |


## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**CompType**](#enum-comptype)  <br>_Enumeration of component types in E3SM._  |
| enum  | [**DataLayout**](#enum-datalayout)  <br>_Memory layout for multi-dimensional data._  |
| typedef [**DataView**](classemulator_1_1DataView.md)&lt; double &gt; | [**DoubleView**](#typedef-doubleview)  <br>_Convenience type aliases._  |
| enum  | [**FileType**](#enum-filetype)  <br>_File type indicator for restart discovery._  |
| typedef [**DataView**](classemulator_1_1DataView.md)&lt; float &gt; | [**FloatView**](#typedef-floatview)  <br> |
| enum  | [**FrequencyUnit**](#enum-frequencyunit)  <br>_Output frequency units._  |
| enum  | [**LogLevel**](#enum-loglevel)  <br>_Log level enumeration._  |
| enum  | [**OutputAvgType**](#enum-outputavgtype)  <br>_Output averaging type._  |
| enum  | [**OutputPrecision**](#enum-outputprecision)  <br>_Output data precision._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  double | [**FILLVALUE**](#variable-fillvalue)   = `1.0e20`<br>_Default fill value for missing data._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  std::string | [**avg\_type\_to\_str**](#function-avg_type_to_str) ([**OutputAvgType**](namespaceemulator.md#enum-outputavgtype) t) <br>_Convert_ [_**OutputAvgType**_](namespaceemulator.md#enum-outputavgtype) _to string._ |
|  void | [**cleanup\_emulator\_context**](#function-cleanup_emulator_context) () <br>_Convenience function to clean up the global emulator context._  |
|  std::unique\_ptr&lt; [**DerivedDiagnostic**](classemulator_1_1DerivedDiagnostic.md) &gt; | [**create\_diagnostic**](#function-create_diagnostic) (const std::string & diag\_name, const [**DiagnosticMetadata**](structemulator_1_1DiagnosticMetadata.md) & metadata) <br>_Parse a field name and create appropriate diagnostic if needed._  |
|  std::string | [**file\_type\_suffix**](#function-file_type_suffix) ([**FileType**](namespaceemulator.md#enum-filetype) t) <br>_Get file type suffix string._  |
|  std::string | [**freq\_unit\_to\_str**](#function-freq_unit_to_str) ([**FrequencyUnit**](namespaceemulator.md#enum-frequencyunit) u) <br>_Convert_ [_**FrequencyUnit**_](namespaceemulator.md#enum-frequencyunit) _to string._ |
|  std::string | [**get\_base\_field\_name**](#function-get_base_field_name) (const std::string & diag\_name) <br>_Extract base field name from diagnostic name._  |
|  std::string | [**get\_comp\_name**](#function-get_comp_name) ([**CompType**](namespaceemulator.md#enum-comptype) type) <br>_Get the short name for a component type._  |
|  std::string | [**get\_export\_prefix**](#function-get_export_prefix) ([**CompType**](namespaceemulator.md#enum-comptype) type) <br>_Get the coupling export prefix for a component type._  |
|  std::string | [**get\_import\_prefix**](#function-get_import_prefix) ([**CompType**](namespaceemulator.md#enum-comptype) type) <br>_Get the coupling import prefix for a component type._  |
|  bool | [**is\_derived\_diagnostic**](#function-is_derived_diagnostic) (const std::string & name) <br>_Check if a field name is a derived diagnostic pattern._  |
|  [**EmulatorConfig**](structemulator_1_1EmulatorConfig.md) | [**parse\_emulator\_config**](#function-parse_emulator_config) (const std::string & yaml\_file, const std::string & section\_name) <br>_Parse configuration from a YAML file._  |
|  [**EmulatorConfig**](structemulator_1_1EmulatorConfig.md) | [**parse\_emulator\_config\_with\_defaults**](#function-parse_emulator_config_with_defaults) (const std::string & yaml\_file, const std::string & section\_name, bool verbose) <br>_Parse configuration with fallback to defaults._  |
|  std::string | [**precision\_to\_str**](#function-precision_to_str) ([**OutputPrecision**](namespaceemulator.md#enum-outputprecision) p) <br>_Convert_ [_**OutputPrecision**_](namespaceemulator.md#enum-outputprecision) _to string._ |
|  [**OutputAvgType**](namespaceemulator.md#enum-outputavgtype) | [**str\_to\_avg\_type**](#function-str_to_avg_type) (const std::string & s) <br>_Convert string to_ [_**OutputAvgType**_](namespaceemulator.md#enum-outputavgtype) _._ |
|  [**FrequencyUnit**](namespaceemulator.md#enum-frequencyunit) | [**str\_to\_freq\_unit**](#function-str_to_freq_unit) (const std::string & s) <br>_Convert string to_ [_**FrequencyUnit**_](namespaceemulator.md#enum-frequencyunit) _._ |
|  [**OutputPrecision**](namespaceemulator.md#enum-outputprecision) | [**str\_to\_precision**](#function-str_to_precision) (const std::string & s) <br>_Convert string to_ [_**OutputPrecision**_](namespaceemulator.md#enum-outputprecision) _._ |




























## Public Types Documentation




### enum CompType 

_Enumeration of component types in E3SM._ 
```C++
enum emulator::CompType {
    ATM = 0,
    OCN = 1,
    ICE = 2,
    LND = 3
};
```




<hr>



### enum DataLayout 

_Memory layout for multi-dimensional data._ 
```C++
enum emulator::DataLayout {
    ROW_MAJOR,
    COLUMN_MAJOR,
    UNKNOWN
};
```




<hr>



### typedef DoubleView 

_Convenience type aliases._ 
```C++
using emulator::DoubleView =  DataView<double>;
```




<hr>



### enum FileType 

_File type indicator for restart discovery._ 
```C++
enum emulator::FileType {
    HISTORY,
    RESTART,
    HISTORY_RESTART
};
```




<hr>



### typedef FloatView 

```C++
using emulator::FloatView =  DataView<float>;
```




<hr>



### enum FrequencyUnit 

_Output frequency units._ 
```C++
enum emulator::FrequencyUnit {
    NSTEPS,
    NSECS,
    NMINS,
    NHOURS,
    NDAYS,
    NMONTHS,
    NYEARS,
    NONE
};
```



Specifies the time unit for output frequency. Follows EAMxx conventions. 


        

<hr>



### enum LogLevel 

_Log level enumeration._ 
```C++
enum emulator::LogLevel {
    DEBUG,
    VERBOSE,
    INFO,
    WARNING,
    ERROR
};
```




<hr>



### enum OutputAvgType 

_Output averaging type._ 
```C++
enum emulator::OutputAvgType {
    INSTANT,
    AVERAGE,
    MIN,
    MAX,
    STD,
    SUM
};
```



Specifies how to combine multiple timesteps in the output window. 


        

<hr>



### enum OutputPrecision 

_Output data precision._ 
```C++
enum emulator::OutputPrecision {
    FLOAT32,
    FLOAT64
};
```




<hr>
## Public Attributes Documentation




### variable FILLVALUE 

_Default fill value for missing data._ 
```C++
double emulator::FILLVALUE;
```




<hr>
## Public Functions Documentation




### function avg\_type\_to\_str 

_Convert_ [_**OutputAvgType**_](namespaceemulator.md#enum-outputavgtype) _to string._
```C++
std::string emulator::avg_type_to_str (
    OutputAvgType t
) 
```




<hr>



### function cleanup\_emulator\_context 

_Convenience function to clean up the global emulator context._ 
```C++
inline void emulator::cleanup_emulator_context () 
```



Equivalent to [**EmulatorContext::singleton()**](classemulator_1_1EmulatorContext.md#function-singleton).clean\_up(). 


        

<hr>



### function create\_diagnostic 

_Parse a field name and create appropriate diagnostic if needed._ 
```C++
std::unique_ptr< DerivedDiagnostic > emulator::create_diagnostic (
    const std::string & diag_name,
    const DiagnosticMetadata & metadata
) 
```



Recognized patterns:
* "{field}\_horiz\_avg" or "{field}\_global\_mean" → [**HorizAvgDiagnostic**](classemulator_1_1HorizAvgDiagnostic.md)
* "{field}\_at\_lev{N}" → [**VertSliceDiagnostic**](classemulator_1_1VertSliceDiagnostic.md) at level N






**Parameters:**


* `diag_name` Requested diagnostic/field name 
* `metadata` Context for creating diagnostics 



**Returns:**

Unique pointer to diagnostic, or nullptr if name is not a diagnostic 





        

<hr>



### function file\_type\_suffix 

_Get file type suffix string._ 
```C++
std::string emulator::file_type_suffix (
    FileType t
) 
```





**Returns:**

".atm.h.", ".atm.r.", or ".atm.rh." 





        

<hr>



### function freq\_unit\_to\_str 

_Convert_ [_**FrequencyUnit**_](namespaceemulator.md#enum-frequencyunit) _to string._
```C++
std::string emulator::freq_unit_to_str (
    FrequencyUnit u
) 
```




<hr>



### function get\_base\_field\_name 

_Extract base field name from diagnostic name._ 
```C++
std::string emulator::get_base_field_name (
    const std::string & diag_name
) 
```



Examples:
* "T\_horiz\_avg" → "T"
* "T\_at\_lev3" → "T"






**Parameters:**


* `diag_name` Diagnostic name 



**Returns:**

Base field name 





        

<hr>



### function get\_comp\_name 

_Get the short name for a component type._ 
```C++
inline std::string emulator::get_comp_name (
    CompType type
) 
```





**Parameters:**


* `type` Component type 



**Returns:**

Short name string (e.g., "atm" for atmosphere) 





        

<hr>



### function get\_export\_prefix 

_Get the coupling export prefix for a component type._ 
```C++
inline std::string emulator::get_export_prefix (
    CompType type
) 
```





**Parameters:**


* `type` Component type 



**Returns:**

Prefix string (e.g., "a2x" for atmosphere) 





        

<hr>



### function get\_import\_prefix 

_Get the coupling import prefix for a component type._ 
```C++
inline std::string emulator::get_import_prefix (
    CompType type
) 
```





**Parameters:**


* `type` Component type 



**Returns:**

Prefix string (e.g., "x2a" for atmosphere) 





        

<hr>



### function is\_derived\_diagnostic 

_Check if a field name is a derived diagnostic pattern._ 
```C++
bool emulator::is_derived_diagnostic (
    const std::string & name
) 
```





**Parameters:**


* `name` Field name to check 



**Returns:**

true if name matches a diagnostic pattern 





        

<hr>



### function parse\_emulator\_config 

_Parse configuration from a YAML file._ 
```C++
EmulatorConfig emulator::parse_emulator_config (
    const std::string & yaml_file,
    const std::string & section_name
) 
```



Parse emulator configuration from a YAML file.


Reads the YAML file and extracts configuration for the specified component section (e.g., "eatm"). All subsections (build, runtime, model\_io, coupling) are parsed if present.


Reads the specified YAML file and extracts the configuration for the given component section.




**Parameters:**


* `yaml_file` Path to YAML configuration file 
* `section_name` Name of the component section (e.g., "eatm", "eocn") 



**Returns:**

Parsed configuration structure 




**Exception:**


* `std::runtime_error` if file cannot be read or section is missing 




        

<hr>



### function parse\_emulator\_config\_with\_defaults 

_Parse configuration with fallback to defaults._ 
```C++
EmulatorConfig emulator::parse_emulator_config_with_defaults (
    const std::string & yaml_file,
    const std::string & section_name,
    bool verbose
) 
```



Parse emulator configuration with fallback to defaults.


Gracefully handles missing or malformed config files by returning default configuration values.


Attempts to parse the YAML file. If the file doesn't exist or parsing fails, returns default configuration values.




**Parameters:**


* `yaml_file` Path to YAML configuration file 
* `section_name` Name of the component section (e.g., "eatm", "eocn") 
* `verbose` Print debug information if true 



**Returns:**

Parsed configuration structure (or defaults on failure) 





        

<hr>



### function precision\_to\_str 

_Convert_ [_**OutputPrecision**_](namespaceemulator.md#enum-outputprecision) _to string._
```C++
std::string emulator::precision_to_str (
    OutputPrecision p
) 
```




<hr>



### function str\_to\_avg\_type 

_Convert string to_ [_**OutputAvgType**_](namespaceemulator.md#enum-outputavgtype) _._
```C++
OutputAvgType emulator::str_to_avg_type (
    const std::string & s
) 
```





**Parameters:**


* `s` String like "instant", "average", "min", "max", "std", "sum" 



**Returns:**

Corresponding [**OutputAvgType**](namespaceemulator.md#enum-outputavgtype) 





        

<hr>



### function str\_to\_freq\_unit 

_Convert string to_ [_**FrequencyUnit**_](namespaceemulator.md#enum-frequencyunit) _._
```C++
FrequencyUnit emulator::str_to_freq_unit (
    const std::string & s
) 
```





**Parameters:**


* `s` String like "nsteps", "ndays", "nmonths", etc. 



**Returns:**

Corresponding [**FrequencyUnit**](namespaceemulator.md#enum-frequencyunit) (NONE if invalid) 





        

<hr>



### function str\_to\_precision 

_Convert string to_ [_**OutputPrecision**_](namespaceemulator.md#enum-outputprecision) _._
```C++
OutputPrecision emulator::str_to_precision (
    const std::string & s
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/coupling_fields.cpp`

