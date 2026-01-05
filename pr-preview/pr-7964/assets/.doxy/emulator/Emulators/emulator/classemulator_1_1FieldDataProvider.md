

# Class emulator::FieldDataProvider



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md)



_Interface for providing field data to output streams._ [More...](#detailed-description)

* `#include <emulator_output_stream.hpp>`





Inherited by the following classes: [emulator::impl::AtmFieldDataProvider](classemulator_1_1impl_1_1AtmFieldDataProvider.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual const std::vector&lt; double &gt; \* | [**get\_field**](#function-get_field) (const std::string & name) const = 0<br>_Get pointer to field data by name._  |
| virtual std::vector&lt; std::string &gt; | [**get\_field\_names**](#function-get_field_names) () const = 0<br>_Get list of all available field names._  |
| virtual int | [**get\_field\_nlevs**](#function-get_field_nlevs) (const std::string & name) const<br>_Check if a field has vertical levels (is 3D)._  |
| virtual int | [**get\_ncols**](#function-get_ncols) () const = 0<br>_Get number of local columns._  |
| virtual  | [**~FieldDataProvider**](#function-fielddataprovider) () = default<br> |




























## Detailed Description


Components implement this to provide access to their field data. 


    
## Public Functions Documentation




### function get\_field 

_Get pointer to field data by name._ 
```C++
virtual const std::vector< double > * emulator::FieldDataProvider::get_field (
    const std::string & name
) const = 0
```





**Parameters:**


* `name` Field name 



**Returns:**

Pointer to data vector, or nullptr if not found 





        

<hr>



### function get\_field\_names 

_Get list of all available field names._ 
```C++
virtual std::vector< std::string > emulator::FieldDataProvider::get_field_names () const = 0
```




<hr>



### function get\_field\_nlevs 

_Check if a field has vertical levels (is 3D)._ 
```C++
inline virtual int emulator::FieldDataProvider::get_field_nlevs (
    const std::string & name
) const
```





**Returns:**

Number of levels, or 1 for 2D fields 





        

<hr>



### function get\_ncols 

_Get number of local columns._ 
```C++
virtual int emulator::FieldDataProvider::get_ncols () const = 0
```




<hr>



### function ~FieldDataProvider 

```C++
virtual emulator::FieldDataProvider::~FieldDataProvider () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_output_stream.hpp`

