

# Class emulator::impl::AtmFieldDataProvider



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**impl**](namespaceemulator_1_1impl.md) **>** [**AtmFieldDataProvider**](classemulator_1_1impl_1_1AtmFieldDataProvider.md)



_Adapter implementing_ [_**FieldDataProvider**_](classemulator_1_1FieldDataProvider.md) _for_[_**AtmFieldManager**_](classemulator_1_1impl_1_1AtmFieldManager.md) _._[More...](#detailed-description)

* `#include <atm_field_data_provider.hpp>`



Inherits the following classes: [emulator::FieldDataProvider](classemulator_1_1FieldDataProvider.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**AtmFieldDataProvider**](#function-atmfielddataprovider) ([**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) & fields, int ncols\_local) <br>_Construct adapter with reference to field manager._  |
|  void | [**detect\_stacked\_fields**](#function-detect_stacked_fields) () <br>_Scan fields and detect stackable slice patterns._  |
| virtual const std::vector&lt; double &gt; \* | [**get\_field**](#function-get_field) (const std::string & name) override const<br>_Get pointer to field data by name._  |
| virtual std::vector&lt; std::string &gt; | [**get\_field\_names**](#function-get_field_names) () override const<br>_Get list of all available field names._  |
| virtual int | [**get\_field\_nlevs**](#function-get_field_nlevs) (const std::string & name) override const<br>_Get number of vertical levels for a field._  |
| virtual int | [**get\_ncols**](#function-get_ncols) () override const<br>_Get number of local columns._  |
|  const std::vector&lt; double &gt; & | [**get\_stacked\_field**](#function-get_stacked_field) (const std::string & basename) const<br>_Get the stacked data for a 3D field._  |
|  bool | [**is\_stacked\_field**](#function-is_stacked_field) (const std::string & name) const<br>_Check if a field is a stacked field._  |
|   | [**~AtmFieldDataProvider**](#function-atmfielddataprovider) () override<br> |


## Public Functions inherited from emulator::FieldDataProvider

See [emulator::FieldDataProvider](classemulator_1_1FieldDataProvider.md)

| Type | Name |
| ---: | :--- |
| virtual const std::vector&lt; double &gt; \* | [**get\_field**](classemulator_1_1FieldDataProvider.md#function-get_field) (const std::string & name) const = 0<br>_Get pointer to field data by name._  |
| virtual std::vector&lt; std::string &gt; | [**get\_field\_names**](classemulator_1_1FieldDataProvider.md#function-get_field_names) () const = 0<br>_Get list of all available field names._  |
| virtual int | [**get\_field\_nlevs**](classemulator_1_1FieldDataProvider.md#function-get_field_nlevs) (const std::string & name) const<br>_Check if a field has vertical levels (is 3D)._  |
| virtual int | [**get\_ncols**](classemulator_1_1FieldDataProvider.md#function-get_ncols) () const = 0<br>_Get number of local columns._  |
| virtual  | [**~FieldDataProvider**](classemulator_1_1FieldDataProvider.md#function-fielddataprovider) () = default<br> |






















































## Detailed Description


## Public Functions Documentation




### function AtmFieldDataProvider 

_Construct adapter with reference to field manager._ 
```C++
emulator::impl::AtmFieldDataProvider::AtmFieldDataProvider (
    AtmFieldManager & fields,
    int ncols_local
) 
```





**Parameters:**


* `fields` Reference to [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) 
* `ncols_local` Number of local columns 




        

<hr>



### function detect\_stacked\_fields 

_Scan fields and detect stackable slice patterns._ 
```C++
void emulator::impl::AtmFieldDataProvider::detect_stacked_fields () 
```



Finds fields matching "basename\_N" pattern and registers them as stackable groups. 


        

<hr>



### function get\_field 

_Get pointer to field data by name._ 
```C++
virtual const std::vector< double > * emulator::impl::AtmFieldDataProvider::get_field (
    const std::string & name
) override const
```



First checks for exact match, then checks for stacked fields.




**Parameters:**


* `name` Field name 



**Returns:**

Pointer to data vector, or nullptr if not found 





        
Implements [*emulator::FieldDataProvider::get\_field*](classemulator_1_1FieldDataProvider.md#function-get_field)


<hr>



### function get\_field\_names 

_Get list of all available field names._ 
```C++
virtual std::vector< std::string > emulator::impl::AtmFieldDataProvider::get_field_names () override const
```



Includes both direct fields and auto-detected stacked fields. 


        
Implements [*emulator::FieldDataProvider::get\_field\_names*](classemulator_1_1FieldDataProvider.md#function-get_field_names)


<hr>



### function get\_field\_nlevs 

_Get number of vertical levels for a field._ 
```C++
virtual int emulator::impl::AtmFieldDataProvider::get_field_nlevs (
    const std::string & name
) override const
```





**Parameters:**


* `name` Field name 



**Returns:**

Number of levels (1 for 2D fields, &gt;1 for stacked fields) 





        
Implements [*emulator::FieldDataProvider::get\_field\_nlevs*](classemulator_1_1FieldDataProvider.md#function-get_field_nlevs)


<hr>



### function get\_ncols 

_Get number of local columns._ 
```C++
inline virtual int emulator::impl::AtmFieldDataProvider::get_ncols () override const
```



Implements [*emulator::FieldDataProvider::get\_ncols*](classemulator_1_1FieldDataProvider.md#function-get_ncols)


<hr>



### function get\_stacked\_field 

_Get the stacked data for a 3D field._ 
```C++
const std::vector< double > & emulator::impl::AtmFieldDataProvider::get_stacked_field (
    const std::string & basename
) const
```



Combines sliced fields [field\_0, field\_1, ...] into a single contiguous buffer with shape [nlevs, ncols].




**Parameters:**


* `basename` Base field name (without \_N suffix) 



**Returns:**

Stacked data buffer 





        

<hr>



### function is\_stacked\_field 

_Check if a field is a stacked field._ 
```C++
bool emulator::impl::AtmFieldDataProvider::is_stacked_field (
    const std::string & name
) const
```




<hr>



### function ~AtmFieldDataProvider 

```C++
emulator::impl::AtmFieldDataProvider::~AtmFieldDataProvider () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/eatm/src/impl/atm_field_data_provider.hpp`

