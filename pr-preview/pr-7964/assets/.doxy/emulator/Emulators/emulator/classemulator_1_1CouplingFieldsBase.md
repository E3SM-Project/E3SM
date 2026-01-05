

# Class emulator::CouplingFieldsBase



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**CouplingFieldsBase**](classemulator_1_1CouplingFieldsBase.md)



_Base utility class for coupling field index management._ [More...](#detailed-description)

* `#include <coupling_fields.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  int | [**num\_exports**](#variable-num_exports)   = `0`<br>_Total number of export fields._  |
|  int | [**num\_imports**](#variable-num_imports)   = `0`<br>_Total number of import fields._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**get\_export\_index**](#function-get_export_index) (const std::string & name) const<br>_Look up the index of an export field by name._  |
|  int | [**get\_import\_index**](#function-get_import_index) (const std::string & name) const<br>_Look up the index of an import field by name._  |
| virtual void | [**initialize**](#function-initialize) (const std::string & export\_fields, const std::string & import\_fields) <br>_Initialize field indices from colon-separated field lists._  |
| virtual  | [**~CouplingFieldsBase**](#function-couplingfieldsbase) () = default<br> |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  std::map&lt; std::string, int &gt; | [**export\_map**](#variable-export_map)  <br>_Export field name to index._  |
|  std::map&lt; std::string, int &gt; | [**import\_map**](#variable-import_map)  <br>_Import field name to index._  |
















## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**parse\_field\_list**](#function-parse_field_list) (const std::string & fields, std::map&lt; std::string, int &gt; & field\_map, int & count) <br>_Parse a colon-separated field list into a name-&gt;index map._  |




## Detailed Description


## Public Attributes Documentation




### variable num\_exports 

_Total number of export fields._ 
```C++
int emulator::CouplingFieldsBase::num_exports;
```




<hr>



### variable num\_imports 

_Total number of import fields._ 
```C++
int emulator::CouplingFieldsBase::num_imports;
```




<hr>
## Public Functions Documentation




### function get\_export\_index 

_Look up the index of an export field by name._ 
```C++
inline int emulator::CouplingFieldsBase::get_export_index (
    const std::string & name
) const
```





**Parameters:**


* `name` Field name to look up 



**Returns:**

Field index (0-based), or -1 if not found 





        

<hr>



### function get\_import\_index 

_Look up the index of an import field by name._ 
```C++
inline int emulator::CouplingFieldsBase::get_import_index (
    const std::string & name
) const
```





**Parameters:**


* `name` Field name to look up 



**Returns:**

Field index (0-based), or -1 if not found 





        

<hr>



### function initialize 

_Initialize field indices from colon-separated field lists._ 
```C++
inline virtual void emulator::CouplingFieldsBase::initialize (
    const std::string & export_fields,
    const std::string & import_fields
) 
```



Parses the export and import field strings and builds internal lookup maps. Field strings use the MCT format: "field1:field2:field3".




**Parameters:**


* `export_fields` Colon-separated list of export (a2x) field names 
* `import_fields` Colon-separated list of import (x2a) field names 




        

<hr>



### function ~CouplingFieldsBase 

```C++
virtual emulator::CouplingFieldsBase::~CouplingFieldsBase () = default
```




<hr>
## Protected Attributes Documentation




### variable export\_map 

_Export field name to index._ 
```C++
std::map<std::string, int> emulator::CouplingFieldsBase::export_map;
```




<hr>



### variable import\_map 

_Import field name to index._ 
```C++
std::map<std::string, int> emulator::CouplingFieldsBase::import_map;
```




<hr>
## Protected Functions Documentation




### function parse\_field\_list 

_Parse a colon-separated field list into a name-&gt;index map._ 
```C++
inline void emulator::CouplingFieldsBase::parse_field_list (
    const std::string & fields,
    std::map< std::string, int > & field_map,
    int & count
) 
```





**Parameters:**


* `fields` Input string of colon-separated field names 
* `field_map` Output map from field name to index 
* `count` Output total count of fields parsed 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/coupling_fields.hpp`

