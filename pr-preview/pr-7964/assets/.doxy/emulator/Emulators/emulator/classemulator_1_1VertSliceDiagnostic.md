

# Class emulator::VertSliceDiagnostic



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**VertSliceDiagnostic**](classemulator_1_1VertSliceDiagnostic.md)



_Extracts a single vertical level from a 3D field._ [More...](#detailed-description)

* `#include <vert_slice_diagnostic.hpp>`



Inherits the following classes: [emulator::DerivedDiagnostic](classemulator_1_1DerivedDiagnostic.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**VertSliceDiagnostic**](#function-vertslicediagnostic) (const std::string & field\_name, int level\_idx, int nlevs) <br>_Construct vertical slicing diagnostic._  |
| virtual void | [**compute**](#function-compute) (const [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) & fields, std::vector&lt; double &gt; & output) override<br>_Compute the diagnostic._  |
| virtual std::string | [**name**](#function-name) () override const<br>_Get the diagnostic name._  |
| virtual int | [**output\_size**](#function-output_size) (int ncols, int nlevs) override const<br>_Get output size._  |
| virtual std::string | [**source\_field**](#function-source_field) () override const<br>_Get the source field name._  |


## Public Functions inherited from emulator::DerivedDiagnostic

See [emulator::DerivedDiagnostic](classemulator_1_1DerivedDiagnostic.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**compute**](classemulator_1_1DerivedDiagnostic.md#function-compute) (const [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) & fields, std::vector&lt; double &gt; & output) = 0<br>_Compute the diagnostic._  |
| virtual std::string | [**name**](classemulator_1_1DerivedDiagnostic.md#function-name) () const = 0<br>_Get the diagnostic name._  |
| virtual int | [**output\_size**](classemulator_1_1DerivedDiagnostic.md#function-output_size) (int ncols, int nlevs) const = 0<br>_Get output size._  |
| virtual std::string | [**source\_field**](classemulator_1_1DerivedDiagnostic.md#function-source_field) () const = 0<br>_Get the source field name._  |
| virtual  | [**~DerivedDiagnostic**](classemulator_1_1DerivedDiagnostic.md#function-deriveddiagnostic) () = default<br> |






















































## Detailed Description


## Public Functions Documentation




### function VertSliceDiagnostic 

_Construct vertical slicing diagnostic._ 
```C++
emulator::VertSliceDiagnostic::VertSliceDiagnostic (
    const std::string & field_name,
    int level_idx,
    int nlevs
) 
```





**Parameters:**


* `field_name` Source field name 
* `level_idx` Level index to extract (0-based) 
* `nlevs` Total number of levels in source field 




        

<hr>



### function compute 

_Compute the diagnostic._ 
```C++
virtual void emulator::VertSliceDiagnostic::compute (
    const FieldDataProvider & fields,
    std::vector< double > & output
) override
```





**Parameters:**


* `fields` Input field provider 
* `output` Output data buffer (caller must allocate) 




        
Implements [*emulator::DerivedDiagnostic::compute*](classemulator_1_1DerivedDiagnostic.md#function-compute)


<hr>



### function name 

_Get the diagnostic name._ 
```C++
inline virtual std::string emulator::VertSliceDiagnostic::name () override const
```



Implements [*emulator::DerivedDiagnostic::name*](classemulator_1_1DerivedDiagnostic.md#function-name)


<hr>



### function output\_size 

_Get output size._ 
```C++
inline virtual int emulator::VertSliceDiagnostic::output_size (
    int ncols,
    int nlevs
) override const
```





**Parameters:**


* `ncols` Number of local columns 
* `nlevs` Number of vertical levels (for source field) 




        
Implements [*emulator::DerivedDiagnostic::output\_size*](classemulator_1_1DerivedDiagnostic.md#function-output_size)


<hr>



### function source\_field 

_Get the source field name._ 
```C++
inline virtual std::string emulator::VertSliceDiagnostic::source_field () override const
```



Implements [*emulator::DerivedDiagnostic::source\_field*](classemulator_1_1DerivedDiagnostic.md#function-source_field)


<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/diagnostics/vert_slice_diagnostic.hpp`

