

# Class emulator::DerivedDiagnostic



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**DerivedDiagnostic**](classemulator_1_1DerivedDiagnostic.md)



_Base class for derived diagnostics._ [More...](#detailed-description)

* `#include <derived_diagnostic.hpp>`





Inherited by the following classes: [emulator::HorizAvgDiagnostic](classemulator_1_1HorizAvgDiagnostic.md),  [emulator::VertSliceDiagnostic](classemulator_1_1VertSliceDiagnostic.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**compute**](#function-compute) (const [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) & fields, std::vector&lt; double &gt; & output) = 0<br>_Compute the diagnostic._  |
| virtual std::string | [**name**](#function-name) () const = 0<br>_Get the diagnostic name._  |
| virtual int | [**output\_size**](#function-output_size) (int ncols, int nlevs) const = 0<br>_Get output size._  |
| virtual std::string | [**source\_field**](#function-source_field) () const = 0<br>_Get the source field name._  |
| virtual  | [**~DerivedDiagnostic**](#function-deriveddiagnostic) () = default<br> |




























## Detailed Description


## Public Functions Documentation




### function compute 

_Compute the diagnostic._ 
```C++
virtual void emulator::DerivedDiagnostic::compute (
    const FieldDataProvider & fields,
    std::vector< double > & output
) = 0
```





**Parameters:**


* `fields` Input field provider 
* `output` Output data buffer (caller must allocate) 




        

<hr>



### function name 

_Get the diagnostic name._ 
```C++
virtual std::string emulator::DerivedDiagnostic::name () const = 0
```




<hr>



### function output\_size 

_Get output size._ 
```C++
virtual int emulator::DerivedDiagnostic::output_size (
    int ncols,
    int nlevs
) const = 0
```





**Parameters:**


* `ncols` Number of local columns 
* `nlevs` Number of vertical levels (for source field) 




        

<hr>



### function source\_field 

_Get the source field name._ 
```C++
virtual std::string emulator::DerivedDiagnostic::source_field () const = 0
```




<hr>



### function ~DerivedDiagnostic 

```C++
virtual emulator::DerivedDiagnostic::~DerivedDiagnostic () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/diagnostics/derived_diagnostic.hpp`

