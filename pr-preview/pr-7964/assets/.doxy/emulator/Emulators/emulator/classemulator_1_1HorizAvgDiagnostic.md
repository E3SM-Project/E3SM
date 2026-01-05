

# Class emulator::HorizAvgDiagnostic



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**HorizAvgDiagnostic**](classemulator_1_1HorizAvgDiagnostic.md)



_Computes area-weighted horizontal average of a field._ [More...](#detailed-description)

* `#include <horiz_avg_diagnostic.hpp>`



Inherits the following classes: [emulator::DerivedDiagnostic](classemulator_1_1DerivedDiagnostic.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**HorizAvgDiagnostic**](#function-horizavgdiagnostic) (const std::string & field\_name, const std::vector&lt; double &gt; & area\_weights, [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) comm) <br>_Construct horizontal averaging diagnostic._  |
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




### function HorizAvgDiagnostic 

_Construct horizontal averaging diagnostic._ 
```C++
emulator::HorizAvgDiagnostic::HorizAvgDiagnostic (
    const std::string & field_name,
    const std::vector< double > & area_weights,
    MPI_Comm comm
) 
```





**Parameters:**


* `field_name` Source field name 
* `area_weights` Area weights for each column (normalized to sum to 1) 
* `comm` MPI communicator for global reduction 




        

<hr>



### function compute 

_Compute the diagnostic._ 
```C++
virtual void emulator::HorizAvgDiagnostic::compute (
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
inline virtual std::string emulator::HorizAvgDiagnostic::name () override const
```



Implements [*emulator::DerivedDiagnostic::name*](classemulator_1_1DerivedDiagnostic.md#function-name)


<hr>



### function output\_size 

_Get output size._ 
```C++
inline virtual int emulator::HorizAvgDiagnostic::output_size (
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
inline virtual std::string emulator::HorizAvgDiagnostic::source_field () override const
```



Implements [*emulator::DerivedDiagnostic::source\_field*](classemulator_1_1DerivedDiagnostic.md#function-source_field)


<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/diagnostics/horiz_avg_diagnostic.hpp`

