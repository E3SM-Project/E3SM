

# Namespace emulator::impl



[**Namespace List**](namespaces.md) **>** [**emulator**](namespaceemulator.md) **>** [**impl**](namespaceemulator_1_1impl.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**AtmCouplingIndices**](structemulator_1_1impl_1_1AtmCouplingIndices.md) <br>_Coupling field indices for atmosphere component._  |
| class | [**AtmFieldDataProvider**](classemulator_1_1impl_1_1AtmFieldDataProvider.md) <br>_Adapter implementing_ [_**FieldDataProvider**_](classemulator_1_1FieldDataProvider.md) _for_[_**AtmFieldManager**_](classemulator_1_1impl_1_1AtmFieldManager.md) _._ |
| class | [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) <br>_Field storage container for atmosphere emulator._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**export\_atm\_fields**](#function-export_atm_fields) (double \* export\_data, int ncols, int nfields, const [**AtmCouplingIndices**](structemulator_1_1impl_1_1AtmCouplingIndices.md) & idx, const [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) & fields) <br>_Export fields from field manager to MCT buffer._  |
|  void | [**import\_atm\_fields**](#function-import_atm_fields) (const double \* import\_data, int ncols, int nfields, const [**AtmCouplingIndices**](structemulator_1_1impl_1_1AtmCouplingIndices.md) & idx, [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) & fields) <br>_Import fields from MCT buffer to field manager._  |
|  bool | [**read\_atm\_initial\_conditions**](#function-read_atm_initial_conditions) (const std::string & filename, int num\_global\_cols, int num\_local\_cols, const std::vector&lt; int &gt; & col\_gids, const std::vector&lt; double &gt; & lat, [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) & fields, const std::vector&lt; std::string &gt; & required\_vars, [**Logger**](classemulator_1_1Logger.md) & logger, bool is\_root) <br>_Read initial conditions from a NetCDF file._  |




























## Public Functions Documentation




### function export\_atm\_fields 

_Export fields from field manager to MCT buffer._ 
```C++
void emulator::impl::export_atm_fields (
    double * export_data,
    int ncols,
    int nfields,
    const AtmCouplingIndices & idx,
    const AtmFieldManager & fields
) 
```



Export fields from local field manager to MCT buffer.


Copies data from the [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md) vectors to the MCT export buffer (a2x fields).




**Parameters:**


* `export_data` Pointer to MCT export buffer (column-major layout) 
* `ncols` Number of local columns 
* `nfields` Number of export fields 
* `idx` Coupling indices 
* `fields` Field manager with data to export 




        

<hr>



### function import\_atm\_fields 

_Import fields from MCT buffer to field manager._ 
```C++
void emulator::impl::import_atm_fields (
    const double * import_data,
    int ncols,
    int nfields,
    const AtmCouplingIndices & idx,
    AtmFieldManager & fields
) 
```



Import fields from MCT buffer to local field manager.


MCT aVect layout is column-major (Fortran): rAttr(nfields, lsize) In C this translates to: data[col \* nfields + field\_idx]


Copies data from the MCT import buffer (x2a fields) to the corresponding vectors in the [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md).




**Parameters:**


* `import_data` Pointer to MCT import buffer (column-major layout) 
* `ncols` Number of local columns 
* `nfields` Number of import fields 
* `idx` Coupling indices 
* `fields` Field manager to populate



**Note:**

MCT uses column-major (Fortran) layout: data[col \* nfields + field\_idx] 





        

<hr>



### function read\_atm\_initial\_conditions 

_Read initial conditions from a NetCDF file._ 
```C++
bool emulator::impl::read_atm_initial_conditions (
    const std::string & filename,
    int num_global_cols,
    int num_local_cols,
    const std::vector< int > & col_gids,
    const std::vector< double > & lat,
    AtmFieldManager & fields,
    const std::vector< std::string > & required_vars,
    Logger & logger,
    bool is_root
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/eatm/src/impl/atm_coupling.cpp`

