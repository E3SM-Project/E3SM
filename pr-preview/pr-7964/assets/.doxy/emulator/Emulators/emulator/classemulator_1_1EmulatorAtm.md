

# Class emulator::EmulatorAtm



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**EmulatorAtm**](classemulator_1_1EmulatorAtm.md)



_Atmosphere emulator component._ [More...](#detailed-description)

* `#include <emulator_atm.hpp>`



Inherits the following classes: [emulator::EmulatorComp](classemulator_1_1EmulatorComp.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**EmulatorAtm**](#function-emulatoratm) () <br> |
|  const [**inference::InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) & | [**get\_inference\_config**](#function-get_inference_config) () const<br>_Get current inference configuration._  |
|  void | [**init\_coupling\_indices**](#function-init_coupling_indices) (const std::string & export\_fields, const std::string & import\_fields) <br>_Initialize coupling field indices from MCT field lists._  |
|  void | [**set\_inference\_config**](#function-set_inference_config) (const [**inference::InferenceConfig**](structemulator_1_1inference_1_1InferenceConfig.md) & config) <br>_Configure inference backend._  |
|  void | [**set\_log\_file**](#function-set_log_file) (const std::string & filename) <br>_Set logging file._  |
|   | [**~EmulatorAtm**](#function-emulatoratm) () override<br> |


## Public Functions inherited from emulator::EmulatorComp

See [emulator::EmulatorComp](classemulator_1_1EmulatorComp.md)

| Type | Name |
| ---: | :--- |
|   | [**EmulatorComp**](classemulator_1_1EmulatorComp.md#function-emulatorcomp) ([**CompType**](namespaceemulator.md#enum-comptype) type) <br>_Construct an emulator component of the given type._  |
|  [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) | [**comm**](classemulator_1_1EmulatorComp.md#function-comm) () const<br>_Get MPI communicator._  |
|  int | [**comp\_id**](classemulator_1_1EmulatorComp.md#function-comp_id) () const<br>_Get component ID._  |
|  void | [**create\_instance**](classemulator_1_1EmulatorComp.md#function-create_instance) ([**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) comm, int comp\_id, const char \* input\_file, int run\_type, int start\_ymd, int start\_tod) <br>_Initialize the component instance._  |
|  void | [**finalize**](classemulator_1_1EmulatorComp.md#function-finalize) () <br>_Finalize and clean up the component._  |
|  void | [**get\_cols\_area**](classemulator_1_1EmulatorComp.md#function-get_cols_area) (double \* area) const<br>_Get area for local columns._  |
|  void | [**get\_cols\_latlon**](classemulator_1_1EmulatorComp.md#function-get_cols_latlon) (double \* lat, double \* lon) const<br>_Get lat/lon for local columns._  |
|  void | [**get\_local\_col\_gids**](classemulator_1_1EmulatorComp.md#function-get_local_col_gids) (int \* gids) const<br>_Get global IDs for local columns._  |
|  int | [**get\_num\_global\_cols**](classemulator_1_1EmulatorComp.md#function-get_num_global_cols) () const<br>_Get total number of columns globally._  |
|  int | [**get\_num\_local\_cols**](classemulator_1_1EmulatorComp.md#function-get_num_local_cols) () const<br>_Get number of columns on this MPI rank._  |
|  int | [**get\_nx**](classemulator_1_1EmulatorComp.md#function-get_nx) () const<br>_Get number of grid points in x-direction (longitude)._  |
|  int | [**get\_ny**](classemulator_1_1EmulatorComp.md#function-get_ny) () const<br>_Get number of grid points in y-direction (latitude)._  |
|  void | [**initialize**](classemulator_1_1EmulatorComp.md#function-initialize) () <br>_Initialize the component (phase 2)._  |
|  bool | [**is\_root**](classemulator_1_1EmulatorComp.md#function-is_root) () const<br>_Check if this is the root rank (rank 0)._  |
|  int | [**rank**](classemulator_1_1EmulatorComp.md#function-rank) () const<br>_Get MPI rank within component communicator._  |
|  void | [**run**](classemulator_1_1EmulatorComp.md#function-run) (int dt) <br>_Execute one time step._  |
|  void | [**set\_grid\_data**](classemulator_1_1EmulatorComp.md#function-set_grid_data) (int nx, int ny, int num\_local\_cols, int num\_global\_cols, const int \* col\_gids, const double \* lat, const double \* lon, const double \* area) <br>_Set grid decomposition data from the driver._  |
|  void | [**setup\_coupling**](classemulator_1_1EmulatorComp.md#function-setup_coupling) (double \* import\_data, double \* export\_data, int num\_imports, int num\_exports, int field\_size) <br>_Set up coupling data buffers._  |
|  [**CompType**](namespaceemulator.md#enum-comptype) | [**type**](classemulator_1_1EmulatorComp.md#function-type) () const<br>_Get component type._  |
| virtual  | [**~EmulatorComp**](classemulator_1_1EmulatorComp.md#function-emulatorcomp) () = default<br> |
















## Protected Attributes inherited from emulator::EmulatorComp

See [emulator::EmulatorComp](classemulator_1_1EmulatorComp.md)

| Type | Name |
| ---: | :--- |
|  std::vector&lt; double &gt; | [**m\_area**](classemulator_1_1EmulatorComp.md#variable-m_area)  <br>_Cell area [km²]._  |
|  std::vector&lt; int &gt; | [**m\_col\_gids**](classemulator_1_1EmulatorComp.md#variable-m_col_gids)  <br>_Global IDs for local columns._  |
|  [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) | [**m\_comm**](classemulator_1_1EmulatorComp.md#variable-m_comm)  <br>_MPI communicator._  |
|  int | [**m\_comp\_id**](classemulator_1_1EmulatorComp.md#variable-m_comp_id)  <br>_Component ID from driver._  |
|  int | [**m\_current\_tod**](classemulator_1_1EmulatorComp.md#variable-m_current_tod)   = `0`<br>_Current time of day [seconds]._  |
|  int | [**m\_current\_ymd**](classemulator_1_1EmulatorComp.md#variable-m_current_ymd)   = `0`<br>_Current date as YYYYMMDD._  |
|  double \* | [**m\_export\_data**](classemulator_1_1EmulatorComp.md#variable-m_export_data)   = `nullptr`<br>_Export buffer pointer (a2x)_  |
|  int | [**m\_field\_size**](classemulator_1_1EmulatorComp.md#variable-m_field_size)   = `0`<br>_Size per field (should = num\_local\_cols)_  |
|  double \* | [**m\_import\_data**](classemulator_1_1EmulatorComp.md#variable-m_import_data)   = `nullptr`<br>_Import buffer pointer (x2a)_  |
|  std::string | [**m\_input\_file**](classemulator_1_1EmulatorComp.md#variable-m_input_file)  <br>_Path to configuration file._  |
|  std::vector&lt; double &gt; | [**m\_lat**](classemulator_1_1EmulatorComp.md#variable-m_lat)  <br>_Latitude [radians]._  |
|  [**Logger**](classemulator_1_1Logger.md) | [**m\_logger**](classemulator_1_1EmulatorComp.md#variable-m_logger)  <br>_Component logger._  |
|  std::vector&lt; double &gt; | [**m\_lon**](classemulator_1_1EmulatorComp.md#variable-m_lon)  <br>_Longitude [radians]._  |
|  int | [**m\_nprocs**](classemulator_1_1EmulatorComp.md#variable-m_nprocs)  <br>_Number of MPI processes._  |
|  int | [**m\_num\_exports**](classemulator_1_1EmulatorComp.md#variable-m_num_exports)   = `0`<br>_Number of export fields._  |
|  int | [**m\_num\_global\_cols**](classemulator_1_1EmulatorComp.md#variable-m_num_global_cols)   = `0`<br>_Global column count._  |
|  int | [**m\_num\_imports**](classemulator_1_1EmulatorComp.md#variable-m_num_imports)   = `0`<br>_Number of import fields._  |
|  int | [**m\_num\_local\_cols**](classemulator_1_1EmulatorComp.md#variable-m_num_local_cols)   = `0`<br>_Local column count._  |
|  int | [**m\_nx**](classemulator_1_1EmulatorComp.md#variable-m_nx)   = `0`<br>_Grid points in x (longitude)_  |
|  int | [**m\_ny**](classemulator_1_1EmulatorComp.md#variable-m_ny)   = `0`<br>_Grid points in y (latitude)_  |
|  int | [**m\_rank**](classemulator_1_1EmulatorComp.md#variable-m_rank)  <br>_MPI rank within component._  |
|  int | [**m\_run\_type**](classemulator_1_1EmulatorComp.md#variable-m_run_type)  <br>_Run type (startup/continue/branch)_  |
|  int | [**m\_step\_count**](classemulator_1_1EmulatorComp.md#variable-m_step_count)   = `0`<br>_Number of steps executed._  |
|  [**CompType**](namespaceemulator.md#enum-comptype) | [**m\_type**](classemulator_1_1EmulatorComp.md#variable-m_type)  <br>_Component type enum._  |






























## Protected Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**final\_impl**](#function-final_impl) () override<br>_Component-specific finalization._  |
| virtual void | [**init\_impl**](#function-init_impl) () override<br>_Component-specific initialization (load model, read ICs)._  |
| virtual void | [**run\_impl**](#function-run_impl) (int dt) override<br>_Execute one time step of inference._  |
| virtual void | [**run\_inference**](#function-run_inference) (const std::vector&lt; double &gt; & inputs, std::vector&lt; double &gt; & outputs) override<br>_Run AI inference via the configured backend._  |


## Protected Functions inherited from emulator::EmulatorComp

See [emulator::EmulatorComp](classemulator_1_1EmulatorComp.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**export\_to\_coupler**](classemulator_1_1EmulatorComp.md#function-export_to_coupler) () <br>_Export fields to coupler (override as needed)._  |
| virtual void | [**final\_impl**](classemulator_1_1EmulatorComp.md#function-final_impl) () = 0<br>_Component-specific finalization._  |
| virtual void | [**import\_from\_coupler**](classemulator_1_1EmulatorComp.md#function-import_from_coupler) () <br>_Import fields from coupler (override as needed)._  |
| virtual void | [**init\_impl**](classemulator_1_1EmulatorComp.md#function-init_impl) () = 0<br>_Component-specific initialization._  |
| virtual void | [**run\_impl**](classemulator_1_1EmulatorComp.md#function-run_impl) (int dt) = 0<br>_Component-specific time step execution._  |
| virtual void | [**run\_inference**](classemulator_1_1EmulatorComp.md#function-run_inference) (const std::vector&lt; double &gt; & inputs, std::vector&lt; double &gt; & outputs) = 0<br>_Run AI inference on packed input/output vectors._  |






## Detailed Description


## Public Functions Documentation




### function EmulatorAtm 

```C++
emulator::EmulatorAtm::EmulatorAtm () 
```




<hr>



### function get\_inference\_config 

_Get current inference configuration._ 
```C++
inline const inference::InferenceConfig & emulator::EmulatorAtm::get_inference_config () const
```





**Returns:**

Reference to the current InferenceConfig 





        

<hr>



### function init\_coupling\_indices 

_Initialize coupling field indices from MCT field lists._ 
```C++
void emulator::EmulatorAtm::init_coupling_indices (
    const std::string & export_fields,
    const std::string & import_fields
) 
```



Parses the colon-separated field lists provided by the MCT layer and sets up internal index mappings for import/export operations.




**Parameters:**


* `export_fields` Colon-separated list of a2x field names 
* `import_fields` Colon-separated list of x2a field names 




        

<hr>



### function set\_inference\_config 

_Configure inference backend._ 
```C++
void emulator::EmulatorAtm::set_inference_config (
    const inference::InferenceConfig & config
) 
```



Must be called before [**initialize()**](classemulator_1_1EmulatorComp.md#function-initialize) to select the inference backend and set model parameters.




**Parameters:**


* `config` Inference configuration (backend type, model path, etc.) 




        

<hr>



### function set\_log\_file 

_Set logging file._ 
```C++
void emulator::EmulatorAtm::set_log_file (
    const std::string & filename
) 
```



Redirects the component logger to the specified file. If empty, logging goes to stdout.




**Parameters:**


* `filename` Path to log file 




        

<hr>



### function ~EmulatorAtm 

```C++
emulator::EmulatorAtm::~EmulatorAtm () override
```




<hr>
## Protected Functions Documentation




### function final\_impl 

_Component-specific finalization._ 
```C++
virtual void emulator::EmulatorAtm::final_impl () override
```



Implements [*emulator::EmulatorComp::final\_impl*](classemulator_1_1EmulatorComp.md#function-final_impl)


<hr>



### function init\_impl 

_Component-specific initialization (load model, read ICs)._ 
```C++
virtual void emulator::EmulatorAtm::init_impl () override
```



Implements [*emulator::EmulatorComp::init\_impl*](classemulator_1_1EmulatorComp.md#function-init_impl)


<hr>



### function run\_impl 

_Execute one time step of inference._ 
```C++
virtual void emulator::EmulatorAtm::run_impl (
    int dt
) override
```



Implements [*emulator::EmulatorComp::run\_impl*](classemulator_1_1EmulatorComp.md#function-run_impl)


<hr>



### function run\_inference 

_Run AI inference via the configured backend._ 
```C++
virtual void emulator::EmulatorAtm::run_inference (
    const std::vector< double > & inputs,
    std::vector< double > & outputs
) override
```



Run inference using the configured backend.


For spatial\_mode, the backend is called with batch\_size=1 and the full [C\*H\*W] flattened tensor. For pointwise mode, batch\_size=H\*W.




**Note:**

If inference fails, this function will abort the simulation with a clear error message. 





        
Implements [*emulator::EmulatorComp::run\_inference*](classemulator_1_1EmulatorComp.md#function-run_inference)


<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/eatm/src/emulator_atm.hpp`

