

# Class emulator::EmulatorComp



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**EmulatorComp**](classemulator_1_1EmulatorComp.md)



_Abstract base class for all emulated E3SM components._ [More...](#detailed-description)

* `#include <emulator_comp.hpp>`





Inherited by the following classes: [emulator::EmulatorAtm](classemulator_1_1EmulatorAtm.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**EmulatorComp**](#function-emulatorcomp) ([**CompType**](namespaceemulator.md#enum-comptype) type) <br>_Construct an emulator component of the given type._  |
|  [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) | [**comm**](#function-comm) () const<br>_Get MPI communicator._  |
|  int | [**comp\_id**](#function-comp_id) () const<br>_Get component ID._  |
|  void | [**create\_instance**](#function-create_instance) ([**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) comm, int comp\_id, const char \* input\_file, int run\_type, int start\_ymd, int start\_tod) <br>_Initialize the component instance._  |
|  void | [**finalize**](#function-finalize) () <br>_Finalize and clean up the component._  |
|  void | [**get\_cols\_area**](#function-get_cols_area) (double \* area) const<br>_Get area for local columns._  |
|  void | [**get\_cols\_latlon**](#function-get_cols_latlon) (double \* lat, double \* lon) const<br>_Get lat/lon for local columns._  |
|  void | [**get\_local\_col\_gids**](#function-get_local_col_gids) (int \* gids) const<br>_Get global IDs for local columns._  |
|  int | [**get\_num\_global\_cols**](#function-get_num_global_cols) () const<br>_Get total number of columns globally._  |
|  int | [**get\_num\_local\_cols**](#function-get_num_local_cols) () const<br>_Get number of columns on this MPI rank._  |
|  int | [**get\_nx**](#function-get_nx) () const<br>_Get number of grid points in x-direction (longitude)._  |
|  int | [**get\_ny**](#function-get_ny) () const<br>_Get number of grid points in y-direction (latitude)._  |
|  void | [**initialize**](#function-initialize) () <br>_Initialize the component (phase 2)._  |
|  bool | [**is\_root**](#function-is_root) () const<br>_Check if this is the root rank (rank 0)._  |
|  int | [**rank**](#function-rank) () const<br>_Get MPI rank within component communicator._  |
|  void | [**run**](#function-run) (int dt) <br>_Execute one time step._  |
|  void | [**set\_grid\_data**](#function-set_grid_data) (int nx, int ny, int num\_local\_cols, int num\_global\_cols, const int \* col\_gids, const double \* lat, const double \* lon, const double \* area) <br>_Set grid decomposition data from the driver._  |
|  void | [**setup\_coupling**](#function-setup_coupling) (double \* import\_data, double \* export\_data, int num\_imports, int num\_exports, int field\_size) <br>_Set up coupling data buffers._  |
|  [**CompType**](namespaceemulator.md#enum-comptype) | [**type**](#function-type) () const<br>_Get component type._  |
| virtual  | [**~EmulatorComp**](#function-emulatorcomp) () = default<br> |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; double &gt; | [**m\_area**](#variable-m_area)  <br>_Cell area [km²]._  |
|  std::vector&lt; int &gt; | [**m\_col\_gids**](#variable-m_col_gids)  <br>_Global IDs for local columns._  |
|  [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) | [**m\_comm**](#variable-m_comm)  <br>_MPI communicator._  |
|  int | [**m\_comp\_id**](#variable-m_comp_id)  <br>_Component ID from driver._  |
|  int | [**m\_current\_tod**](#variable-m_current_tod)   = `0`<br>_Current time of day [seconds]._  |
|  int | [**m\_current\_ymd**](#variable-m_current_ymd)   = `0`<br>_Current date as YYYYMMDD._  |
|  double \* | [**m\_export\_data**](#variable-m_export_data)   = `nullptr`<br>_Export buffer pointer (a2x)_  |
|  int | [**m\_field\_size**](#variable-m_field_size)   = `0`<br>_Size per field (should = num\_local\_cols)_  |
|  double \* | [**m\_import\_data**](#variable-m_import_data)   = `nullptr`<br>_Import buffer pointer (x2a)_  |
|  std::string | [**m\_input\_file**](#variable-m_input_file)  <br>_Path to configuration file._  |
|  std::vector&lt; double &gt; | [**m\_lat**](#variable-m_lat)  <br>_Latitude [radians]._  |
|  [**Logger**](classemulator_1_1Logger.md) | [**m\_logger**](#variable-m_logger)  <br>_Component logger._  |
|  std::vector&lt; double &gt; | [**m\_lon**](#variable-m_lon)  <br>_Longitude [radians]._  |
|  int | [**m\_nprocs**](#variable-m_nprocs)  <br>_Number of MPI processes._  |
|  int | [**m\_num\_exports**](#variable-m_num_exports)   = `0`<br>_Number of export fields._  |
|  int | [**m\_num\_global\_cols**](#variable-m_num_global_cols)   = `0`<br>_Global column count._  |
|  int | [**m\_num\_imports**](#variable-m_num_imports)   = `0`<br>_Number of import fields._  |
|  int | [**m\_num\_local\_cols**](#variable-m_num_local_cols)   = `0`<br>_Local column count._  |
|  int | [**m\_nx**](#variable-m_nx)   = `0`<br>_Grid points in x (longitude)_  |
|  int | [**m\_ny**](#variable-m_ny)   = `0`<br>_Grid points in y (latitude)_  |
|  int | [**m\_rank**](#variable-m_rank)  <br>_MPI rank within component._  |
|  int | [**m\_run\_type**](#variable-m_run_type)  <br>_Run type (startup/continue/branch)_  |
|  int | [**m\_step\_count**](#variable-m_step_count)   = `0`<br>_Number of steps executed._  |
|  [**CompType**](namespaceemulator.md#enum-comptype) | [**m\_type**](#variable-m_type)  <br>_Component type enum._  |
















## Protected Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**export\_to\_coupler**](#function-export_to_coupler) () <br>_Export fields to coupler (override as needed)._  |
| virtual void | [**final\_impl**](#function-final_impl) () = 0<br>_Component-specific finalization._  |
| virtual void | [**import\_from\_coupler**](#function-import_from_coupler) () <br>_Import fields from coupler (override as needed)._  |
| virtual void | [**init\_impl**](#function-init_impl) () = 0<br>_Component-specific initialization._  |
| virtual void | [**run\_impl**](#function-run_impl) (int dt) = 0<br>_Component-specific time step execution._  |
| virtual void | [**run\_inference**](#function-run_inference) (const std::vector&lt; double &gt; & inputs, std::vector&lt; double &gt; & outputs) = 0<br>_Run AI inference on packed input/output vectors._  |




## Detailed Description


## Public Functions Documentation




### function EmulatorComp 

_Construct an emulator component of the given type._ 
```C++
explicit emulator::EmulatorComp::EmulatorComp (
    CompType type
) 
```





**Parameters:**


* `type` Component type (ATM, OCN, ICE, LND) 




        

<hr>



### function comm 

_Get MPI communicator._ 
```C++
inline MPI_Comm emulator::EmulatorComp::comm () const
```




<hr>



### function comp\_id 

_Get component ID._ 
```C++
inline int emulator::EmulatorComp::comp_id () const
```




<hr>



### function create\_instance 

_Initialize the component instance._ 
```C++
void emulator::EmulatorComp::create_instance (
    MPI_Comm comm,
    int comp_id,
    const char * input_file,
    int run_type,
    int start_ymd,
    int start_tod
) 
```



Initialize the component instance with MPI and metadata.


Sets up MPI communicator, stores component metadata, and prepares for grid data and coupling setup.




**Parameters:**


* `comm` MPI communicator for this component 
* `comp_id` Component ID assigned by the driver 
* `input_file` Path to component configuration file 
* `run_type` Run type (startup, continue, branch) 
* `start_ymd` Start date as YYYYMMDD 
* `start_tod` Start time of day in seconds 




        

<hr>



### function finalize 

_Finalize and clean up the component._ 
```C++
void emulator::EmulatorComp::finalize () 
```



Finalize the component and release resources.


Releases resources, finalizes inference backend, and performs any necessary cleanup. 


        

<hr>



### function get\_cols\_area 

_Get area for local columns._ 
```C++
void emulator::EmulatorComp::get_cols_area (
    double * area
) const
```





**Parameters:**


* `area` Output area array [km²] 




        

<hr>



### function get\_cols\_latlon 

_Get lat/lon for local columns._ 
```C++
void emulator::EmulatorComp::get_cols_latlon (
    double * lat,
    double * lon
) const
```





**Parameters:**


* `lat` Output latitude array [radians] 
* `lon` Output longitude array [radians] 




        

<hr>



### function get\_local\_col\_gids 

_Get global IDs for local columns._ 
```C++
void emulator::EmulatorComp::get_local_col_gids (
    int * gids
) const
```





**Parameters:**


* `gids` Output array (must be pre-allocated with num\_local\_cols) 




        

<hr>



### function get\_num\_global\_cols 

_Get total number of columns globally._ 
```C++
inline int emulator::EmulatorComp::get_num_global_cols () const
```




<hr>



### function get\_num\_local\_cols 

_Get number of columns on this MPI rank._ 
```C++
inline int emulator::EmulatorComp::get_num_local_cols () const
```




<hr>



### function get\_nx 

_Get number of grid points in x-direction (longitude)._ 
```C++
inline int emulator::EmulatorComp::get_nx () const
```




<hr>



### function get\_ny 

_Get number of grid points in y-direction (latitude)._ 
```C++
inline int emulator::EmulatorComp::get_ny () const
```




<hr>



### function initialize 

_Initialize the component (phase 2)._ 
```C++
void emulator::EmulatorComp::initialize () 
```



Initialize the component by calling the derived class implementation.


Called after grid and coupling setup. Loads AI model, reads initial conditions, and prepares for time stepping. 


        

<hr>



### function is\_root 

_Check if this is the root rank (rank 0)._ 
```C++
inline bool emulator::EmulatorComp::is_root () const
```




<hr>



### function rank 

_Get MPI rank within component communicator._ 
```C++
inline int emulator::EmulatorComp::rank () const
```




<hr>



### function run 

_Execute one time step._ 
```C++
void emulator::EmulatorComp::run (
    int dt
) 
```



Execute one time step: import → run → export.


Performs the main computational cycle:



* Import fields from coupler
* Run AI inference
* Export fields to coupler






**Parameters:**


* `dt` Time step size in seconds 




        

<hr>



### function set\_grid\_data 

_Set grid decomposition data from the driver._ 
```C++
void emulator::EmulatorComp::set_grid_data (
    int nx,
    int ny,
    int num_local_cols,
    int num_global_cols,
    const int * col_gids,
    const double * lat,
    const double * lon,
    const double * area
) 
```



Override grid data with values from Fortran driver.


Receives the local portion of the global grid as distributed by the E3SM driver's domain decomposition.




**Parameters:**


* `nx` Number of grid points in x-direction (longitude) 
* `ny` Number of grid points in y-direction (latitude) 
* `num_local_cols` Number of columns on this MPI rank 
* `num_global_cols` Total number of columns globally 
* `col_gids` Global IDs for each local column 
* `lat` Latitude values for each local column [radians] 
* `lon` Longitude values for each local column [radians] 
* `area` Grid cell areas for each local column [km²]

Used when the driver provides grid decomposition directly instead of reading from a file. 


        

<hr>



### function setup\_coupling 

_Set up coupling data buffers._ 
```C++
void emulator::EmulatorComp::setup_coupling (
    double * import_data,
    double * export_data,
    int num_imports,
    int num_exports,
    int field_size
) 
```



Set up coupling buffer pointers from the MCT layer.


Receives pointers to the MCT attribute vector buffers for data exchange with other components via the coupler.




**Parameters:**


* `import_data` Pointer to import buffer (x2a fields) 
* `export_data` Pointer to export buffer (a2x fields) 
* `num_imports` Number of import fields 
* `num_exports` Number of export fields 
* `field_size` Size of each field (should equal num\_local\_cols) 




        

<hr>



### function type 

_Get component type._ 
```C++
inline CompType emulator::EmulatorComp::type () const
```




<hr>



### function ~EmulatorComp 

```C++
virtual emulator::EmulatorComp::~EmulatorComp () = default
```




<hr>
## Protected Attributes Documentation




### variable m\_area 

_Cell area [km²]._ 
```C++
std::vector<double> emulator::EmulatorComp::m_area;
```




<hr>



### variable m\_col\_gids 

_Global IDs for local columns._ 
```C++
std::vector<int> emulator::EmulatorComp::m_col_gids;
```




<hr>



### variable m\_comm 

_MPI communicator._ 
```C++
MPI_Comm emulator::EmulatorComp::m_comm;
```




<hr>



### variable m\_comp\_id 

_Component ID from driver._ 
```C++
int emulator::EmulatorComp::m_comp_id;
```




<hr>



### variable m\_current\_tod 

_Current time of day [seconds]._ 
```C++
int emulator::EmulatorComp::m_current_tod;
```




<hr>



### variable m\_current\_ymd 

_Current date as YYYYMMDD._ 
```C++
int emulator::EmulatorComp::m_current_ymd;
```




<hr>



### variable m\_export\_data 

_Export buffer pointer (a2x)_ 
```C++
double* emulator::EmulatorComp::m_export_data;
```




<hr>



### variable m\_field\_size 

_Size per field (should = num\_local\_cols)_ 
```C++
int emulator::EmulatorComp::m_field_size;
```




<hr>



### variable m\_import\_data 

_Import buffer pointer (x2a)_ 
```C++
double* emulator::EmulatorComp::m_import_data;
```




<hr>



### variable m\_input\_file 

_Path to configuration file._ 
```C++
std::string emulator::EmulatorComp::m_input_file;
```




<hr>



### variable m\_lat 

_Latitude [radians]._ 
```C++
std::vector<double> emulator::EmulatorComp::m_lat;
```




<hr>



### variable m\_logger 

_Component logger._ 
```C++
Logger emulator::EmulatorComp::m_logger;
```




<hr>



### variable m\_lon 

_Longitude [radians]._ 
```C++
std::vector<double> emulator::EmulatorComp::m_lon;
```




<hr>



### variable m\_nprocs 

_Number of MPI processes._ 
```C++
int emulator::EmulatorComp::m_nprocs;
```




<hr>



### variable m\_num\_exports 

_Number of export fields._ 
```C++
int emulator::EmulatorComp::m_num_exports;
```




<hr>



### variable m\_num\_global\_cols 

_Global column count._ 
```C++
int emulator::EmulatorComp::m_num_global_cols;
```




<hr>



### variable m\_num\_imports 

_Number of import fields._ 
```C++
int emulator::EmulatorComp::m_num_imports;
```




<hr>



### variable m\_num\_local\_cols 

_Local column count._ 
```C++
int emulator::EmulatorComp::m_num_local_cols;
```




<hr>



### variable m\_nx 

_Grid points in x (longitude)_ 
```C++
int emulator::EmulatorComp::m_nx;
```




<hr>



### variable m\_ny 

_Grid points in y (latitude)_ 
```C++
int emulator::EmulatorComp::m_ny;
```




<hr>



### variable m\_rank 

_MPI rank within component._ 
```C++
int emulator::EmulatorComp::m_rank;
```




<hr>



### variable m\_run\_type 

_Run type (startup/continue/branch)_ 
```C++
int emulator::EmulatorComp::m_run_type;
```




<hr>



### variable m\_step\_count 

_Number of steps executed._ 
```C++
int emulator::EmulatorComp::m_step_count;
```




<hr>



### variable m\_type 

_Component type enum._ 
```C++
CompType emulator::EmulatorComp::m_type;
```




<hr>
## Protected Functions Documentation




### function export\_to\_coupler 

_Export fields to coupler (override as needed)._ 
```C++
inline virtual void emulator::EmulatorComp::export_to_coupler () 
```




<hr>



### function final\_impl 

_Component-specific finalization._ 
```C++
virtual void emulator::EmulatorComp::final_impl () = 0
```




<hr>



### function import\_from\_coupler 

_Import fields from coupler (override as needed)._ 
```C++
inline virtual void emulator::EmulatorComp::import_from_coupler () 
```




<hr>



### function init\_impl 

_Component-specific initialization._ 
```C++
virtual void emulator::EmulatorComp::init_impl () = 0
```




<hr>



### function run\_impl 

_Component-specific time step execution._ 
```C++
virtual void emulator::EmulatorComp::run_impl (
    int dt
) = 0
```




<hr>



### function run\_inference 

_Run AI inference on packed input/output vectors._ 
```C++
virtual void emulator::EmulatorComp::run_inference (
    const std::vector< double > & inputs,
    std::vector< double > & outputs
) = 0
```





**Parameters:**


* `inputs` Packed input features 
* `outputs` Packed output features (will be resized as needed) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_comp.hpp`

