

# File emulator\_atm\_interface.cpp



[**FileList**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**emulator\_atm\_interface.cpp**](emulator__atm__interface_8cpp.md)

[Go to the source code of this file](emulator__atm__interface_8cpp_source.md)

_C interface for the atmosphere emulator (Fortran-callable)._ [More...](#detailed-description)

* `#include "../../common/src/emulator_comp.hpp"`
* `#include "../../common/src/emulator_context.hpp"`
* `#include "emulator_atm.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**emulator\_atm\_create\_instance**](#function-emulator_atm_create_instance) (const MPI\_Fint f\_comm, const int comp\_id, const char \* input\_file, const char \* log\_file, const int run\_type, const int start\_ymd, const int start\_tod) <br>_Create and initialize the atmosphere emulator instance._  |
|  void | [**emulator\_atm\_finalize**](#function-emulator_atm_finalize) () <br>_Finalize and clean up the atmosphere emulator._  |
|  void | [**emulator\_atm\_get\_cols\_area**](#function-emulator_atm_get_cols_area) (double \* area) <br>_Get cell areas for local columns._  |
|  void | [**emulator\_atm\_get\_cols\_latlon**](#function-emulator_atm_get_cols_latlon) (double \* lat, double \* lon) <br>_Get latitude/longitude for local columns._  |
|  void | [**emulator\_atm\_get\_local\_cols\_gids**](#function-emulator_atm_get_local_cols_gids) (int \* gids) <br>_Get global IDs for local columns._  |
|  int | [**emulator\_atm\_get\_num\_global\_cols**](#function-emulator_atm_get_num_global_cols) () <br>_Get total number of global columns._  |
|  int | [**emulator\_atm\_get\_num\_local\_cols**](#function-emulator_atm_get_num_local_cols) () <br>_Get number of local columns on this rank._  |
|  int | [**emulator\_atm\_get\_nx**](#function-emulator_atm_get_nx) () <br>_Get grid size in x-direction._  |
|  int | [**emulator\_atm\_get\_ny**](#function-emulator_atm_get_ny) () <br>_Get grid size in y-direction._  |
|  void | [**emulator\_atm\_init**](#function-emulator_atm_init) () <br>_Initialize the atmosphere emulator (phase 2)._  |
|  void | [**emulator\_atm\_init\_coupling\_indices**](#function-emulator_atm_init_coupling_indices) (const char \* export\_fields, const char \* import\_fields) <br>_Initialize coupling field index mappings._  |
|  void | [**emulator\_atm\_run**](#function-emulator_atm_run) (const int dt) <br>_Execute one time step._  |
|  void | [**emulator\_atm\_set\_grid\_data**](#function-emulator_atm_set_grid_data) (const int nx, const int ny, const int num\_local\_cols, const int num\_global\_cols, const int \* col\_gids, const double \* lat, const double \* lon, const double \* area) <br>_Set grid decomposition data from driver._  |
|  void | [**emulator\_atm\_setup\_coupling**](#function-emulator_atm_setup_coupling) (double \* import\_data, double \* export\_data, const int num\_imports, const int num\_exports, const int field\_size) <br>_Set up coupling buffer pointers from MCT._  |




























## Detailed Description


Provides extern "C" functions that can be called from Fortran via ISO\_C\_BINDING. These functions wrap the EmulatorAtm C++ class methods.


The interface follows the E3SM component lifecycle:
* create\_instance - Initialize MPI and component
* set\_grid\_data - Set grid decomposition
* init\_coupling\_indices - Parse MCT field lists
* setup\_coupling - Set coupling buffer pointers
* init - Load model and prepare for time stepping
* run - Execute time steps
* finalize - Clean up resources






**See also:** EmulatorAtm for the C++ implementation 


**See also:** emulator\_atm\_f2c.F90 for the Fortran interface module 



    
## Public Functions Documentation




### function emulator\_atm\_create\_instance 

_Create and initialize the atmosphere emulator instance._ 
```C++
void emulator_atm_create_instance (
    const MPI_Fint f_comm,
    const int comp_id,
    const char * input_file,
    const char * log_file,
    const int run_type,
    const int start_ymd,
    const int start_tod
) 
```





**Parameters:**


* `f_comm` Fortran MPI communicator 
* `comp_id` Component ID from driver 
* `input_file` Path to atm\_in configuration file 
* `log_file` Path to log file (NULL for stdout) 
* `run_type` Run type (startup, continue, branch) 
* `start_ymd` Start date as YYYYMMDD 
* `start_tod` Start time of day in seconds 




        

<hr>



### function emulator\_atm\_finalize 

_Finalize and clean up the atmosphere emulator._ 
```C++
void emulator_atm_finalize () 
```




<hr>



### function emulator\_atm\_get\_cols\_area 

_Get cell areas for local columns._ 
```C++
void emulator_atm_get_cols_area (
    double * area
) 
```





**Parameters:**


* `area` Output area array 




        

<hr>



### function emulator\_atm\_get\_cols\_latlon 

_Get latitude/longitude for local columns._ 
```C++
void emulator_atm_get_cols_latlon (
    double * lat,
    double * lon
) 
```





**Parameters:**


* `lat` Output latitude array [radians] 
* `lon` Output longitude array [radians] 




        

<hr>



### function emulator\_atm\_get\_local\_cols\_gids 

_Get global IDs for local columns._ 
```C++
void emulator_atm_get_local_cols_gids (
    int * gids
) 
```





**Parameters:**


* `gids` Output array (must be pre-allocated) 




        

<hr>



### function emulator\_atm\_get\_num\_global\_cols 

_Get total number of global columns._ 
```C++
int emulator_atm_get_num_global_cols () 
```




<hr>



### function emulator\_atm\_get\_num\_local\_cols 

_Get number of local columns on this rank._ 
```C++
int emulator_atm_get_num_local_cols () 
```




<hr>



### function emulator\_atm\_get\_nx 

_Get grid size in x-direction._ 
```C++
int emulator_atm_get_nx () 
```




<hr>



### function emulator\_atm\_get\_ny 

_Get grid size in y-direction._ 
```C++
int emulator_atm_get_ny () 
```




<hr>



### function emulator\_atm\_init 

_Initialize the atmosphere emulator (phase 2)._ 
```C++
void emulator_atm_init () 
```




<hr>



### function emulator\_atm\_init\_coupling\_indices 

_Initialize coupling field index mappings._ 
```C++
void emulator_atm_init_coupling_indices (
    const char * export_fields,
    const char * import_fields
) 
```




<hr>



### function emulator\_atm\_run 

_Execute one time step._ 
```C++
void emulator_atm_run (
    const int dt
) 
```





**Parameters:**


* `dt` Time step size in seconds 




        

<hr>



### function emulator\_atm\_set\_grid\_data 

_Set grid decomposition data from driver._ 
```C++
void emulator_atm_set_grid_data (
    const int nx,
    const int ny,
    const int num_local_cols,
    const int num_global_cols,
    const int * col_gids,
    const double * lat,
    const double * lon,
    const double * area
) 
```




<hr>



### function emulator\_atm\_setup\_coupling 

_Set up coupling buffer pointers from MCT._ 
```C++
void emulator_atm_setup_coupling (
    double * import_data,
    double * export_data,
    const int num_imports,
    const int num_exports,
    const int field_size
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/eatm/src/emulator_atm_interface.cpp`

