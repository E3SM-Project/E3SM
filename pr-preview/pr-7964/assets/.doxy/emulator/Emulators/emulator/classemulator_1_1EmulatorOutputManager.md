

# Class emulator::EmulatorOutputManager



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**EmulatorOutputManager**](classemulator_1_1EmulatorOutputManager.md)



_Manages all diagnostic output for an emulator component._ [More...](#detailed-description)

* `#include <emulator_output_manager.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**EmulatorOutputManager**](#function-emulatoroutputmanager) () = default<br> |
|  void | [**finalize**](#function-finalize) () <br>_Finalize and close all files._  |
|  std::string | [**find\_restart\_file**](#function-find_restart_file) (const std::string & rpointer\_dir, [**FileType**](namespaceemulator.md#enum-filetype) file\_type) const<br>_Find restart file from rpointer._  |
|  void | [**init\_timestep**](#function-init_timestep) (int current\_step, double dt) <br>_Called at the start of each timestep._  |
|  void | [**initialize**](#function-initialize) (const [**DiagnosticConfig**](structemulator_1_1DiagnosticConfig.md) & config, [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) comm, const std::vector&lt; int &gt; & col\_gids, int nlat, int nlon, const std::string & case\_name, const std::string & run\_dir, [**Logger**](classemulator_1_1Logger.md) & logger) <br>_Initialize the output manager._  |
|  bool | [**is\_initialized**](#function-is_initialized) () const<br>_Check if output manager is initialized._  |
|  bool | [**is\_restart\_step**](#function-is_restart_step) (int step) const<br>_Check if restart should be written at this step._  |
|  size\_t | [**num\_history\_streams**](#function-num_history_streams) () const<br>_Get number of history streams._  |
|  bool | [**read\_history\_restart**](#function-read_history_restart) (const std::string & filename) <br>_Read history restart from file._  |
|  bool | [**read\_restart**](#function-read_restart) (const std::string & filename) <br>_Read model restart file._  |
|  void | [**run**](#function-run) (int current\_step, const [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) & fields) <br>_Process output at current timestep._  |
|  void | [**setup**](#function-setup) (const [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) & fields) <br>_Set up output streams with field information._  |
|  bool | [**write\_history\_restart**](#function-write_history_restart) (int step) <br>_Write history restart for all non-instant streams._  |
|  bool | [**write\_restart**](#function-write_restart) (const [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) & fields, int step) <br>_Write model restart file._  |
|   | [**~EmulatorOutputManager**](#function-emulatoroutputmanager) () = default<br> |




























## Detailed Description


## Public Functions Documentation




### function EmulatorOutputManager 

```C++
emulator::EmulatorOutputManager::EmulatorOutputManager () = default
```




<hr>



### function finalize 

_Finalize and close all files._ 
```C++
void emulator::EmulatorOutputManager::finalize () 
```




<hr>



### function find\_restart\_file 

_Find restart file from rpointer._ 
```C++
std::string emulator::EmulatorOutputManager::find_restart_file (
    const std::string & rpointer_dir,
    FileType file_type
) const
```





**Parameters:**


* `rpointer_dir` Directory containing rpointer.atm 
* `file_type` RESTART or HISTORY\_RESTART 



**Returns:**

Filename if found, empty string otherwise 





        

<hr>



### function init\_timestep 

_Called at the start of each timestep._ 
```C++
void emulator::EmulatorOutputManager::init_timestep (
    int current_step,
    double dt
) 
```





**Parameters:**


* `current_step` Current model step 
* `dt` Timestep in seconds 




        

<hr>



### function initialize 

_Initialize the output manager._ 
```C++
void emulator::EmulatorOutputManager::initialize (
    const DiagnosticConfig & config,
    MPI_Comm comm,
    const std::vector< int > & col_gids,
    int nlat,
    int nlon,
    const std::string & case_name,
    const std::string & run_dir,
    Logger & logger
) 
```





**Parameters:**


* `config` Diagnostic configuration 
* `comm` MPI communicator 
* `col_gids` Global IDs of local columns 
* `nlat` Number of latitudes 
* `nlon` Number of longitudes 
* `case_name` Case name for filenames 
* `run_dir` Run directory for output files 
* `logger` [**Logger**](classemulator_1_1Logger.md) reference 




        

<hr>



### function is\_initialized 

_Check if output manager is initialized._ 
```C++
inline bool emulator::EmulatorOutputManager::is_initialized () const
```




<hr>



### function is\_restart\_step 

_Check if restart should be written at this step._ 
```C++
bool emulator::EmulatorOutputManager::is_restart_step (
    int step
) const
```




<hr>



### function num\_history\_streams 

_Get number of history streams._ 
```C++
inline size_t emulator::EmulatorOutputManager::num_history_streams () const
```




<hr>



### function read\_history\_restart 

_Read history restart from file._ 
```C++
bool emulator::EmulatorOutputManager::read_history_restart (
    const std::string & filename
) 
```





**Parameters:**


* `filename` History restart file path 



**Returns:**

true if successful 





        

<hr>



### function read\_restart 

_Read model restart file._ 
```C++
bool emulator::EmulatorOutputManager::read_restart (
    const std::string & filename
) 
```





**Parameters:**


* `filename` Restart file path 
* `fields` Field data receiver (mutable) 



**Returns:**

true if successful 





        

<hr>



### function run 

_Process output at current timestep._ 
```C++
void emulator::EmulatorOutputManager::run (
    int current_step,
    const FieldDataProvider & fields
) 
```



Updates averaging buffers and writes output if needed.




**Parameters:**


* `current_step` Current model step 
* `fields` Field data provider 




        

<hr>



### function setup 

_Set up output streams with field information._ 
```C++
void emulator::EmulatorOutputManager::setup (
    const FieldDataProvider & fields
) 
```





**Parameters:**


* `fields` Field data provider 




        

<hr>



### function write\_history\_restart 

_Write history restart for all non-instant streams._ 
```C++
bool emulator::EmulatorOutputManager::write_history_restart (
    int step
) 
```





**Parameters:**


* `step` Current step (for filename) 



**Returns:**

true if successful 





        

<hr>



### function write\_restart 

_Write model restart file._ 
```C++
bool emulator::EmulatorOutputManager::write_restart (
    const FieldDataProvider & fields,
    int step
) 
```



Writes all prognostic fields to a restart file and updates rpointer.atm.




**Parameters:**


* `fields` Field data provider 
* `step` Current step (for filename) 



**Returns:**

true if successful 





        

<hr>



### function ~EmulatorOutputManager 

```C++
emulator::EmulatorOutputManager::~EmulatorOutputManager () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_output_manager.hpp`

