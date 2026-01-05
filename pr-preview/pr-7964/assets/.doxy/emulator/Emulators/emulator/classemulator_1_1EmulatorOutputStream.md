

# Class emulator::EmulatorOutputStream



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**EmulatorOutputStream**](classemulator_1_1EmulatorOutputStream.md)



_Manages a single diagnostic output stream._ [More...](#detailed-description)

* `#include <emulator_output_stream.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**EmulatorOutputStream**](#function-emulatoroutputstream) () = default<br> |
|  const [**OutputStreamConfig**](structemulator_1_1OutputStreamConfig.md) & | [**config**](#function-config) () const<br>_Get the stream configuration._  |
|  void | [**finalize**](#function-finalize) () <br>_Finalize and close any open files._  |
|  void | [**init\_timestep**](#function-init_timestep) (int current\_step, double dt) <br>_Called at the start of each timestep._  |
|  void | [**initialize**](#function-initialize) (const [**OutputStreamConfig**](structemulator_1_1OutputStreamConfig.md) & config, [**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) comm, const std::vector&lt; int &gt; & col\_gids, int nlat, int nlon, [**Logger**](classemulator_1_1Logger.md) & logger) <br>_Initialize the output stream._  |
|  bool | [**is\_write\_step**](#function-is_write_step) () const<br>_Check if this step is a write step._  |
|  bool | [**needs\_history\_restart**](#function-needs_history_restart) () const<br>_Check if stream needs history restart (non-instant averaging)._  |
|  bool | [**read\_history\_restart**](#function-read_history_restart) (const std::string & filename) <br>_Read averaging state from history restart file._  |
|  void | [**run**](#function-run) (int current\_step, const [**FieldDataProvider**](classemulator_1_1FieldDataProvider.md) & fields, const std::string & case\_name) <br>_Process fields at current timestep._  |
|  bool | [**write\_history\_restart**](#function-write_history_restart) (const std::string & filename) <br>_Write averaging state to history restart file._  |
|   | [**~EmulatorOutputStream**](#function-emulatoroutputstream) () = default<br> |




























## Detailed Description


Handles:



* NetCDF file creation and management
* Accumulation buffer for averaging (AVERAGE, MIN, MAX, STD, SUM)
* Write timing based on configured frequency
* Field stacking for sliced variables 




    
## Public Functions Documentation




### function EmulatorOutputStream 

```C++
emulator::EmulatorOutputStream::EmulatorOutputStream () = default
```




<hr>



### function config 

_Get the stream configuration._ 
```C++
inline const OutputStreamConfig & emulator::EmulatorOutputStream::config () const
```




<hr>



### function finalize 

_Finalize and close any open files._ 
```C++
void emulator::EmulatorOutputStream::finalize () 
```




<hr>



### function init\_timestep 

_Called at the start of each timestep._ 
```C++
void emulator::EmulatorOutputStream::init_timestep (
    int current_step,
    double dt
) 
```





**Parameters:**


* `current_step` Current model step number 
* `dt` Timestep in seconds 




        

<hr>



### function initialize 

_Initialize the output stream._ 
```C++
void emulator::EmulatorOutputStream::initialize (
    const OutputStreamConfig & config,
    MPI_Comm comm,
    const std::vector< int > & col_gids,
    int nlat,
    int nlon,
    Logger & logger
) 
```





**Parameters:**


* `config` Stream configuration 
* `comm` MPI communicator 
* `col_gids` Global IDs of local columns 
* `nlat` Number of latitudes (for global grid) 
* `nlon` Number of longitudes 
* `logger` [**Logger**](classemulator_1_1Logger.md) reference 




        

<hr>



### function is\_write\_step 

_Check if this step is a write step._ 
```C++
inline bool emulator::EmulatorOutputStream::is_write_step () const
```




<hr>



### function needs\_history\_restart 

_Check if stream needs history restart (non-instant averaging)._ 
```C++
inline bool emulator::EmulatorOutputStream::needs_history_restart () const
```




<hr>



### function read\_history\_restart 

_Read averaging state from history restart file._ 
```C++
bool emulator::EmulatorOutputStream::read_history_restart (
    const std::string & filename
) 
```





**Parameters:**


* `filename` Input filename 



**Returns:**

true if successful 





        

<hr>



### function run 

_Process fields at current timestep._ 
```C++
void emulator::EmulatorOutputStream::run (
    int current_step,
    const FieldDataProvider & fields,
    const std::string & case_name
) 
```



Updates averaging buffers if not an instant stream. Writes output if this is a write step.




**Parameters:**


* `current_step` Current step number 
* `fields` Field data provider 
* `case_name` Case name for filename generation 




        

<hr>



### function write\_history\_restart 

_Write averaging state to history restart file._ 
```C++
bool emulator::EmulatorOutputStream::write_history_restart (
    const std::string & filename
) 
```





**Parameters:**


* `filename` Output filename 



**Returns:**

true if successful 





        

<hr>



### function ~EmulatorOutputStream 

```C++
emulator::EmulatorOutputStream::~EmulatorOutputStream () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_output_stream.hpp`

