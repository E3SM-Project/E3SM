

# Class emulator::EmulatorIO



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**EmulatorIO**](classemulator_1_1EmulatorIO.md)



_Static class providing parallel I/O using SCORPIO/PIO._ [More...](#detailed-description)

* `#include <emulator_io.hpp>`







































## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**close\_file**](#function-close_file) (int ncid) <br>_Close an open file._  |
|  int | [**create\_file**](#function-create_file) (const std::string & filename) <br>_Create a new NetCDF file for writing._  |
|  int | [**define\_dim**](#function-define_dim) (int ncid, const std::string & dimname, int length) <br>_Define a new dimension in a NetCDF file._  |
|  int | [**define\_var**](#function-define_var) (int ncid, const std::string & varname, int nctype, const std::vector&lt; int &gt; & dimids) <br>_Define a new variable in a NetCDF file._  |
|  bool | [**end\_def**](#function-end_def) (int ncid) <br>_End define mode and switch to data mode._  |
|  void | [**finalize**](#function-finalize) () <br>_Finalize the PIO subsystem._  |
|  int | [**get\_dim\_size**](#function-get_dim_size) (int ncid, const std::string & dimname) <br>_Get the size of a dimension._  |
|  bool | [**has\_var**](#function-has_var) (int ncid, const std::string & varname) <br>_Check if a variable exists in the file._  |
|  void | [**initialize**](#function-initialize) ([**MPI\_Comm**](inference__backend_8hpp.md#typedef-mpi_comm) comm, const std::string & comp\_name) <br>_Initialize the PIO subsystem._  |
|  bool | [**is\_initialized**](#function-is_initialized) () <br>_Check if PIO is initialized._  |
|  int | [**open\_file**](#function-open_file) (const std::string & filename) <br>_Open an existing NetCDF file for reading._  |
|  bool | [**read\_var\_1d**](#function-read_var_1d) (int ncid, const std::string & varname, double \* data, int size) <br>_Read a 1D double array from a NetCDF variable._  |
|  bool | [**read\_var\_1d\_int**](#function-read_var_1d_int) (int ncid, const std::string & varname, int \* data, int size) <br>_Read a 1D integer array from a NetCDF variable._  |
|  bool | [**read\_var\_2d**](#function-read_var_2d) (int ncid, const std::string & varname, double \* data, int nx, int ny) <br>_Read a 2D double array from a NetCDF variable._  |
|  bool | [**read\_var\_3d\_slice**](#function-read_var_3d_slice) (int ncid, const std::string & varname, double \* data, int nx, int ny, int time\_idx) <br>_Read a 2D slice from a 3D variable (time, lat, lon)._  |
|  void | [**sync\_file**](#function-sync_file) (int ncid) <br>_Synchronize (flush) file contents to disk._  |
|  bool | [**write\_var\_1d**](#function-write_var_1d) (int ncid, const std::string & varname, const double \* data, int size) <br>_Write a 1D double array to a NetCDF variable._  |
|  bool | [**write\_var\_2d**](#function-write_var_2d) (int ncid, const std::string & varname, const double \* data, int nx, int ny) <br>_Write a 2D double array to a NetCDF variable._  |


























## Detailed Description


## Public Static Functions Documentation




### function close\_file 

_Close an open file._ 
```C++
static void emulator::EmulatorIO::close_file (
    int ncid
) 
```





**Parameters:**


* `ncid` File ID from [**open\_file()**](classemulator_1_1EmulatorIO.md#function-open_file) or [**create\_file()**](classemulator_1_1EmulatorIO.md#function-create_file) 




        

<hr>



### function create\_file 

_Create a new NetCDF file for writing._ 
```C++
static int emulator::EmulatorIO::create_file (
    const std::string & filename
) 
```





**Parameters:**


* `filename` Path to the new file 



**Returns:**

NetCDF file ID (ncid), or -1 on failure 





        

<hr>



### function define\_dim 

_Define a new dimension in a NetCDF file._ 
```C++
static int emulator::EmulatorIO::define_dim (
    int ncid,
    const std::string & dimname,
    int length
) 
```





**Parameters:**


* `ncid` File ID (must be in define mode) 
* `dimname` Dimension name 
* `length` Dimension length 



**Returns:**

Dimension ID, or -1 on failure 





        

<hr>



### function define\_var 

_Define a new variable in a NetCDF file._ 
```C++
static int emulator::EmulatorIO::define_var (
    int ncid,
    const std::string & varname,
    int nctype,
    const std::vector< int > & dimids
) 
```





**Parameters:**


* `ncid` File ID (must be in define mode) 
* `varname` Variable name 
* `nctype` NetCDF type (e.g., NC\_DOUBLE) 
* `dimids` Vector of dimension IDs 



**Returns:**

Variable ID, or -1 on failure 





        

<hr>



### function end\_def 

_End define mode and switch to data mode._ 
```C++
static bool emulator::EmulatorIO::end_def (
    int ncid
) 
```



Must be called after defining dimensions and variables, before writing data. 

**Parameters:**


* `ncid` File ID 



**Returns:**

true on success 





        

<hr>



### function finalize 

_Finalize the PIO subsystem._ 
```C++
static void emulator::EmulatorIO::finalize () 
```



Should be called before MPI\_Finalize to release resources. 


        

<hr>



### function get\_dim\_size 

_Get the size of a dimension._ 
```C++
static int emulator::EmulatorIO::get_dim_size (
    int ncid,
    const std::string & dimname
) 
```





**Parameters:**


* `ncid` File ID 
* `dimname` Dimension name 



**Returns:**

Dimension length, or -1 if not found 





        

<hr>



### function has\_var 

_Check if a variable exists in the file._ 
```C++
static bool emulator::EmulatorIO::has_var (
    int ncid,
    const std::string & varname
) 
```





**Parameters:**


* `ncid` File ID 
* `varname` Variable name 



**Returns:**

true if variable exists 





        

<hr>



### function initialize 

_Initialize the PIO subsystem._ 
```C++
static void emulator::EmulatorIO::initialize (
    MPI_Comm comm,
    const std::string & comp_name
) 
```



Must be called once before any file operations. Typically called from [**EmulatorComp::create\_instance()**](classemulator_1_1EmulatorComp.md#function-create_instance).




**Parameters:**


* `comm` MPI communicator for parallel I/O 
* `comp_name` Component name (for logging) 




        

<hr>



### function is\_initialized 

_Check if PIO is initialized._ 
```C++
static inline bool emulator::EmulatorIO::is_initialized () 
```





**Returns:**

true if [**initialize()**](classemulator_1_1EmulatorIO.md#function-initialize) has been called successfully 





        

<hr>



### function open\_file 

_Open an existing NetCDF file for reading._ 
```C++
static int emulator::EmulatorIO::open_file (
    const std::string & filename
) 
```





**Parameters:**


* `filename` Path to the file 



**Returns:**

NetCDF file ID (ncid), or -1 on failure 





        

<hr>



### function read\_var\_1d 

_Read a 1D double array from a NetCDF variable._ 
```C++
static bool emulator::EmulatorIO::read_var_1d (
    int ncid,
    const std::string & varname,
    double * data,
    int size
) 
```





**Parameters:**


* `ncid` File ID 
* `varname` Variable name 
* `data` Output buffer (must be pre-allocated) 
* `size` Number of elements to read 



**Returns:**

true on success, false if variable not found or read failed 





        

<hr>



### function read\_var\_1d\_int 

_Read a 1D integer array from a NetCDF variable._ 
```C++
static bool emulator::EmulatorIO::read_var_1d_int (
    int ncid,
    const std::string & varname,
    int * data,
    int size
) 
```





**Parameters:**


* `ncid` File ID 
* `varname` Variable name 
* `data` Output buffer 
* `size` Number of elements to read 



**Returns:**

true on success 





        

<hr>



### function read\_var\_2d 

_Read a 2D double array from a NetCDF variable._ 
```C++
static bool emulator::EmulatorIO::read_var_2d (
    int ncid,
    const std::string & varname,
    double * data,
    int nx,
    int ny
) 
```





**Parameters:**


* `ncid` File ID 
* `varname` Variable name 
* `data` Output buffer (must be pre-allocated, row-major) 
* `nx` Size in x dimension 
* `ny` Size in y dimension 



**Returns:**

true on success 





        

<hr>



### function read\_var\_3d\_slice 

_Read a 2D slice from a 3D variable (time, lat, lon)._ 
```C++
static bool emulator::EmulatorIO::read_var_3d_slice (
    int ncid,
    const std::string & varname,
    double * data,
    int nx,
    int ny,
    int time_idx
) 
```





**Parameters:**


* `ncid` File ID 
* `varname` Variable name 
* `data` Output buffer 
* `nx` Size in x dimension 
* `ny` Size in y dimension 
* `time_idx` Time index to read 



**Returns:**

true on success 





        

<hr>



### function sync\_file 

_Synchronize (flush) file contents to disk._ 
```C++
static void emulator::EmulatorIO::sync_file (
    int ncid
) 
```





**Parameters:**


* `ncid` File ID 




        

<hr>



### function write\_var\_1d 

_Write a 1D double array to a NetCDF variable._ 
```C++
static bool emulator::EmulatorIO::write_var_1d (
    int ncid,
    const std::string & varname,
    const double * data,
    int size
) 
```





**Parameters:**


* `ncid` File ID 
* `varname` Variable name (must already exist) 
* `data` Input data buffer 
* `size` Number of elements to write 



**Returns:**

true on success 





        

<hr>



### function write\_var\_2d 

_Write a 2D double array to a NetCDF variable._ 
```C++
static bool emulator::EmulatorIO::write_var_2d (
    int ncid,
    const std::string & varname,
    const double * data,
    int nx,
    int ny
) 
```





**Parameters:**


* `ncid` File ID 
* `varname` Variable name (must already exist) 
* `data` Input data buffer (row-major) 
* `nx` Size in x dimension 
* `ny` Size in y dimension 



**Returns:**

true on success 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_io.hpp`

