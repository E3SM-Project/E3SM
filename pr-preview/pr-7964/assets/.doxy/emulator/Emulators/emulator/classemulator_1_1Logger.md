

# Class emulator::Logger



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**Logger**](classemulator_1_1Logger.md)



_Simple logger with optional file output._ [More...](#detailed-description)

* `#include <emulator_logger.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Logger**](#function-logger) () = default<br> |
|  void | [**debug**](#function-debug) (const std::string & message) <br>_Log a debug message._  |
|  void | [**error**](#function-error) (const std::string & message) <br>_Log an error message._  |
|  void | [**info**](#function-info) (const std::string & message) <br>_Log an informational message._  |
|  void | [**log**](#function-log) ([**LogLevel**](namespaceemulator.md#enum-loglevel) level, const std::string & message) <br>_Log a message at the specified level._  |
|  void | [**set\_file**](#function-set_file) (const std::string & filename) <br>_Set the log file path._  |
|  void | [**verbose**](#function-verbose) (const std::string & message) <br>_Log a verbose message._  |
|  void | [**warn**](#function-warn) (const std::string & message) <br>_Log a warning message._  |
|   | [**~Logger**](#function-logger) () <br> |




























## Detailed Description


## Public Functions Documentation




### function Logger 

```C++
emulator::Logger::Logger () = default
```




<hr>



### function debug 

_Log a debug message._ 
```C++
void emulator::Logger::debug (
    const std::string & message
) 
```




<hr>



### function error 

_Log an error message._ 
```C++
void emulator::Logger::error (
    const std::string & message
) 
```




<hr>



### function info 

_Log an informational message._ 
```C++
void emulator::Logger::info (
    const std::string & message
) 
```




<hr>



### function log 

_Log a message at the specified level._ 
```C++
void emulator::Logger::log (
    LogLevel level,
    const std::string & message
) 
```





**Parameters:**


* `level` Log level 
* `message` Message to log 




        

<hr>



### function set\_file 

_Set the log file path._ 
```C++
void emulator::Logger::set_file (
    const std::string & filename
) 
```



If a file is already open, it is closed before opening the new file. If filename is empty, logging reverts to stdout.




**Parameters:**


* `filename` Path to log file (empty string for stdout) 




        

<hr>



### function verbose 

_Log a verbose message._ 
```C++
void emulator::Logger::verbose (
    const std::string & message
) 
```




<hr>



### function warn 

_Log a warning message._ 
```C++
void emulator::Logger::warn (
    const std::string & message
) 
```




<hr>



### function ~Logger 

```C++
emulator::Logger::~Logger () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_logger.hpp`

