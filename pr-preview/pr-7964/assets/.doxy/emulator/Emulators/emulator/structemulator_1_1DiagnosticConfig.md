

# Struct emulator::DiagnosticConfig



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**DiagnosticConfig**](structemulator_1_1DiagnosticConfig.md)



_Complete diagnostic output configuration._ [More...](#detailed-description)

* `#include <emulator_diagnostics.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**HistoryRestartConfig**](structemulator_1_1HistoryRestartConfig.md) | [**history\_restart**](#variable-history_restart)  <br>_.rh. file config_  |
|  std::vector&lt; [**OutputStreamConfig**](structemulator_1_1OutputStreamConfig.md) &gt; | [**history\_streams**](#variable-history_streams)  <br>_.h. file streams_  |
|  [**RestartConfig**](structemulator_1_1RestartConfig.md) | [**restart**](#variable-restart)  <br>_.r. file config_  |












































## Detailed Description


## Public Attributes Documentation




### variable history\_restart 

_.rh. file config_ 
```C++
HistoryRestartConfig emulator::DiagnosticConfig::history_restart;
```




<hr>



### variable history\_streams 

_.h. file streams_ 
```C++
std::vector<OutputStreamConfig> emulator::DiagnosticConfig::history_streams;
```




<hr>



### variable restart 

_.r. file config_ 
```C++
RestartConfig emulator::DiagnosticConfig::restart;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_diagnostics.hpp`

