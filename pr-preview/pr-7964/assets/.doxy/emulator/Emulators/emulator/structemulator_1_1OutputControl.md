

# Struct emulator::OutputControl



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**OutputControl**](structemulator_1_1OutputControl.md)



_Controls output timing and tracks averaging state._ 

* `#include <emulator_output_stream.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  int | [**current\_step**](#variable-current_step)   = `0`<br> |
|  double | [**dt**](#variable-dt)   = `0.0`<br> |
|  int | [**frequency**](#variable-frequency)   = `1`<br> |
|  [**FrequencyUnit**](namespaceemulator.md#enum-frequencyunit) | [**frequency\_unit**](#variable-frequency_unit)   = `FrequencyUnit::NDAYS`<br> |
|  int | [**last\_write\_step**](#variable-last_write_step)   = `-1`<br> |
|  int | [**next\_write\_step**](#variable-next_write_step)   = `0`<br> |
|  int | [**nsamples\_since\_last\_write**](#variable-nsamples_since_last_write)   = `0`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**compute\_next\_write\_step**](#function-compute_next_write_step) (int current\_step, double dt) <br>_Compute next write step from frequency and current step._  |
|  bool | [**is\_write\_step**](#function-is_write_step) (int step) const<br>_Check if current step is a write step._  |
|  bool | [**output\_enabled**](#function-output_enabled) () const<br>_Check if output is enabled._  |
|  double | [**seconds\_per\_unit**](#function-seconds_per_unit) () const<br>_Get number of seconds per frequency unit._  |




























## Public Attributes Documentation




### variable current\_step 

```C++
int emulator::OutputControl::current_step;
```




<hr>



### variable dt 

```C++
double emulator::OutputControl::dt;
```




<hr>



### variable frequency 

```C++
int emulator::OutputControl::frequency;
```




<hr>



### variable frequency\_unit 

```C++
FrequencyUnit emulator::OutputControl::frequency_unit;
```




<hr>



### variable last\_write\_step 

```C++
int emulator::OutputControl::last_write_step;
```




<hr>



### variable next\_write\_step 

```C++
int emulator::OutputControl::next_write_step;
```




<hr>



### variable nsamples\_since\_last\_write 

```C++
int emulator::OutputControl::nsamples_since_last_write;
```




<hr>
## Public Functions Documentation




### function compute\_next\_write\_step 

_Compute next write step from frequency and current step._ 
```C++
void emulator::OutputControl::compute_next_write_step (
    int current_step,
    double dt
) 
```




<hr>



### function is\_write\_step 

_Check if current step is a write step._ 
```C++
bool emulator::OutputControl::is_write_step (
    int step
) const
```




<hr>



### function output\_enabled 

_Check if output is enabled._ 
```C++
inline bool emulator::OutputControl::output_enabled () const
```




<hr>



### function seconds\_per\_unit 

_Get number of seconds per frequency unit._ 
```C++
double emulator::OutputControl::seconds_per_unit () const
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_output_stream.hpp`

