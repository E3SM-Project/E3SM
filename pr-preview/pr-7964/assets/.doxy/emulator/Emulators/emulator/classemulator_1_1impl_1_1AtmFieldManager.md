

# Class emulator::impl::AtmFieldManager



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**impl**](namespaceemulator_1_1impl.md) **>** [**AtmFieldManager**](classemulator_1_1impl_1_1AtmFieldManager.md)



_Field storage container for atmosphere emulator._ [More...](#detailed-description)

* `#include <atm_field_manager.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; double &gt; | [**aldif**](#variable-aldif)  <br>_NIR diffuse albedo [-]._  |
|  std::vector&lt; double &gt; | [**aldir**](#variable-aldir)  <br>_NIR direct albedo [-]._  |
|  std::vector&lt; double &gt; | [**asdif**](#variable-asdif)  <br>_Visible diffuse albedo [-]._  |
|  std::vector&lt; double &gt; | [**asdir**](#variable-asdir)  <br>_Visible direct albedo [-]._  |
|  std::vector&lt; double &gt; | [**cflx**](#variable-cflx)  <br>_CO2 flux [kg/m²/s]._  |
|  std::vector&lt; double &gt; | [**dens**](#variable-dens)  <br>_Air density [kg/m³]._  |
|  std::map&lt; std::string, std::vector&lt; double &gt; &gt; | [**dynamic\_fields**](#variable-dynamic_fields)  <br>_Runtime-registered fields for configurable I/O._  |
|  std::vector&lt; double &gt; | [**icefrac**](#variable-icefrac)  <br>_Sea ice fraction [-]._  |
|  std::vector&lt; double &gt; | [**lhf**](#variable-lhf)  <br>_Latent heat flux [W/m²]._  |
|  std::vector&lt; double &gt; | [**lndfrac**](#variable-lndfrac)  <br>_Land fraction [-]._  |
|  std::vector&lt; double &gt; | [**lwdn**](#variable-lwdn)  <br>_Downward longwave radiation [W/m²]._  |
|  std::vector&lt; double &gt; | [**lwup**](#variable-lwup)  <br>_Upward longwave radiation [W/m²]._  |
|  std::vector&lt; double &gt; | [**net\_inputs**](#variable-net_inputs)  <br>_Packed input tensor for inference._  |
|  std::vector&lt; double &gt; | [**net\_outputs**](#variable-net_outputs)  <br>_Packed output tensor from inference._  |
|  std::vector&lt; double &gt; | [**ocnfrac**](#variable-ocnfrac)  <br>_Ocean fraction [-]._  |
|  std::vector&lt; double &gt; | [**pbot**](#variable-pbot)  <br>_Pressure at bottom level [Pa]._  |
|  std::vector&lt; double &gt; | [**pslv**](#variable-pslv)  <br>_Sea level pressure [Pa]._  |
|  std::vector&lt; double &gt; | [**ptem**](#variable-ptem)  <br>_Potential temperature [K]._  |
|  std::vector&lt; double &gt; | [**qref**](#variable-qref)  <br>_Reference specific humidity [kg/kg]._  |
|  std::vector&lt; double &gt; | [**rainc**](#variable-rainc)  <br>_Convective precipitation [kg/m²/s]._  |
|  std::vector&lt; double &gt; | [**rainl**](#variable-rainl)  <br>_Large-scale precipitation [kg/m²/s]._  |
|  std::vector&lt; double &gt; | [**shf**](#variable-shf)  <br>_Sensible heat flux [W/m²]._  |
|  std::vector&lt; double &gt; | [**shum**](#variable-shum)  <br>_Specific humidity [kg/kg]._  |
|  std::vector&lt; double &gt; | [**snowc**](#variable-snowc)  <br>_Convective snowfall [kg/m²/s]._  |
|  std::vector&lt; double &gt; | [**snowhice**](#variable-snowhice)  <br>_Snow height over ice [m]._  |
|  std::vector&lt; double &gt; | [**snowhland**](#variable-snowhland)  <br>_Snow height over land [m]._  |
|  std::vector&lt; double &gt; | [**snowl**](#variable-snowl)  <br>_Large-scale snowfall [kg/m²/s]._  |
|  std::vector&lt; double &gt; | [**sst**](#variable-sst)  <br>_Sea surface temperature [K]._  |
|  std::vector&lt; double &gt; | [**swndf**](#variable-swndf)  <br>_NIR diffuse shortwave [W/m²]._  |
|  std::vector&lt; double &gt; | [**swndr**](#variable-swndr)  <br>_NIR direct shortwave [W/m²]._  |
|  std::vector&lt; double &gt; | [**swnet**](#variable-swnet)  <br>_Net shortwave radiation [W/m²]._  |
|  std::vector&lt; double &gt; | [**swvdf**](#variable-swvdf)  <br>_Visible diffuse shortwave [W/m²]._  |
|  std::vector&lt; double &gt; | [**swvdr**](#variable-swvdr)  <br>_Visible direct shortwave [W/m²]._  |
|  std::vector&lt; double &gt; | [**tbot**](#variable-tbot)  <br>_Temperature at bottom level [K]._  |
|  std::vector&lt; double &gt; | [**tref**](#variable-tref)  <br>_Reference temperature [K]._  |
|  std::vector&lt; double &gt; | [**ts**](#variable-ts)  <br>_Surface temperature [K]._  |
|  std::vector&lt; double &gt; | [**u10**](#variable-u10)  <br>_10m wind speed [m/s]_  |
|  std::vector&lt; double &gt; | [**u10withgusts**](#variable-u10withgusts)  <br>_10m wind with gusts [m/s]_  |
|  std::vector&lt; double &gt; | [**ubot**](#variable-ubot)  <br>_Zonal wind at bottom level [m/s]._  |
|  std::vector&lt; double &gt; | [**vbot**](#variable-vbot)  <br>_Meridional wind at bottom level [m/s]._  |
|  std::vector&lt; double &gt; | [**wsx**](#variable-wsx)  <br>_Zonal wind stress [N/m²]._  |
|  std::vector&lt; double &gt; | [**wsy**](#variable-wsy)  <br>_Meridional wind stress [N/m²]._  |
|  std::vector&lt; double &gt; | [**zbot**](#variable-zbot)  <br>_Height at bottom level [m]._  |


## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  int | [**N\_INPUT\_CHANNELS**](#variable-n_input_channels)   = `39`<br>_Default number of input channels (legacy, use config instead)._  |
|  int | [**N\_OUTPUT\_CHANNELS**](#variable-n_output_channels)   = `44`<br>_Default number of output channels (legacy, use config instead)._  |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**AtmFieldManager**](#function-atmfieldmanager) () = default<br> |
|  void | [**allocate**](#function-allocate) (int ncols) <br>_Allocate all field vectors for given number of columns._  |
|  void | [**deallocate**](#function-deallocate) () <br>_Deallocate all field vectors._  |
|  std::vector&lt; double &gt; \* | [**get\_field\_ptr**](#function-get_field_ptr) (const std::string & name) <br>_Get pointer to a field vector by name._  |
|  bool | [**is\_allocated**](#function-is_allocated) () const<br>_Check if fields are allocated._  |
|  void | [**register\_dynamic\_field**](#function-register_dynamic_field) (const std::string & name) <br>_Register or create a dynamic field._  |
|  void | [**set\_defaults**](#function-set_defaults) (int ncols) <br>_Set default climatological values for export fields._  |
|   | [**~AtmFieldManager**](#function-atmfieldmanager) () = default<br> |




























## Detailed Description


## Public Attributes Documentation




### variable aldif 

_NIR diffuse albedo [-]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::aldif;
```




<hr>



### variable aldir 

_NIR direct albedo [-]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::aldir;
```




<hr>



### variable asdif 

_Visible diffuse albedo [-]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::asdif;
```




<hr>



### variable asdir 

_Visible direct albedo [-]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::asdir;
```




<hr>



### variable cflx 

_CO2 flux [kg/m²/s]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::cflx;
```




<hr>



### variable dens 

_Air density [kg/m³]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::dens;
```




<hr>



### variable dynamic\_fields 

_Runtime-registered fields for configurable I/O._ 
```C++
std::map<std::string, std::vector<double> > emulator::impl::AtmFieldManager::dynamic_fields;
```




<hr>



### variable icefrac 

_Sea ice fraction [-]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::icefrac;
```




<hr>



### variable lhf 

_Latent heat flux [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::lhf;
```




<hr>



### variable lndfrac 

_Land fraction [-]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::lndfrac;
```




<hr>



### variable lwdn 

_Downward longwave radiation [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::lwdn;
```




<hr>



### variable lwup 

_Upward longwave radiation [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::lwup;
```




<hr>



### variable net\_inputs 

_Packed input tensor for inference._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::net_inputs;
```




<hr>



### variable net\_outputs 

_Packed output tensor from inference._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::net_outputs;
```




<hr>



### variable ocnfrac 

_Ocean fraction [-]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::ocnfrac;
```




<hr>



### variable pbot 

_Pressure at bottom level [Pa]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::pbot;
```




<hr>



### variable pslv 

_Sea level pressure [Pa]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::pslv;
```




<hr>



### variable ptem 

_Potential temperature [K]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::ptem;
```




<hr>



### variable qref 

_Reference specific humidity [kg/kg]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::qref;
```




<hr>



### variable rainc 

_Convective precipitation [kg/m²/s]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::rainc;
```




<hr>



### variable rainl 

_Large-scale precipitation [kg/m²/s]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::rainl;
```




<hr>



### variable shf 

_Sensible heat flux [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::shf;
```




<hr>



### variable shum 

_Specific humidity [kg/kg]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::shum;
```




<hr>



### variable snowc 

_Convective snowfall [kg/m²/s]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::snowc;
```




<hr>



### variable snowhice 

_Snow height over ice [m]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::snowhice;
```




<hr>



### variable snowhland 

_Snow height over land [m]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::snowhland;
```




<hr>



### variable snowl 

_Large-scale snowfall [kg/m²/s]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::snowl;
```




<hr>



### variable sst 

_Sea surface temperature [K]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::sst;
```




<hr>



### variable swndf 

_NIR diffuse shortwave [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::swndf;
```




<hr>



### variable swndr 

_NIR direct shortwave [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::swndr;
```




<hr>



### variable swnet 

_Net shortwave radiation [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::swnet;
```




<hr>



### variable swvdf 

_Visible diffuse shortwave [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::swvdf;
```




<hr>



### variable swvdr 

_Visible direct shortwave [W/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::swvdr;
```




<hr>



### variable tbot 

_Temperature at bottom level [K]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::tbot;
```




<hr>



### variable tref 

_Reference temperature [K]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::tref;
```




<hr>



### variable ts 

_Surface temperature [K]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::ts;
```




<hr>



### variable u10 

_10m wind speed [m/s]_ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::u10;
```




<hr>



### variable u10withgusts 

_10m wind with gusts [m/s]_ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::u10withgusts;
```




<hr>



### variable ubot 

_Zonal wind at bottom level [m/s]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::ubot;
```




<hr>



### variable vbot 

_Meridional wind at bottom level [m/s]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::vbot;
```




<hr>



### variable wsx 

_Zonal wind stress [N/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::wsx;
```




<hr>



### variable wsy 

_Meridional wind stress [N/m²]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::wsy;
```




<hr>



### variable zbot 

_Height at bottom level [m]._ 
```C++
std::vector<double> emulator::impl::AtmFieldManager::zbot;
```




<hr>
## Public Static Attributes Documentation




### variable N\_INPUT\_CHANNELS 

_Default number of input channels (legacy, use config instead)._ 
```C++
int emulator::impl::AtmFieldManager::N_INPUT_CHANNELS;
```




<hr>



### variable N\_OUTPUT\_CHANNELS 

_Default number of output channels (legacy, use config instead)._ 
```C++
int emulator::impl::AtmFieldManager::N_OUTPUT_CHANNELS;
```




<hr>
## Public Functions Documentation




### function AtmFieldManager 

```C++
emulator::impl::AtmFieldManager::AtmFieldManager () = default
```




<hr>



### function allocate 

_Allocate all field vectors for given number of columns._ 
```C++
void emulator::impl::AtmFieldManager::allocate (
    int ncols
) 
```



Allocates import and export field vectors. The `net_inputs` and `net_outputs` vectors are cleared but not pre-allocated (they are sized dynamically in `prepare_inputs()`).




**Parameters:**


* `ncols` Number of local columns 




        

<hr>



### function deallocate 

_Deallocate all field vectors._ 
```C++
void emulator::impl::AtmFieldManager::deallocate () 
```




<hr>



### function get\_field\_ptr 

_Get pointer to a field vector by name._ 
```C++
std::vector< double > * emulator::impl::AtmFieldManager::get_field_ptr (
    const std::string & name
) 
```



Looks up the field in the hardcoded field map, then in dynamic fields.




**Parameters:**


* `name` Field name (e.g., "ts", "tbot", "pbot") 



**Returns:**

Pointer to the field vector, or nullptr if not found 





        

<hr>



### function is\_allocated 

_Check if fields are allocated._ 
```C++
inline bool emulator::impl::AtmFieldManager::is_allocated () const
```





**Returns:**

true if [**allocate()**](classemulator_1_1impl_1_1AtmFieldManager.md#function-allocate) has been called 





        

<hr>



### function register\_dynamic\_field 

_Register or create a dynamic field._ 
```C++
void emulator::impl::AtmFieldManager::register_dynamic_field (
    const std::string & name
) 
```



If the field doesn't exist (not in hardcoded or dynamic fields), creates a new entry in dynamic\_fields. If already allocated, resizes the new field to match ncols.




**Parameters:**


* `name` Field name to register 




        

<hr>



### function set\_defaults 

_Set default climatological values for export fields._ 
```C++
void emulator::impl::AtmFieldManager::set_defaults (
    int ncols
) 
```



Initializes export fields with reasonable default values for testing.




**Parameters:**


* `ncols` Number of columns 




        

<hr>



### function ~AtmFieldManager 

```C++
emulator::impl::AtmFieldManager::~AtmFieldManager () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/eatm/src/impl/atm_field_manager.hpp`

