

# Struct emulator::impl::AtmCouplingIndices



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**impl**](namespaceemulator_1_1impl.md) **>** [**AtmCouplingIndices**](structemulator_1_1impl_1_1AtmCouplingIndices.md)



_Coupling field indices for atmosphere component._ [More...](#detailed-description)

* `#include <atm_coupling.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  int | [**Faxa\_lwdn**](#variable-faxa_lwdn)   = `-1`<br>_Downward longwave radiation [W/m²]._  |
|  int | [**Faxa\_rainc**](#variable-faxa_rainc)   = `-1`<br>_Convective precipitation [kg/m²/s]._  |
|  int | [**Faxa\_rainl**](#variable-faxa_rainl)   = `-1`<br>_Large-scale precipitation [kg/m²/s]._  |
|  int | [**Faxa\_snowc**](#variable-faxa_snowc)   = `-1`<br>_Convective snowfall [kg/m²/s]._  |
|  int | [**Faxa\_snowl**](#variable-faxa_snowl)   = `-1`<br>_Large-scale snowfall [kg/m²/s]._  |
|  int | [**Faxa\_swndf**](#variable-faxa_swndf)   = `-1`<br>_NIR diffuse shortwave [W/m²]._  |
|  int | [**Faxa\_swndr**](#variable-faxa_swndr)   = `-1`<br>_NIR direct shortwave [W/m²]._  |
|  int | [**Faxa\_swnet**](#variable-faxa_swnet)   = `-1`<br>_Net shortwave radiation [W/m²]._  |
|  int | [**Faxa\_swvdf**](#variable-faxa_swvdf)   = `-1`<br>_Visible diffuse shortwave [W/m²]._  |
|  int | [**Faxa\_swvdr**](#variable-faxa_swvdr)   = `-1`<br>_Visible direct shortwave [W/m²]._  |
|  int | [**Faxx\_evap**](#variable-faxx_evap)   = `-1`<br>_Evaporative flux [kg/m²/s]._  |
|  int | [**Faxx\_lat**](#variable-faxx_lat)   = `-1`<br>_Latent heat flux [W/m²]._  |
|  int | [**Faxx\_lwup**](#variable-faxx_lwup)   = `-1`<br>_Upward longwave radiation [W/m²]._  |
|  int | [**Faxx\_sen**](#variable-faxx_sen)   = `-1`<br>_Sensible heat flux [W/m²]._  |
|  int | [**Faxx\_taux**](#variable-faxx_taux)   = `-1`<br>_Zonal wind stress [N/m²]._  |
|  int | [**Faxx\_tauy**](#variable-faxx_tauy)   = `-1`<br>_Meridional wind stress [N/m²]._  |
|  int | [**Sa\_dens**](#variable-sa_dens)   = `-1`<br>_Air density [kg/m³]._  |
|  int | [**Sa\_pbot**](#variable-sa_pbot)   = `-1`<br>_Pressure at bottom level [Pa]._  |
|  int | [**Sa\_pslv**](#variable-sa_pslv)   = `-1`<br>_Sea level pressure [Pa]._  |
|  int | [**Sa\_ptem**](#variable-sa_ptem)   = `-1`<br>_Potential temperature [K]._  |
|  int | [**Sa\_shum**](#variable-sa_shum)   = `-1`<br>_Specific humidity [kg/kg]._  |
|  int | [**Sa\_tbot**](#variable-sa_tbot)   = `-1`<br>_Temperature at bottom level [K]._  |
|  int | [**Sa\_u**](#variable-sa_u)   = `-1`<br>_Zonal wind at bottom level [m/s]._  |
|  int | [**Sa\_v**](#variable-sa_v)   = `-1`<br>_Meridional wind at bottom level [m/s]._  |
|  int | [**Sa\_z**](#variable-sa_z)   = `-1`<br>_Atmospheric height at bottom level [m]._  |
|  int | [**Sf\_ifrac**](#variable-sf_ifrac)   = `-1`<br>_Sea ice fraction [-]._  |
|  int | [**Sf\_lfrac**](#variable-sf_lfrac)   = `-1`<br>_Land fraction [-]._  |
|  int | [**Sf\_ofrac**](#variable-sf_ofrac)   = `-1`<br>_Ocean fraction [-]._  |
|  int | [**Si\_snowh**](#variable-si_snowh)   = `-1`<br>_Snow height over ice [m]._  |
|  int | [**Sl\_snowh**](#variable-sl_snowh)   = `-1`<br>_Snow height over land [m]._  |
|  int | [**So\_t**](#variable-so_t)   = `-1`<br>_Ocean surface temperature [K]._  |
|  int | [**Sx\_anidf**](#variable-sx_anidf)   = `-1`<br>_NIR diffuse albedo [-]._  |
|  int | [**Sx\_anidr**](#variable-sx_anidr)   = `-1`<br>_NIR direct albedo [-]._  |
|  int | [**Sx\_avsdf**](#variable-sx_avsdf)   = `-1`<br>_Visible diffuse albedo [-]._  |
|  int | [**Sx\_avsdr**](#variable-sx_avsdr)   = `-1`<br>_Visible direct albedo [-]._  |
|  int | [**Sx\_qref**](#variable-sx_qref)   = `-1`<br>_Reference specific humidity [kg/kg]._  |
|  int | [**Sx\_t**](#variable-sx_t)   = `-1`<br>_Merged surface temperature [K]._  |
|  int | [**Sx\_tref**](#variable-sx_tref)   = `-1`<br>_Reference temperature [K]._  |
|  int | [**Sx\_u10**](#variable-sx_u10)   = `-1`<br>_10m wind speed [m/s]_  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**initialize**](#function-initialize) ([**CouplingFieldsBase**](classemulator_1_1CouplingFieldsBase.md) & fields) <br>_Initialize indices from coupling field base._  |




























## Detailed Description


## Public Attributes Documentation




### variable Faxa\_lwdn 

_Downward longwave radiation [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_lwdn;
```




<hr>



### variable Faxa\_rainc 

_Convective precipitation [kg/m²/s]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_rainc;
```




<hr>



### variable Faxa\_rainl 

_Large-scale precipitation [kg/m²/s]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_rainl;
```




<hr>



### variable Faxa\_snowc 

_Convective snowfall [kg/m²/s]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_snowc;
```




<hr>



### variable Faxa\_snowl 

_Large-scale snowfall [kg/m²/s]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_snowl;
```




<hr>



### variable Faxa\_swndf 

_NIR diffuse shortwave [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_swndf;
```




<hr>



### variable Faxa\_swndr 

_NIR direct shortwave [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_swndr;
```




<hr>



### variable Faxa\_swnet 

_Net shortwave radiation [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_swnet;
```




<hr>



### variable Faxa\_swvdf 

_Visible diffuse shortwave [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_swvdf;
```




<hr>



### variable Faxa\_swvdr 

_Visible direct shortwave [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxa_swvdr;
```




<hr>



### variable Faxx\_evap 

_Evaporative flux [kg/m²/s]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxx_evap;
```




<hr>



### variable Faxx\_lat 

_Latent heat flux [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxx_lat;
```




<hr>



### variable Faxx\_lwup 

_Upward longwave radiation [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxx_lwup;
```




<hr>



### variable Faxx\_sen 

_Sensible heat flux [W/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxx_sen;
```




<hr>



### variable Faxx\_taux 

_Zonal wind stress [N/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxx_taux;
```




<hr>



### variable Faxx\_tauy 

_Meridional wind stress [N/m²]._ 
```C++
int emulator::impl::AtmCouplingIndices::Faxx_tauy;
```




<hr>



### variable Sa\_dens 

_Air density [kg/m³]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_dens;
```




<hr>



### variable Sa\_pbot 

_Pressure at bottom level [Pa]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_pbot;
```




<hr>



### variable Sa\_pslv 

_Sea level pressure [Pa]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_pslv;
```




<hr>



### variable Sa\_ptem 

_Potential temperature [K]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_ptem;
```




<hr>



### variable Sa\_shum 

_Specific humidity [kg/kg]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_shum;
```




<hr>



### variable Sa\_tbot 

_Temperature at bottom level [K]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_tbot;
```




<hr>



### variable Sa\_u 

_Zonal wind at bottom level [m/s]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_u;
```




<hr>



### variable Sa\_v 

_Meridional wind at bottom level [m/s]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_v;
```




<hr>



### variable Sa\_z 

_Atmospheric height at bottom level [m]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sa_z;
```




<hr>



### variable Sf\_ifrac 

_Sea ice fraction [-]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sf_ifrac;
```




<hr>



### variable Sf\_lfrac 

_Land fraction [-]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sf_lfrac;
```




<hr>



### variable Sf\_ofrac 

_Ocean fraction [-]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sf_ofrac;
```




<hr>



### variable Si\_snowh 

_Snow height over ice [m]._ 
```C++
int emulator::impl::AtmCouplingIndices::Si_snowh;
```




<hr>



### variable Sl\_snowh 

_Snow height over land [m]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sl_snowh;
```




<hr>



### variable So\_t 

_Ocean surface temperature [K]._ 
```C++
int emulator::impl::AtmCouplingIndices::So_t;
```




<hr>



### variable Sx\_anidf 

_NIR diffuse albedo [-]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_anidf;
```




<hr>



### variable Sx\_anidr 

_NIR direct albedo [-]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_anidr;
```




<hr>



### variable Sx\_avsdf 

_Visible diffuse albedo [-]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_avsdf;
```




<hr>



### variable Sx\_avsdr 

_Visible direct albedo [-]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_avsdr;
```




<hr>



### variable Sx\_qref 

_Reference specific humidity [kg/kg]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_qref;
```




<hr>



### variable Sx\_t 

_Merged surface temperature [K]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_t;
```




<hr>



### variable Sx\_tref 

_Reference temperature [K]._ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_tref;
```




<hr>



### variable Sx\_u10 

_10m wind speed [m/s]_ 
```C++
int emulator::impl::AtmCouplingIndices::Sx_u10;
```




<hr>
## Public Functions Documentation




### function initialize 

_Initialize indices from coupling field base._ 
```C++
void emulator::impl::AtmCouplingIndices::initialize (
    CouplingFieldsBase & fields
) 
```



Looks up each field name in the [**CouplingFieldsBase**](classemulator_1_1CouplingFieldsBase.md) maps and stores the corresponding index.




**Parameters:**


* `fields` Initialized [**CouplingFieldsBase**](classemulator_1_1CouplingFieldsBase.md) with parsed field lists 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/eatm/src/impl/atm_coupling.hpp`

