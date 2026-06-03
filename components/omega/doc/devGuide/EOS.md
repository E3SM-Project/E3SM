(omega-dev-eos) =

# Equation of State (EOS)

Omega includes an `Eos` class that provides functions that compute `SpecVol`, `SpecVolDisplaced`,
and `BruntVaisalaFreqSq`. Current EOS options are a linear EOS, a constant EOS,
or an EOS computed using the TEOS-10 75 term expansion from
[Roquet et al. 2015](https://www.sciencedirect.com/science/article/pii/S1463500315000566).
If `SpecVolDisplaced` is calculated with the linear or constant EOS option,
it will be equal to `SpecVol` as there is no pressure/depth dependence for
those EOS options. For the constant EOS option, `SpecVol` is set to `1/RhoSw`
for all active cells/layers. `SpecVolDisplaced` computes specific volume
adiabatically displaced to `K + KDisp` (where `K` counted positive downward, ie `K+1` is one layer below `K`). Note: `SpecVol` must be calculated before `BruntVaisalaFreqSq`, as
`SpecVol` is an input for the `BruntVaisalaFreqSq` calculation. If the linear EOS option is used, then the `BruntVaisalaFreqSq`
is calculated using linear coefficients. If the TEOS-10 option is used, the `BruntVaisalaFreqSq` is calculated with non-linear
coefficients according to the [TEOS-10 toolbox](https://www.teos-10.org/software.htm). Note: two assumption for ease of computation and efficiency have been made
for the `BruntVaisalaFreqSq` TEOS-10 option that differ from how it is calculated in the TEOS-10 toolbox:
(1) gravity is assumed to be constant and not a function of depth and latitude, and (2) the interface value of the specific volume is
calculated as the average between two layer values, rather than being recalculated using the interface values of temperature,
salinity, and pressure. Both of these assumptions incur less than a 1% error.
For the constant EOS option, `BruntVaisalaFreqSq` is identically zero.

## Eos type

An enumeration listing all implemented schemes is provided. It needs to be extended every time an
EOS is added. It is used to identify which EOS method is to be used at run time.

```c++
enum class EosType { LinearEos, Teos10Eos, ConstantEos };
```

## Initialization

An instance of the `Eos` class requires a [`HorzMesh`](#omega-dev-horz-mesh), so the mesh class
and all of its dependencies need to be initialized before the `Eos` class can be. The static method:

```c++
OMEGA::Eos::init();
```

initializes the default `Eos`. A pointer to it can be retrieved at any time using:

```c++
OMEGA::Eos* DefEos = OMEGA::Eos::getInstance();
```

## Computation of Eos

To compute `SpecVol` for a particular set of temperature, salinity, and pressure arrays, do

```c++
Eos.computeSpecVol(ConsrvTemp, AbsSalinity, Pressure);
```

`SpecVolDisplaced` is calculated using local temperature and salinity values, but a pressure
value at `K + KDisp`. To compute `SpecVolDisplaced` for a particular set of temperature, salinity,
and pressure arrays and displaced vertical index level, do

```c++
Eos.computeSpecVolDisp(ConsrvTemp, AbsSalinity, Pressure, KDisp);
```

where `KDisp` is the number of `k` layers you want to displace each specific volume layer to.
For example, to displace each level to one below, set `KDisp = 1`.

To compute `BruntVaisalaFreqSq` for a particular set of temperature, salinity, pressure, and specific
volume arrays, do

```c++
Eos.computeBruntVaisalaFreqSq(ConservTemp, AbsSalinity, Pressure, SpecVol);
```

## Helper functions for conversion

The TEOS-10 implementation includes helper functions for temperature
conversions and freezing-point calculations.

To compute conservative freezing temperature from absolute salinity and
pressure, use

```c++
ComputeSpecVolTeos10.calcCtFreezing(Sa, P, SaturationFract);
```

This helper follows the TEOS-10 polynomial approximation used by
`gsw_ct_freezing_poly` in the GSW toolbox.

To convert Conservative Temperature to potential temperature through the EOS
interface, use

```c++
Eos.calcPtFromCt(Sa, Ct);
```

To convert potential temperature back to Conservative Temperature through the
EOS interface, use

```c++
Eos.calcCtFromPt(Sa, Pt);
```

For `EosType::Teos10Eos`, these wrappers dispatch to TEOS-10 helper formulas.
For non-TEOS options (`LinearEos` and `ConstantEos`), the wrappers return the
input temperature unchanged.

## Removal of Eos

To clear the Eos instance do:

```c++
OMEGA::Eos::destroyInstance();
```
