(omega-dev-eos) =

# Equation of State (EOS)

Omega includes an `Eos` class that provides functions that compute `SpecVol`, `SpecVolDisplaced`,
and `BruntVaisalaFreq`. Current EOS options are a linear EOS or an EOS computed using the TEOS-10
75 term expansion from [Roquet et al. 2015](https://www.sciencedirect.com/science/article/pii/S1463500315000566).
If `SpecVolDisplaced` is calculated with the linear EOS option, it will be equal to `SpecVol` as there
is no pressure/depth dependence for the linear EOS. `SpecVolDisplaced` computes specific volume
adiabatically displaced to `K + KDisp`. Note: `SpecVol` must be calculated before `BruntVaisalaFreq`, as
`SpecVol` is an input for the `BruntVaisalaFreq` calculation.

## Eos type

An enumeration listing all implemented schemes is provided. It needs to be extended every time an
EOS is added. It is used to identify which EOS method is to be used at run time.

```c++
enum class EosType { Linear, Teos10Poly75t };
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

To compute `BruntVaisalaFreq` for a particular set of temperature, salinity, pressure, and specific
volume arrays, do

```c++
Eos.computeBruntVaisalaFreq(ConservTemp, AbsSalinity, Pressure, SpecVol);
```

## Removal of Eos

To clear the Eos instance do:

```c++
OMEGA::Eos::destroyInstance();
