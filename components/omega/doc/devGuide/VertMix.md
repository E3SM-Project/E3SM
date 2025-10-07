(omega-dev-vertmix) =

## Vertical Mixing Coefficients

Omega includes an `VertMix` class that provides functions that compute `VertDiff`, `VertVisc`, and
`BruntVaisalaFreq`.


Current EOS options are a linear EOS or an EOS computed using the TEOS-10 75 term expansion from
Roquet et al. 2015. If `SpecVolDisplaced` is calculated with the linear EOS option, it will be equal
to `SpecVol` as there is no pressure/depth dependence for the linear EOS. `SpecVolDisplaced` has two
options: `"absolute"` computes specific volume adiabatically displaced to `KDisp` for all `k` values
and `"relative"` computes specific volume adiabatically displaced to `K + KDisp`.

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
OMEGA::Eos* DefEos = OMEGA::Eos::getDefault();
```
The create method:
```c++
OMEGA::Eos::create(Name, Mesh, NVertLevels);
```
allocates the `SpecVol` and `SpecVolDisplaced` arrays for the mesh and vertical level dimensions.

## Computation of Eos

To compute `SpecVol` for a particular set of temperature, salinity, and pressure arrays, do
```c++
Eos.computeSpecVol(SpecVol, ConsrvTemp, AbsSalinity, Pressure);
```
`SpecVolDisplaced` is calculated using local temperature and salinity values, but a pressure
value at `KDisp` or `K + KDisp` depending on whether `"absolute"` or `"relative"` is provided
for `DispType`. To compute `SpecVolDisplaced` for a particular set of temperature, salinity,
and pressure arrays and displaced vertical index level, do
```c++
Eos.computeSpecVolDisp(SpecVol, ConsrvTemp, AbsSalinity, Pressure, KDisp, DispType);
```
where `KDisp` is the vertical `k` index of the level you want to displace all of the specific
volume to if `DispType = "absolute"` and is the number of `k` levels you want to displace each
specific volume level to if `DispType = "relative"`. For example, to displace the entire specific
volume column to the surface, set `KDisp = 0` and `DispType = "absolute"`, and to displace each
level to one below, set `KDisp = 1` and `DispType = "relative"`.

## Removal of Eos

To erase a specific named Eos instance use `erase`
```c++
OMEGA::Eos::erase(Name);
```
To clear all instances do:
```c++
OMEGA::Eos::clear();
```
