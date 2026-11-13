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

## First derivatives of specific volume

The `Eos` class can also compute the first derivatives of the specific volume
with respect to conservative temperature, absolute salinity, and pressure,
together with the specific volume itself:

```c++
Eos.computeSpecVolAndDerivs(ConservTemp, AbsSalinity, Pressure);
```

`Pressure` is the relative pressure (gauge pressure in Pa) as elsewhere in
`Eos`, and the derivatives are returned per `degC`, per `(g/kg)`, and per `Pa`
respectively. Note the pressure derivative is per Pascal, not per decibar.

The results are stored in the `SpecVolDCt`, `SpecVolDSa` and `SpecVolDP`
members alongside `SpecVol`, and all three are registered as fields in the
`Eos` group so they can be written to a stream. Because `SpecVol` is computed
here as well, `computeSpecVolAndDerivs` replaces a call to `computeSpecVol`
rather than accompanying one; calling both would evaluate the equation of state
twice. The valid range of the derivative fields spans the full range of `Real`
rather than starting at zero, since the salinity derivative is negative
everywhere and the temperature derivative is negative in cold, nearly fresh
water.

The two methods are kept separate rather than always computing the derivatives
because the derivatives roughly double the TEOS-10 arithmetic per cell and
layer, and not every call needs them. `AuxiliaryState::computeMomVertAux` is
the only place that calls `computeSpecVol`; everything else consumes the
`Eos::SpecVol` array rather than recomputing it, including
`computeBruntVaisalaFreqSq`, which takes the specific volume as an argument.
That one call site is reached once per time stepper stage through
`computeMomAux` and `computeAll` in the tendency calculation, and once more
from `VertMix::VertMixImplicit`, which refreshes the pressure and specific
volume before the vertical mixing coefficients are formed.

A run using the higher-order pressure gradient therefore needs the derivatives
at every time step, and at those call sites `computeSpecVolAndDerivs` takes the
place of the `computeSpecVol` call that would otherwise be made, leaving one
evaluation of the equation of state where there was one before. Two things
still call for the plain `computeSpecVol`. First, `PressureGradType` is a
runtime option that defaults to `Centered`, so a run may never need the
derivatives at all. Second, even with the higher-order pressure gradient
selected, the `VertMix::VertMixImplicit` update feeds only
`computeGeomZHeight` and `computeBruntVaisalaFreqSq`, neither of which reads
the derivatives, so computing them there would be work that nothing consumes.
Which method to call is thus a decision for each call site, not one the `Eos`
class should make for it.

There is no displaced counterpart to `computeSpecVolAndDerivs`. The pressure
gradient needs the derivatives at the in-situ pressure of the layer, whereas
`computeSpecVolDisp` exists to evaluate the specific volume at the pressure of
a displaced layer. Nothing about the derivatives prevents an adiabatically
displaced version: it would take the same `KDisp` argument as
`computeSpecVolDisp` and evaluate the same coefficients at the displaced
pressure, with no new polynomial. It is left out here only because no caller
needs it yet, and it would mean carrying three more model-sized arrays and
fields.

All four values come from a single pass over the equation of state. For
`EosType::Teos10Eos` the derivatives are the analytic derivatives of the same
75-term polynomial used for the specific volume, evaluated at the same
normalized state, so no second call to the equation of state is made. The
pressure derivative reuses the pressure coefficients already assembled for the
specific volume; the temperature and salinity derivatives need coefficient sets
of their own but share the normalization and the square root. For
`EosType::LinearEos` the derivatives are `-DRhoDT` and `-DRhoDS` times the
square of the specific volume, with no pressure dependence, and for
`EosType::ConstantEos` all three vanish.

The thermal expansion and haline contraction coefficients used by the
`BruntVaisalaFreqSq` calculation are formed from these same derivatives,
`alpha = SpecVolDCt / SpecVol` and `beta = -SpecVolDSa / SpecVol`, so the
polynomial coefficients exist in only one place.

### A note on GSW-C

The GSW toolbox may be redistributed only without modification, so the
derivative routines in GSW-C are not ported or adapted here; they also could
not be called from a Kokkos device kernel. The implementation instead
differentiates the published Roquet et al. 2015 polynomial that `Teos10Eos`
already carries. GSW-C is used unmodified, through its public API, as an
independent check in the unit test.

That test compares against `gsw_specvol_first_derivatives` over a range of
states and finds agreement of order `1e-14` for the temperature and salinity
derivatives. The pressure derivative agrees only to about `2e-12`, and the
difference is on the GSW-C side: its `v_P` is evaluated from coefficients that
have been pre-multiplied by their pressure exponents and rounded, whereas the
Omega implementation differentiates the full-precision coefficients and matches
the exact derivative to roughly `1e-16`.

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
