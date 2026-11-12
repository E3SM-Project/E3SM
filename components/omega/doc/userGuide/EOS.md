(omega-user-eos)=

# Equation of State (EOS)

The equation of state (EOS) for the ocean describes the relationship between specific volume of seawater (in $\textrm{m}^3/\textrm{kg}$; the reciprocal of density) and temperature (in $^{\circ}\textrm{C}$), salinity (in $\textrm{g/kg}$), and pressure (in $\textrm{dbar}$). Through the hydrostatic balance (which relates density/specific volume gradients to pressure gradients), the equation of state provides a connection between active tracers (temperature and salinity) and the fluid dynamics.

Three choices of EOS are provided by Omega: a linear EOS, a constant EOS, and a TEOS-10 EOS. The linear EOS simplifies the relationship by excluding the influence of pressure and using constant expansion/contraction coefficients, making the specific volume a simple linear function of temperature and salinity. The constant EOS is intended for idealized/debug-style configurations and sets the specific volume to `1/RhoSw` for all active cells/layers (with `RhoSw` set to `1026.0` in `GlobalConstants.h`). However, both simplified options should only be used for simpler idealized test cases and are not appropriate for realistic ocean simulations. The TEOS-10 EOS is a 75-term polynomial expression from [Roquet et al. 2015](https://www.sciencedirect.com/science/article/pii/S1463500315000566) that approximates the [Thermodynamic Equation of Seawater 2010](https://www.teos-10.org/pubs/TEOS-10_Manual.pdf), but in a less complex and more computationally efficient manner, and is the preferred EOS for real ocean simulations in Omega.

The user-configurable option `EosType` can be set to `linear`, `teos10`, or `constant`. Parameters in the `Linear` subsection are used only when `EosType` is set to `linear`.

```yaml
Eos:
   EosType: teos10
   Linear:
      DRhoDT: -0.2
      DRhoDS: 0.8
      RhoT0S0: 1000.0
```

where `DRhoDT` is the thermal expansion coefficient ($\textrm{kg}/(\textrm{m}^3 \cdot ^{\circ}\textrm{C})$), `DRhoDS` is the saline contraction coefficient ($\textrm{kg}/\textrm{m}^3$), and `RhoT0S0` is the reference density at (T,S)=(0,0) (in $\textrm{kg}/\textrm{m}^3$).

In addition to `SpecVol`, the displaced specific volume `SpecVolDisplaced`, the squared Brunt-Vaisala frequency `BruntVaisalaFreqSq` and the first derivatives of specific volume `SpecVolDCt`, `SpecVolDSa` and `SpecVolDP` are also calculated by the EOS.

## TEOS-10 Helper Conversions

When using the `teos10` EOS option, helper methods are available for common
temperature conversions used by ocean thermodynamics workflows:

- `calcCtFreezing(Sa, P, SaturationFract)` computes the Conservative
   Temperature freezing point from absolute salinity and pressure.
- `calcPtFromCt(Sa, Ct)` converts Conservative Temperature to potential
   temperature.
- `calcCtFromPt(Sa, Pt)` converts potential temperature to Conservative
   Temperature.

Negative absolute salinity is clamped to 0 g/kg in the TEOS-10 temperature
conversion and freezing helpers, and in the linear EOS freezing calculation.
The clamp is local to these helpers; it does not modify tracer salinity or
apply to specific-volume and derivative calculations. This clamp is to avoid
producing NaNs in the tracer tendencies (teos-10 functions) or produce a
positive freezing temperature (linear case) that would trigger frazil formation
if the ocean state reaches negative salinity.

These helper methods are available through the EOS implementation but do not
replace the standard `computeSpecVol`, `computeSpecVolDisp`, or
`computeBruntVaisalaFreqSq` calculations.

## First Derivatives of Specific Volume

The `Eos` class can also compute the first derivatives of the specific volume with respect to conservative temperature (in $\textrm{m}^3\textrm{kg}^{-1}\,^{\circ}\textrm{C}^{-1}$), absolute salinity (in $\textrm{m}^3\textrm{g}^{-1}$), and pressure (in $\textrm{m}^3\textrm{kg}^{-1}\textrm{Pa}^{-1}$). These are needed by the higher-order horizontal pressure gradient, which expands the specific volume about a reference state within each layer instead of evaluating the full equation of state at every quadrature point.

There is no user-configurable option associated with the derivatives. They are computed on request by the parts of the model that need them, in the same pass that computes the specific volume, so selecting `teos10` does not make the model slower unless a scheme that uses them is enabled. They are available for all three `EosType` choices: for `teos10` they are the analytic derivatives of the same 75-term polynomial, for `linear` they follow from the configured `DRhoDT` and `DRhoDS` and have no pressure dependence, and for `constant` they are zero.

The derivatives are stored in the `SpecVolDCt`, `SpecVolDSa` and `SpecVolDP` fields of the `Eos` field group and can be requested in a stream's contents just like `SpecVol`.

The thermal expansion and haline contraction coefficients that enter the squared Brunt-Vaisala frequency are computed from these same derivatives.

## Displaced Specific Volume

The `Eos` class calculates the density of a parcel of fluid that is adiabatically displaced by a relative `k` levels (`k` counted positive downward), capturing the effects of pressure/depth changes. This is primarily used to calculate quantities for determining the water column stability (i.e. the stratification) and the vertical mixing coefficients (viscosity and diffusivity). Note: when using the `Linear` or `constant` EOS option, `SpecVolDisplaced` will be the same as `SpecVol` since the specific volume calculation is independent of pressure/depth.

## Squared Brunt Vaisala Frequency

The `Eos` class also calculates the squared Brunt Vaisala Frequency, which is used by the `VertMix` class to determine water column stability for both the convective adjustment and shear-instability driven mixing. When using the `Linear` EOS option, the `BruntVaisalaFreqSq` is calculated using linear coefficients. When using the `Teos-10` EOS option, the `BruntVaisalaFreqSq` is calculated with non-linear coefficients according to the TOES-10. When using the `constant` EOS option, `BruntVaisalaFreqSq` is zero in this implementation. For all options, `SpecVol` must be calculated prior to calculating `BruntVaisalaFreqSq`, as it is an input to the `BruntVaisalaFreqSq` calculation.
