(omega-dev-vertmix) =

# Vertical Mixing Coefficients

Omega includes a `VertMix` class that provides functions that compute `VertDiff` and `VertVisc`, the
vertical diffusivity and viscosity, where both are defined at the top of cell centers. Currently the
values of `VertDiff` and `VertVisc` are calculated using the combination of three options: (1) a
constant background mixing value, (2) a convective instability mixing value, and (3) a Richardson
number dependent shear mixing value based upon Pacanowski and Philander (1981). These options are
additive. For both the convective and shear mixing values `BruntVaisalaFreq` is needed, which
is calculated by the `EOS` class. Currently, the `VertMix` class is set up to have two enumerations,
`PP` and `KPP`, where the first only uses the three options listed above and the second uses the K Profile Parameterization [(KPP; Large et al., 1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872). The `KPP` option is not yet functional (it will throw an error if used), but the framework for its development have been added here.

## VertMix Type

An enumeration listing of the schemes is provided (`KPP` is not yet functional). It needs to be extended every time a `VertMix` type is added. It is used to identify which VertMix method is to be used at run time.

```c++
enum class VertMixType { PP, KPP };
```

## Initialization and Usage

The primary class `VertMix` is implemented using the Singleton pattern to ensure a single instance manages all vertical mixing calculations.

```c++
// Initialize VertMix
VertMix* vertMix = VertMix::getInstance("Default", mesh, nLevels);

// Compute mixing coefficients
vertMix->computeVertMix(VertDiff, VertVisc, NormalVelocity,
                       TangentialVelocity, BruntVaisalaFreq);
```

## Configuration

The initialization process reads parameters from the yaml configuration file with the following structure and
default values:

```yaml
VertMix:
  VertMixType: PP
  Background:
    Viscosity: 1e-4
    Diffusivity: 1e-5
  Convective:
    Enable: true
    Diffusivity: 1.0
    TriggerBVF: 0.0
  Shear:
    Enable: true
    NuZero: 0.005
    Alpha: 5.0
    Exponent: 2.0
```

## Class Structure

### Core Data Members

- `VertDiff`: 2D array storing vertical diffusivity coefficients (m²/s)
- `VertVisc`: 2D array storing vertical viscosity coefficients (m²/s)

### Mixing Parameters

1. Background Mixing:
   - `BackDiff`: Background vertical diffusivity (m²/s; Default: 1e-5)
   - `BackVisc`: Background vertical viscosity (m²/s: Default: 1e-4)

2. Convective Mixing:
   - `EnableConvMix`: Flag to enable/disable convective mixing (Default: True)
   - `ConvDiff`: Convective mixing coefficient (m²/s; Default: 1.0)
   - `ConvTriggerBVF`: Trigger threshold for convective mixing (Default: < 0.0)

3. Shear Mixing:
   - `EnableShearMix`: Flag to enable/disable shear mixing (Default: True)
   - `ShearNuZero`: Base coefficient for Pacanowski-Philander scheme (Default: 0.005)
   - `ShearAlpha`: Alpha parameter for P-P scheme (Default: 5.0)
   - `ShearExponent`: Exponent parameter for P-P scheme (Default: 2.0)

## Core Functionality (Vertical Mixing Coefficient Calculation)

The main computation is handled by:

```cpp
void computeVertMix(Array2DReal VertDiff,
                   Array2DReal VertVisc,
                   const Array2DReal &NormalVelocity,
                   const Array2DReal &TangentialVelocity,
                   const Array2DReal &BruntVaisalaFreq);
```

This method combines the effects of:
- Background mixing (constant coefficients)
- Convective mixing (triggered by static instability)
- Shear mixing (Pacanowski-Philander scheme)
