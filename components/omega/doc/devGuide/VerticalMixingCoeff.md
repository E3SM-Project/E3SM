(omega-dev-vertmix) =

# Vertical Mixing Coefficients

Omega includes a `VertMix` class that provides functions that compute `VertDiff` and `VertVisc`, the
vertical diffusivity and viscosity, where both are defined at the center of the cell and top of the layer.
Currently the values of `VertDiff` and `VertVisc` are calculated using the linear combination of three options: (1) a
constant background mixing value, (2) a convective instability mixing value, and (3) a Richardson
number dependent shear mixing value from the [Pacanowski and Philander (1981)](https://journals.ametsoc.org/view/journals/phoc/11/11/1520-0485_1981_011_1443_povmin_2_0_co_2.xml) parameterization. These options are linearly additive. In the future, additional additive options will be implemented, such as the K Profile Parameterization [(KPP; Large et al., 1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872). For both the convective and shear mixing values `BruntVaisalaFreqSq` is needed, which
is calculated by the `EOS` class.

## Initialization and Usage

The primary class `VertMix` is implemented using the Singleton pattern to ensure a single instance manages all vertical mixing calculations.

```c++
// Initialize VertMix
VertMix* VMix = VertMix::getInstance();

// Compute mixing coefficients
VMix->computeVertMix(NormalVelocity, TangentialVelocity, BruntVaisalaFreqSq);
```

## Configuration

The initialization process reads parameters from the yaml configuration file with the following structure and
default values:

```yaml
VertMix:
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
   - `ConvTriggerBVF`: Trigger threshold for convective mixing (Default: 0.0)

3. Shear Mixing:
   - `EnableShearMix`: Flag to enable/disable shear mixing (Default: True)
   - `ShearNuZero`: Base coefficient for Pacanowski-Philander scheme (Default: 0.005)
   - `ShearAlpha`: Alpha parameter for P-P scheme (Default: 5.0)
   - `ShearExponent`: Exponent parameter for P-P scheme (Default: 2.0)

## Core Functionality (Vertical Mixing Coefficient Calculation)

The main computation is handled by:

```cpp
void computeVertMix(const Array2DReal &NormalVelocity,
                   const Array2DReal &TangentialVelocity,
                   const Array2DReal &BruntVaisalaFreqSq);
```

This method combines the effects of:
- Background mixing (constant coefficients)
- Convective mixing (triggered by static instability)
- Shear instability driven mixing (Pacanowski-Philander scheme; to be changed to [Large et al., 1994](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872) shear mixing scheme in a later development)
