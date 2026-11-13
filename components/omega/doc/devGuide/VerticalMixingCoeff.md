(omega-dev-vertmix) =

# Vertical Mixing Coefficients

Omega includes a `VertMix` class that provides functions that compute `VertDiff`, `VertVisc`, `GradRichNum`, and `GradRichNumSmoothed`, the
vertical diffusivity and viscosity, the gradient Richardson number, a smoothed gradient Richardson number, where all are defined at the center of the cell and top of the layer.
Currently the values of `VertDiff` and `VertVisc` are calculated using the linear combination of three options: (1) a
constant background mixing value, (2) a convective instability mixing value, and (3) a Richardson
number dependent shear instability driven mixing value from the [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/94RG01872) or LMD94 interior shear instability driven mixing parameterization. These options are linearly additive. In the future, additional additive options will be implemented, such as the K Profile Parameterization [(KPP; Large et al., 1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872). For both the convective and shear instability driven mixing values `BruntVaisalaFreqSq` is needed, which
is calculated by the `EOS` class. `GradRichNum` is smoothed using a 1-2-1 (vertical) filter to produce `GradRichNumSmoothed` which is used by the LMD94 shear instability driven mixing formulation.

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
    BaseShearValue: 0.005
    RiCrit: 0.7
    Exponent: 3.0
    RiSmoothLoops: 2
```

## Class Structure

### Core Data Members

- `VertDiff`: 2D array storing vertical diffusivity coefficients (m²/s)
- `VertVisc`: 2D array storing vertical viscosity coefficients (m²/s)
- `GradRichNum`: 2D array storing the gradient Richardson number
- `GradRichNumSmoothed`: 2D array storing the smoothed gradient Richardson number

### Mixing Parameters

1. Background Mixing:
   - `BackDiff`: Background vertical diffusivity (m²/s; Default: 1e-5)
   - `BackVisc`: Background vertical viscosity (m²/s: Default: 1e-4)

2. Convective Mixing:
   - `EnableConvMix`: Flag to enable/disable convective mixing (Default: True)
   - `ConvDiff`: Convective mixing coefficient (m²/s; Default: 1.0)
   - `ConvTriggerBVF`: Trigger threshold for convective mixing (Default: 0.0)

3. Shear Instability Driven Mixing:
   - `EnableShearMix`: Flag to enable/disable shear instability driven mixing (Default: True)
   - `BaseShearValue`: Base values of maximum viscosity and diffusivity for the LMD94 interior shear instability driven mixing formulation (Default: 0.005)
   - `RiCrit`: Critical Richardson number for the LMD94 formulation (Default: 0.7)
   - `ShearExponent`: Exponent parameter number for the LMD94 formulation (Default: 3.0)
   - `RiSmoothLoops`: Number of 1-2-1 filter passes to apply to the gradient Richardson number smoothing (Default: 2)

## Core Functionality (Vertical Mixing Coefficient Calculation)

The main computation is handled by:

```cpp
void computeVertMix(const Array2DReal &NormalVelocity,
                   const Array2DReal &TangentialVelocity,
                   const Array2DReal &BruntVaisalaFreqSq);
```
