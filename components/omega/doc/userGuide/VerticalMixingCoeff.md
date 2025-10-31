(omega-user-vertmix)=

# Vertical Mixing Coefficients

The vertical mixing module in Omega handles the parameterization of unresolved vertical mixing
processes in the ocean. It calculates vertical diffusivity and viscosity coefficients that
determine how properties (like momentum, heat, salt, and biogeochemical tracers) mix vertically
in the ocean model. Currently, Omega offers three different mixing processes within the water column: (1) a constant
background mixing value, (2) a convective instability mixing value, and (3) a Richardson number
dependent shear instability driven mixing value from the [Pacanowski and Philander (1981)](https://journals.ametsoc.org/view/journals/phoc/11/11/1520-0485_1981_011_1443_povmin_2_0_co_2.xml) parameterization. These are linearly additive and are describe a bit
more in detail below. Other mixing processes and parameterizations, such as the the K Profile Parameterization [(KPP; Large et al., 1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872) will be added in the future.

The user-configurable options are the following parameters in the yaml configuration file:

```yaml
VertMix:
  Background:
    Viscosity: 1.0e-4    # Background vertical viscosity (m²/s)
    Diffusivity: 1.0e-5  # Background vertical diffusivity (m²/s)
  Convective:
    Enable: true         # Enables the convective-induced mixing option
    Diffusivity: 1.0     # Convective mixing coefficient (m²/s)
    TriggerBVF: -1e-4    # Brunt-Vaisala frequency threshold
  Shear:
    Enable: true         # Enables the shear-instability driven mixing option
    NuZero: 1.0e-2       # Base viscosity coefficient
    Alpha: 5.0           # Stability parameter
    Exponent: 2          # Richardson number exponent
```

## Vertical Mixing Processes/Types

### 1. Background Mixing

A constant background mixing value that represents small-scale mixing processes not explicitly resolved by the model. Typically, this is chosen to represent low values of vertical mixing
happening in the ocean's stratified interior.

### 2. Convective Mixing

Enhanced convective adjustment mixing that occurs in statically unstable regions of the water column to parameterize convection and homogenize properties. In Omega this is mixing is defaulted to occur when the Brunt Vaisala Frequency is less than 0.0 (unstable), and is off when equal to (neutral) or greater than (stable) 0.0.

$$
\kappa =
\begin{cases}
\kappa_{b} + \kappa_{conv} \quad \text{ if } N^2 < N^2_{crit}\\
\kappa_{b} \quad \text{ if } N^2 \geq N^2_{crit}
\end{cases}
$$

This is different than some current implementations (i.e. in MPAS-Ocean and the CVMix package), where convective adjustment occurs both with unstable and neutral conditions ($N^2 \leq N^2_{crit}$).

### 3. Shear Mixing

Mixing induced by vertical velocity shear, implemented using the Pacanowski-Philander scheme, through the gradient Richardson number (ratio of buoyancy to shear).

$$
\nu = \frac{\nu_o}{(1+\alpha Ri)^n} + \nu_b\,,
$$

$$
\kappa = \frac{\nu}{(1+\alpha Ri)} + \kappa_b\,.
$$

where $Ri$ is defined as:

$$
Ri = \frac{N^2}{\left|\frac{\partial \mathbf{U}}{\partial z}\right|^2}\,,
$$
