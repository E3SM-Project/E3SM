(omega-user-vertmix)=

# Vertical Mixing Coefficients

The vertical mixing module in Omega handles the parameterization of unresolved vertical mixing
processes in the ocean. It calculates vertical diffusivity and viscosity coefficients that
determine how properties (like momentum, heat, salt, and biogeochemical tracers) mix vertically
in the ocean model. Currently, Omega offers three different mixing processes within the water column: (1) a constant
background mixing value, (2) a convective instability mixing value, and (3) a Richardson number-dependent shear-instability-driven mixing value from the [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/94RG01872) or LMD94 shear instability driven mixing parameterization. These are linearly additive and are describe a bit
more in detail below. Other mixing processes and parameterizations, such as the the K Profile Parameterization [(KPP; Large et al., 1994)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94rg01872) will be added in the future. In addition to diffusivity and viscosity coefficients, the vertical mixing module calculates the gradient Richardson number and smooths that gradient Richardson number using a 1-2-1 filter before using it in the shear instability driven mixing calculation.

The user-configurable options are the following parameters in the yaml configuration file:

```yaml
VertMix:
  Background:
    Viscosity: 1.0e-4    # Background vertical viscosity (m²/s)
    Diffusivity: 1.0e-5  # Background vertical diffusivity (m²/s)
  Convective:
    Enable: true         # Enables the convective-induced mixing option
    Diffusivity: 1.0     # Convective mixing coefficient (m²/s)
    TriggerBVF: 0.0      # Squared Brunt-Vaisala frequency threshold
  Shear:
    Enable: true         # Enables the shear-instability driven mixing option
    BaseShearValue: 0.005 # Base viscosity/diffusivity value
    RiCrit: 0.7          # Critical Richardson number
    Exponent: 3.0        # Richardson number exponent
    RiSmoothLoops: 2     # Number of Richardson number smoothing loops
```

## Vertical Mixing Processes/Types

### 1. Background Mixing

A constant background mixing value that represents small-scale mixing processes not explicitly resolved or modeled. Typically, this is chosen to represent low values of vertical mixing happening in the ocean's stratified interior and is assumed roughly equivalent to the globally averaged interior mixing from all sources (e.g., from internal wave breaking).

### 2. Convective Mixing

Enhanced convective adjustment mixing that occurs in statically unstable regions of the water column to parameterize convection and homogenize properties. In Omega this mixing is defaulted to occur when the squared Brunt Vaisala Frequency is less than 0.0 (unstable), and is off when equal to (neutral) or greater than (stable) 0.0. The convective diffusivity is thus added to the background diffusivity ($\kappa_b$) to form the total diffusivity coefficient $\kappa$.

$$
\kappa =
\begin{cases}
\kappa_b + \kappa_{conv} \quad \text{ if } N^2 < N^2_{crit}\\
\kappa_b \quad \text{ if } N^2 \geq N^2_{crit}
\end{cases}
$$

This is different than some current implementations (i.e. in MPAS-Ocean and the CVMix package), where convective adjustment occurs both with unstable and neutral conditions ($N^2 \leq N^2_{crit}$). $\kappa_{conv}$ and $N^2_{crit}$ are constant parameters set in the `VertMix` section of the yaml file (`Diffusivity` and `TriggerBVF` under the `Convective` header).

### 3. Shear-Instability-Driven Mixing

Mixing induced by vertical pseudo-velocity shear, implemented using the LMD94 scheme, through the gradient Richardson number (ratio of buoyancy to shear).

$$
\nu_{shear} = \kappa_{shear} =
\begin{cases}
\nu_o \quad \text{ if } Ri_g < 0\\
\nu_o \left[1 - \left( \frac{Ri_g}{Ri_{crit}} \right)^2 \right]^p \text{ if } 0 \leq Ri_g < Ri_{crit}\\
0.0 \quad \text{ if } Ri_{crit} \leq Ri_g
\end{cases}
$$

where $\nu_o$, $Ri_{crit}$, and $p$ are constant parameters set in the `VertMix` section of the yaml file (`BaseShearValue`, `RiCrit`, and `Exponent` under the `Shear` header). $Ri_g$ is defined as:

$$
Ri_g = \frac{N^2}{\left|\frac{\partial \mathbf{U}}{\partial z}\right|^2}\,,
$$

where $N^2$ is calculated by the EOS based on the ocean state and $\mathbf{U}$ is the magnitude of the horizontal velocity. $Ri_g$ is calculated by the vertical mixing module and then smoothed with a 1-2-1 (vertical) filter before being used to calculate the shear-instability-driven mixing. $N^2$, $\left| \frac{\partial \mathbf{U}}{\partial z}\right|^2$ and $Ri_g$ of `k` are all defined at the cell center, top interface of layer `k`. $N^2$, $\nu$ and $\kappa$ are set to zero at the surface.
