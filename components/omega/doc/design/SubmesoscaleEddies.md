# Submesoscale Eddy Parameterization

## 1. Overview

This design document provides a summary of the evolution of the submesoscale eddy parameterization, commonly used in coarse-resolution ocean models to represent the effects of mixed-layer instabilities.

### The Fox-Kemper et al. (2008) Parameterization

The foundational parameterization for submesoscale eddies was developed by [Fox-Kemper et al. (2008)](https://doi.org/10.1175/2007JPO3792.1), hereafter FK08. It is designed to represent the restratification of the ocean mixed layer driven by mixed-layer instabilities (MLIs). These instabilities arise from horizontal buoyancy gradients within the mixed layer and act to flatten isopycnals, thereby converting available potential energy into eddy kinetic energy and stratifying the water column.

The core of the FK08 scheme is the introduction of an eddy-induced velocity, $\vec{v}_e$, which advects buoyancy. The vertical component of the eddy-induced overturning streamfunction, $\Psi_e$, is parameterized as:

$$
\Psi_e = C_e \frac{H^2}{|\vec{f}|} |\nabla_h b| \cdot S(z)
$$

Where:
- $C_e$ is a non-dimensional tuning coefficient, typically O(0.1).
- $H$ is the mixed-layer depth.
- $\vec{f}$ is the Coriolis parameter.
- $\nabla_h b$ is the horizontal buoyancy gradient.
- $S(z)$ is a shape function that defines the vertical structure of the overturning.

This streamfunction produces an eddy-induced vertical velocity, $w_e$, that restratifies the mixed layer by slumping the horizontal density fronts. The shape function $S(z)$ is typically chosen to be zero at the surface and the base of the mixed layer, with a maximum in between, ensuring conservation of mass.

### Implementation in Ocean Models

The implementation of the submesoscale eddy parameterization in coarse-resolution models like MOM6 is detailed in [Fox-Kemper et al. (2011)](https://doi.org/10.1016/j.ocemod.2010.09.002) (hereafter FK11). This version adapts the original theory to account for the fact that global models do not resolve the small-scale fronts where instabilities grow.

#### Parameterization Summary

The parameterization is centered on an eddy-induced overturning streamfunction, $\Psi$. For coarse-resolution models, this is formulated as (based on Eq. 6 of [Fox-Kemper et al., 2011](https://doi.org/10.1016/j.ocemod.2010.09.002)):

$$
\Psi = C_e \frac{\Delta_s}{L_f} \frac{H^2 \nabla_h \overline{b}^z}{\sqrt{f^2 + \tau^{-2}}} \mu(z)
$$ (submeso-definition)

The key parameters and variables are:
- $\Psi$: The eddy-induced overturning streamfunction.
- $C_e$: A non-dimensional efficiency coefficient, typically hard-coded to a value like 0.0625.
- $\Gamma_{\Delta}$: An additional tuning factor, often combined with other terms into a single coefficient.
- $\Delta_s$: The minimum of the local grid scale and the deformation radius.
- $L_f$: The characteristic width of submesoscale fronts, often treated as a constant or parameterized.
- $H$: The mixed-layer depth, often time-filtered to represent the memory of the instabilities.
- $|\nabla_h \bar{b}|$: The magnitude of the horizontal buoyancy gradient, averaged over the mixed layer.
- $f$: The Coriolis parameter.
- $\tau$: A timescale for momentum mixing across the mixed layer, often calculated from the surface friction velocity $u_*$.
- $\mu(z)$: A vertical shape function that is zero at the surface and the base of the mixed layer.

In practice, the upscaling factor $\frac{\Delta_s}{l_f}$ can be absorbed into a single tunable parameter, `FOX_KEMPER_ML_RESTRAT` in MOM6, simplifying the equation. The term $\sqrt{f^2 + \tau^{-2}}$ in the denominator provides a continuous and well-behaved formulation near the equator.

### Modifications by Bodner et al.

Subsequent studies revealed that the original FK08 parameterization could produce excessive and unrealistic restratification in certain regimes, particularly in regions of deep winter mixed layers and weak horizontal buoyancy gradients, such as the Southern Ocean.

[Bodner et al. (2023)](https://doi.org/10.1175/JPO-D-21-0297.1) (hereafter B23) proposed modifications to address this issue. Their primary change was to make the efficiency coefficient, $C_e$, a dynamic variable dependent on the local dynamical state of the mixed layer. The goal was to activate the parameterization only when the conditions are favorable for the growth of MLIs.

The modified coefficient, $C_e'$, is scaled by a factor that depends on the ratio of the mixed-layer Rossby number ($Ro = U/fL$) to a critical Rossby number. This effectively dampens the parameterization's strength when the vertical shear is weak compared to the planetary vorticity and stratification, preventing spurious restratification in dynamically quiescent regions.

### Impacts of Modifications (Uchida et al., 2024 Preprint)

A recent preprint by [Uchida et al. (2024)](https://essopenarchive.org/users/846346/articles/1298145-representation-of-surface-mixed-layer-eddies-affects-the-large-scale-ventilation-of-the-global-ocean) investigates the global impacts of the Bodner et al. modifications within a climate model. Their findings highlight significant improvements in model fidelity:

1.  **Reduced Spurious Restratification:** The modified scheme substantially reduces the excessive restratification seen with the original FK08 scheme, particularly in the high latitudes of the Southern Hemisphere during winter.
2.  **Improved Mixed-Layer Depths:** As a result, the model simulates deeper and more realistic winter mixed-layer depths in the Southern Ocean, which is a long-standing bias in many climate models.
3.  **Impacts on Water Mass Formation:** The deeper mixed layers enhance the formation of Subantarctic Mode Water (SAMW), improving the representation of water mass properties and the transport of heat and carbon into the ocean interior.
4.  **Teleconnections:** The changes in the Southern Ocean have far-reaching impacts, including a strengthening of the Atlantic Meridional Overturning Circulation (AMOC) due to modified water mass pathways.

Overall, the Uchida et al. study demonstrates that the physically-based constraints introduced by Bodner et al. lead to a more realistic and stable climate simulation.

### Eddy-Induced Velocity Shape Function

The parameterization assumes a specific vertical structure for the eddy-induced overturning, encapsulated in the shape function $S(z)$. This function dictates the profile of the eddy-induced vertical velocity, $w_e$, within the mixed layer (from the surface at $z=0$ to the base at $z=-H$). A common choice is a parabolic profile, which is zero at the boundaries and maximal in the middle of the mixed layer. This profile creates convergence of the eddy flow in the upper mixed layer and divergence in the lower mixed layer, driving restratification.

The schematic below illustrates the vertical profile of the eddy-induced streamfunction $\Psi_e$ and the resulting vertical velocity $w_e$.

```text
+----------------------------------------+
|              Mixed Layer               |
|                                        |
| z = 0 (Surface)                        |
|   |- w_e = 0 ---------> max Ψ_e        |
|   |               at z = -H/2          |
|   |- w_e > 0 (Upwelling ⇑)             |
|   v                                    |
| z = -H (ML Base)                       |
|   L- w_e < 0 (Downwelling ⇓)           |
+----------------------------------------+
```

## 2. Requirements

### 2.1 Requirement: Buoyancy calculation must be modular

To allow for easier addition of a future mesoscale eddy closure ([e.g. Mak et al 2018](https://journals.ametsoc.org/view/journals/phoc/48/10/jpo-d-18-0017.1.xml) and/or [Korn et al 2018](https://www.sciencedirect.com/science/article/pii/S1463500318301859)) the horizontal buoyancy calculation is required for the submesoscale eddy parameterization ([](#submeso-definition)) and for mesoscale eddy closures.  The implementation herein should be easily interfaced with each parameterization.  The key difference is that the mesoscale eddy parameterization is based on buoyancy gradients at every level, whereas the submesoscale parameterization considers the mixed layer average buoyancy.  For a non-boussinesq model, the buoyancy is defined as

$$
b \equiv g \frac{\rho - \rho_0}{\rho}
$$

where $\rho_0$ is a constant reference density.  The vertical mixed layer average ($H$) in pseudo-height coordinates is given by

$$
b_{ml} = \frac{1}{H} \int_{-H}^0 g\frac{\rho - \rho_0}{\rho} d \tilde{z}
$$

This can be rewritten as

$$
b_{ml} = \frac{1}{H}\int_{-H}^0 g d\tilde{z} - \frac{1}{H}\int_{-H}^0 g \frac{\rho_0}{\rho} d \tilde{z}
$$

since $g$ is constant, we can write

$$
b_{ml} = g  - \frac{1}{H}\int_{-H}^0 g \frac{\rho_0}{\rho} d \tilde{z} = g - \frac{1}{H}\int_{-H}^0 g \rho_0 \alpha
$$

where $\alpha$ is the specific volume.

### 2.2 Requirement: Buoyancy calculation must account for tilted pressure coordinates

As in MPAS-Ocean, the horizontal buoyancy gradient calculation must supported tilted coordinates.  In the psuedo height coordinate for Omega, the horizontal gradient of buoyancy for a tilting surface ($\tilde{z}$) is

$$
\nabla b = \nabla_h b + \frac{\partial b}{\partial \tilde{z}} \nabla_h \tilde{z}
$$

where $\nabla_h$ is a gradient along a fixed model layer in Omega.

### 2.3 Requirement: The parameterization must be configurable in the traditiona FK11 form and the B23 update

The parameterization of submesoscale eddies must accept fixed parameters from the input Yaml file for things like the efficiency ($C_e$), the frontal width ($L_f$) and the inverse timescale ($\tau$) to be consistent with FK.  Via a configuration flag the parameterization should use the B23 form with dynamic predictions of these key variables.

## 3. Algorithmic Formulation

MAYBE move the FK stuff and B23 here?  Need to update B23 as it's too terse.
