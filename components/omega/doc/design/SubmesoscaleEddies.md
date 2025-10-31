(omega-design-submesoscale-eddies)=
# Submesoscale Eddy Parameterization

## 1. Overview

This design document provides a summary of the evolution of the submesoscale eddy parameterization, commonly used in coarse-resolution ocean models to represent the effects of mixed-layer instabilities.

### The Fox-Kemper et al. (2008) Parameterization

The foundational parameterization for submesoscale eddies was developed by [Fox-Kemper et al. (2008)](https://doi.org/10.1175/2007JPO3792.1), hereafter FK08. It is designed to represent the restratification of the ocean mixed layer driven by mixed-layer instabilities (MLIs). These instabilities arise from horizontal buoyancy gradients within the mixed layer and act to flatten isopycnals, thereby converting available potential energy into eddy kinetic energy and stratifying the water column.

The core of the FK08 scheme is the introduction of an eddy-induced velocity, $\vec{v}_e$, which advects buoyancy. The vertical component of the eddy-induced overturning streamfunction, $\Psi_e$, is parameterized as:

$$
\Psi_e = C_e \frac{H^2}{|\vec{f}|} |\nabla_h b| \cdot \mu(z)
$$

Where:
- $C_e$ is a non-dimensional tuning coefficient, typically O(0.1).
- $H$ is the mixed-layer depth.
- $\vec{f}$ is the Coriolis parameter.
- $\nabla_h b$ is the horizontal buoyancy gradient.
- $\mu(z)$ is a shape function that defines the vertical structure of the overturning.

This streamfunction produces an eddy-induced vertical velocity, $w_e$, that restratifies the mixed layer by slumping the horizontal density fronts. The shape function $\mu(z)$ is typically chosen to be zero at the surface and the base of the mixed layer, with a maximum in between, ensuring conservation of mass.  It's form is discussed in Section [](#33-eddy-induced-velocity-shape-function).

### Implementation in Ocean Models

The implementation of the submesoscale eddy parameterization in coarse-resolution models like MOM6 is detailed in [Fox-Kemper et al. (2011)](https://doi.org/10.1016/j.ocemod.2010.09.002) (hereafter FK11). This version adapts the original theory to account for the fact that global models do not resolve the small-scale fronts where instabilities grow.

#### Parameterization Summary

The parameterization is centered on an eddy-induced overturning streamfunction, $\Psi$. For coarse-resolution models, this is formulated as (based on Eq. 6 of [Fox-Kemper et al., 2011](https://doi.org/10.1016/j.ocemod.2010.09.002)):


$$
\Psi(z) = C_e \frac{\Delta_s}{L_f} \frac{H^2 \nabla_h \overline{b_{ml}\times \hat{z}}}{\sqrt{f^2 + \tau^{-2}}} \mu(z)
$$ (submeso-definition)

The key parameters and variables are:
- $\Psi$: The eddy-induced overturning streamfunction.
- $C_e$: A non-dimensional efficiency coefficient, typically hard-coded to a value like 0.06.
- $\Delta_s$: The model grid length parameter, defined in Section [](#31-the-original-fk11-formulation).
- $L_f$: The characteristic width of submesoscale fronts, often treated as a constant or parameterized.
- $H$: The mixed-layer depth, often time-filtered to represent the memory of the instabilities.
- $|\nabla_h \bar{b}|$: The magnitude of the horizontal buoyancy gradient, averaged over the mixed layer..
- $\mu(z)$: A vertical shape function that is zero at the surface and the base of the mixed layer.

### Modifications by Bodner et al.

Subsequent studies revealed that the original FK08 parameterization could produce excessive and unrealistic restratification in certain regimes, particularly in regions of deep winter mixed layers and weak horizontal buoyancy gradients, such as the Southern Ocean.

[Bodner et al. (2023)](https://doi.org/10.1175/JPO-D-21-0297.1) (hereafter B23) proposed modifications to address this issue. Their primary change was to make the parameters above dynamic. The secondary goal was to activate the parameterization only when the conditions are favorable for the growth of MLIs.  Further details are discussed in Section 3 below.

### Impacts of Modifications (Uchida et al., 2024 Preprint)

A recent preprint by [Uchida et al. (2024)](https://essopenarchive.org/users/846346/articles/1298145-representation-of-surface-mixed-layer-eddies-affects-the-large-scale-ventilation-of-the-global-ocean) investigates the global impacts of the Bodner et al. modifications within a climate model. Their findings highlight significant improvements in model fidelity:

1.  **Reduced Spurious Restratification:** The modified scheme substantially reduces the excessive restratification seen with the original FK08 scheme, particularly in the high latitudes of the Southern Hemisphere during winter.
2.  **Improved Mixed-Layer Depths:** As a result, the model simulates deeper and more realistic winter mixed-layer depths in the Southern Ocean, which is a long-standing bias in many climate models.
3.  **Impacts on Water Mass Formation:** The deeper mixed layers enhance the formation of Subantarctic Mode Water (SAMW), improving the representation of water mass properties and the transport of heat and carbon into the ocean interior.
4.  **Teleconnections:** The changes in the Southern Ocean have far-reaching impacts, including a strengthening of the Atlantic Meridional Overturning Circulation (AMOC) due to modified water mass pathways.

Overall, the Uchida et al. study demonstrates that the physically-based constraints introduced by Bodner et al. lead to a more realistic and stable climate simulation.

## 2. Requirements

### 2.1 Requirement: Buoyancy calculation must be modular

To allow for easier addition of a future mesoscale eddy closure ([e.g. Mak et al 2018](https://journals.ametsoc.org/view/journals/phoc/48/10/jpo-d-18-0017.1.xml) and/or [Korn et al 2018](https://www.sciencedirect.com/science/article/pii/S1463500318301859)) the horizontal buoyancy calculation is required for the submesoscale eddy parameterization ([](#submeso-definition)) and for mesoscale eddy closures.  The implementation herein should be easily interfaced with each parameterization.  For a non-boussinesq model, the buoyancy is defined as

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

In Omega two forms of the submesoscale eddy parameterization will be included: FK11 and B23.

### 3.1 The original FK11 formulation

The basic definition of the submesoscale parameterization is given in [](#submeso-definition).  In this equation, $L_f$, $C_e$, and $\tau$ are parameters.  Some models, such as MOM6 merge $\frac{\Delta_s}{L_f}$ into a single parameter.  This is not done in Omega to allow for variations in the grid size ($\Delta_s$) common to many target unstructured meshes.  As in MPAS-Ocean, the characteristic frontal width ($L_f$) is given by

$$
L_f = \text{max}\left(L_f^{min}, \frac{NH}{\tau}, \frac{\nabla b_{ml} H}{\tau^2} \right)
$$

where $L_f^{min}$ is a chosen minimum width, typically taken between 1 and 5 km and $\tau$ is an inverse timescale, which is defined as

$$
\tau = \sqrt{f^2 + \tau_{min}^{-2}}
$$

Here, $f$ is the Coriolis parameter and $\tau_{min}$ is a chosen minimum frontal slumping timescale to prevent excessive values of $\tau$ near the equator.

In the frontal width definition, the second length scale is the mixed layer deformation length scale.  Given prevalent model biases in simulations of $N$, especially in the mixed layer, the third length scale estimate was introduced.

Finally, the $\Delta_s$ parameter is taken as the minimum of `dcEdge` and `dsMax` where dsMax is a parameter usually set to 100km.

### 3.2 B23 updates

Leveraging turbulent thermal wind balance [McWilliams et al, 2015](https://journals.ametsoc.org/view/journals/phoc/45/8/jpo-d-14-0211.1.xml) and leveraging the ePBL mixing closure [Reichl and Hallberg, 2018](https://www.sciencedirect.com/science/article/pii/S1463500318301069), [](#submeso-definition) can be rewritten as

$$
\Psi(z) = C_r \frac{\Delta_s \left|f\right| h H^2 \nabla b_{ml} \times \hat{z} }{(m_* u_u^3 + n_* w_*^3)^{2/3}}\mu(z)
$$ (bodner23-definition)

While this formula introduces new parameters ($m_*$ and $n_*$) it has a number of important advantages over the FK11 formulation.  First, the ad-hoc $\tau_{min}$ factor to prevent infinite timescale values at the equator does not appear. Second, this form automatically tapers to zero at the equator ($f \rightarrow 0$), which is physically expected.  Finally, the need to compute an uncertain $N$ in the mixed layer is no longer necessary.  The two new parameters can be thought of as tunable parameters to match any parameterization of the turbulent momentum flux or can be determined through spatically dependent formulas such as those in [Reichl and Hallberg, 2018](https://www.sciencedirect.com/science/article/pii/S1463500318301069), their (25)-(28).  The parameters in [](#bodner23-definition) are defined in the following table

| Parameter | Definition |
|---|:---:|
| $h$ | boundary layer depth |
| $H$ | mixed layer depth |
| $u_*$ | surface friction velocity ($u_* = \sqrt{\tau/\rho}$) |
| $w_*$ | convective velocity scale ($w_* = \left(B_o h\right)^{1/3}$) |
| $\Delta_s$ | grid length limiter ($\Delta_s = \text{min}(dcEdge,100km)) |
| $b_{ml}$ | buoyancy averaged over the mixed layer depth |
| $m_*$, $n_*$ | tunable parameters to match turbulent momentum flux |
| $C_r$ | efficiency parameter $C_r \equiv C_e / Ri_{crit}$ |
| $Ri_{crit}$ | critical Richardson number, commonly set to 0.25 |


### 3.3 Eddy-Induced Velocity Shape Function

Both the FK11 and B23 parameterizations assume the same vertical structure for the eddy-induced overturning, encapsulated in the shape function $\mu(z)$. This function dictates the profile of the eddy-induced vertical velocity, $w_e$, within the mixed layer (from the surface at $z=0$ to the base at $z=-H$). A common choice is a parabolic profile, which is zero at the boundaries and maximal in the middle of the mixed layer. This profile creates convergence of the eddy flow in the upper mixed layer and divergence in the lower mixed layer, driving restratification.

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

The function is given as

$$
\mu(z) = \text{max}\left\{0,\left[1-\left(\frac{2z}{H}+1\right)^2\right]\left[1+\frac{5}{21}\left(\frac{2z}{H}+1\right)^2\right]\right\}
$$

### 3.4 Mixed layer depth calculation

Both the FK11 and B23 formulations are critically dependent on determination of the mixed layer depth. For this parameterization, the mixed layer depth will use the density threshold criterion ([Holte and Talley, 2009](https://journals.ametsoc.org/view/journals/atot/26/9/2009jtecho543_1.xml), their option 1).  The mixed layer depth is defined as the pressure where

$$
\left|\rho(p) - \rho(p_0)\right| \geq \Delta \rho_t
$$ (mld-definition)

In [](#mld-definition), $p_0$ is a reference pressure, which is set to $10^5$ Pa ($\approx 10$m), and $\Delta \rho_t$ is a density threshold criterion, which is typically set to 0.03 kg m$^{-3}$.  When the density threshold occurs between layers, the mixed layer depth is found via a linear fit between the bounding layers.

In future versions of Omega, ice shelf cavities slightly complicate the determination of the mixed layer depth.  The pressure must be referenced to the land ice pressure.  This document will be updated once ice shelf cavities are implemented in Omega.

## 4. Design

The most complex part of this design is the spatial dependence implicit in the averaging of the buoyancy.  The necessary `parallelReduce` calculation requires reduction over a spatially variable inner range.  This will be accomplished with the hierarchical parallelism wrappers described in CITE HP design doc.  In the initial implementation, only the FK11 form will be included.  The B23 formulation requires the boundary layer depth, which can not be included until KPP [Large et al, 1994](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/94RG01872) is implemented for vertical mixing.

### 4.1 Data types and parameters

The following options are necessary for the submesoscale eddy parameterization:
- `EnableBodner23` if true, the modifications from Bodner et al. (2023) are utilized, if false, the original FK11 form is used.
- `LfMin` controls the lower limit for the diagnosed frontal width (500m to 5km)
- `TauMin` minimum timescale for spindown
- `Mstar` nondimensional parameter controlling the contribution of surface friction velocity in the B23 form
- `Nstar` nondimensional parameter controlling the contribution of convective velocity scale in the B23 form
- `Ce` submesoscale eddy efficiency (0.06-0.08)
- `RiCrit` critical Richardson number for the B23 form (~0.25)
- `DsMax` maximum grid length to include in the submesoscale calculation (~100 km)

### 4.2 Methods
