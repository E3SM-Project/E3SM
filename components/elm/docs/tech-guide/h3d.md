# Hybrid-3D hillslope hydrological model

## Overview
H3D represents subsurface lateral flow and groundwater dynamics within an idealized hillslope unit inside each land 
grid cell. Rather than treating each soil column as an isolated vertical profile, 
h3D connects multiple columns along a topographic gradient to explicitly simulate downslope drainage, 
water table redistribution, and interaction with the channel.

This approach is conceptually one-dimensional along the hillslope, 
but retains three-dimensional realism by including slope, width, and drainage area variation.
The model solves the Dupuit–Boussinesq groundwater flow equation for saturated thickness using an implicit finite-difference method.

## Conceptual Representation

Each land unit is subdivided into several h3D columns positioned along an idealized hillslope:

- Lower boundary corresponds to the stream outlet.
- Upper boundary represents the local topographic divide.
- Intermediate nodes represent soil columns along the slope.

For each node, the state variable is the saturated thickness, $h_{sat}(x,t)$, measured from the bedrock to the water table. 
The model tracks how $h_{sat}(x,t)$ evolves in time due to:

- Local recharge (vertical infiltration)
- Downslope lateral flow driven by topographic gradients
- Variable transmissivity and drainable porosity along the slope


Each land unit represents a single hillslope, which has consistent topographic and geometric properties:
- **Overall slope angle:** $\theta$ — mean hillslope angle (rad)
- **Width function:** $w(x)$ — lateral width distribution along the hillslope (m)
- **Distance function:** $x(i)$ — distance from the stream outlet to node $i$ (m)
- **Total hillslope area:** $A_{hs} = \displaystyle \int_0^L w(x)\,dx$ (m²)

Each column is a point along hillslope, representing a cross-section at distance x along the hillslope. Each has local properties, connected by lateral flow to adjacent columns:
- Soil hydraulic properties (K_sat, porosity, etc.)
- Local area contribution (hs_dA)
- Bedrock depth

Land unit captures the hillsclope-scale connectivity, while columns capture spatial variation along the flow path. Each column contributes proportionally to the area within the land unit, conserves water mass at hillslope scale, having realistic flow convergence/divergence. Column 1 = stream/outlet, column N = hillslope divide.


## Govergning Equation

The fundamental PDE is a Dupuit-style Boussinesq groundwater flow equation for saturated flow along the slope:

$$
\frac{\partial h}{\partial t} = \frac{1}{f_{\text{drain}}} \frac{\partial}{\partial x} \left[ T(x,h)\left(\frac{\partial h}{\partial }\cos\theta + \sin\theta\right)\right]+ \frac{R}{f_{\text{drain}}}
$$

where $h(x,t)$ is the saturated thickness [m],
$f_{\text{drain}}$ is the drainable porosity [–],
$T(x,h)$ is the transmissivity [m² s⁻¹],
$\theta$ is the hillslope angle [rad],
and $R$ is the recharge rate [m s⁻¹].

## Boundary Conditions

| Boundary       | Condition                             | Implementation           |
| :------------- | :------------------------------------ | :----------------------- |
| Lower (stream) | Head-dependent outflow                | Robin-type term in $r_1$ |
| Upper (divide) | No-flow ($\partial h/\partial x = 0$) | Set $c_N = 0$            |

## Constitutive Relationships

### Transmissivity

The local transmissivity depends on saturated thickness and anisotropy:

$$
T(x,h)= \frac{K_{\text{aniso}}\,K_{\text{sat}}(z_{wt})\,h\,w(x)}{1000}
$$

where $K_{\text{aniso}}=100$ is the horizontal/vertical anisotropy factor,
$K_{\text{sat}}(z_{wt})$ is the saturated hydraulic conductivity at the
water-table depth [mm s⁻¹], and $w(x)$ is the hillslope width [m].
Division by 1000 converts from mm s⁻¹ to m s⁻¹.

### Variable Drainable Porosity

The specific yield varies with depth following a Brooks–Corey relation:

$$
f_{\text{drain}}
= \theta_{\text{sat}}
\left[
1 -
\left(
1 +
\frac{1000\,\max(0,\,z_{\text{bed}}-h)}{\psi_{\text{sat}}}
\right)^{-1/b}
\right],
\qquad
f_{\text{drain}}\ge0.02
$$

where $\theta_{\text{sat}}$ is porosity [–],
$z_{\text{bed}}$ is bedrock depth [m],
$\psi_{\text{sat}}$ is air-entry suction [mm],
and $b$ is the Brooks–Corey pore-size index [–].

This function allows for a smooth transition between unsaturated and fully saturated conditions and ensures stability under variable soil thickness.


## Numerical Implementation

### Spatial Discretization

The PDE is solved implicitly in space and time using a tridiagonal
system for $h_i^{n+1}$ at each node $i$:

$$
a_i h_{i-1}^{n+1} + b_i h_i^{n+1} + c_i h_{i+1}^{n+1} = r_i
$$


Interior nodes ($i=2,\dots,N-1$)

$$
\begin{aligned}
a_i &= -\frac{T_{i-\frac12}^n \cos\theta \,\Delta t}
           {\Delta x_{i-\frac12}\,\Delta x_i\,w_i}, \\
c_i &= -\frac{T_{i+\frac12}^n \cos\theta \,\Delta t}
           {\Delta x_{i+\frac12}\,\Delta x_i\,w_i}, \\
b_i &= f_{\text{drain},i} - (a_i + c_i), \\
r_i &= f_{\text{drain},i} h_i^n
      + \frac{\Delta t\sin\theta}{w_i\Delta x_i}
        (T_{i+\frac12}^n - T_{i-\frac12}^n)
      + \Delta t R_i.
\end{aligned}
$$

Upper boundary ($i=N$, divide)

$$
\begin{aligned}
a_N &= -\frac{T_{N-\frac12}^n \cos\theta \,\Delta t}
            {\Delta x_{N-\frac12}\,\Delta x_N\,w_N}, \\
c_N &= 0, \\
b_N &= f_{\text{drain},N} - a_N, \\
r_N &= f_{\text{drain},N} h_N^n
      - \frac{\Delta t\sin\theta}{w_N\Delta x_N} T_{N-\frac12}^n
      + \Delta t R_N.
\end{aligned}
$$

### Temporal Discretization

A backward-Euler time step is used for stability. Nonlinear terms in transmissivity and porosity are treated by Picard iteration, updating T and $f_{drain} until:

$$
\max_i |h_i^{k+1} - h_i^{k}| < 10^{-4}\ \mathrm{m}
$$

If convergence fails, the time step is halved adaptively. 

$$
\Delta t_{h3d}^{new} = 0.5\,\Delta t_{h3d}^{old}
$$

Sub-steps are accumulated until the total integration time equals the parent ELM time step:

$$
\sum \Delta t_{h3d} = \Delta t_{ELM}
$$

### Subsurface Runoff and Storage Change

$$
\Delta S_{\text{sat},i}
= f_{\text{drain},i}(h_i^{n+1}-h_i^n),
\qquad
R_{\text{sub},i} = -\Delta S_{\text{sat},i},
\qquad
Q_{\text{sub},i} = \frac{R_{\text{sub},i}}{\Delta t}\times1000
$$

Water-Table Depth

$$
z_{wt,i} = z_{\text{bed},i} - h_i
$$

## Subroutines and workflow

### DrainageH3D

Top-level routine for h3D hydrology.
Responsible for preparing column-level variables (layer thickness, water table index, slope, conductivity), 
computing area-weighted parameters, and calling the h3D solver (H3D_DRI).

### H3D_DRI

Performs the iterative time stepping of the hillslope system:

- Initializes the saturated thickness $h_{sat}$ for each h3D column.

- Computes area-weighted inputs (mean slope, width, transmissivity, decay factor).

- Advances $h_{sat}$ over sub-steps by calling LateralResponse.

- Converts changes in saturated storage to drainage flux:

_{sat} = f_{drain}\,(h^{n+1} - h^{n}),
\qquad
Q_{sub} = -\frac{\Delta S_{sat}}{\Delta t}
$$

Outputs updated water-table depth and drainage rates (qflx_drain_h3d).

### LateralResponse

Solves the implicit Dupuit–Boussinesq system for all h3D nodes in a landunit.
Constructs a tridiagonal matrix from the finite-difference discretization:

$$
a_i h_{i-1}^{n+1} + b_i h_i^{n+1} + c_i h_{i+1}^{n+1} = r_i
$$

Steps:
- Computes node-specific yield $f_{drain}(c)$ and transmissivity $wK!H(k)$.

- Applies slope-dependent flux terms and boundary conditions.

- Solves using the Thomas algorithm (Tridiagonal_h3D).

- Iterates until the solution converges (see § 5.2).

## Outputs

After the h3D solve, the model provides:

- qflx_drain_h3d — subsurface (baseflow) drainage [mm s⁻¹]

- qflx_rsub_sat_h3d — saturation-excess runoff [mm s⁻¹]

- zwt_h3d — updated water-table depth [m]

- f_drain — variable specific yield [–]

- ΔS_sat — change in saturated storage [m]

These outputs replace or augment SIMTOP drainage for h3D-active columns and are fed into the land surface,
 biogeochemical, and river-routing components of ELM.
