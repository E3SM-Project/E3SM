# Hybrid-3D hillslope hydrological model

The model represents subsurface flow and groundwater discharge along an
idealized hillslope within each grid cell. It solves the Dupuit–Boussinesq
form of the lateral flow equation for the saturated thickness $h(x,t)$,
using a finite‐difference implicit solver. The formulation includes
variable transmissivity, slope‐driven gradients, variable drainable
porosity, and recharge coupling.

$$
\frac{\partial h}{\partial t} = \frac{1}{f_{\text{drain}}} \frac{\partial}{\partial x} \left[ T(x,h)\left(\frac{\partial h}{\partial }\cos\theta + \sin\theta\right)\right]+ \frac{R}{f_{\text{drain}}}
$$

where $h(x,t)$ is the saturated thickness [m],
$f_{\text{drain}}$ is the drainable porosity [–],
$T(x,h)$ is the transmissivity [m² s⁻¹],
$\theta$ is the hillslope angle [rad],
and $R$ is the recharge rate [m s⁻¹].

## Transmissivity

$$
T(x,h)= \frac{K_{\text{aniso}}\,K_{\text{sat}}(z_{wt})\,h\,w(x)}{1000}
$$

where $K_{\text{aniso}}=100$ is the horizontal/vertical anisotropy factor,
$K_{\text{sat}}(z_{wt})$ is the saturated hydraulic conductivity at the
water-table depth [mm s⁻¹], and $w(x)$ is the hillslope width [m].
Division by 1000 converts from mm s⁻¹ to m s⁻¹.

## Variable Drainable Porosity

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
