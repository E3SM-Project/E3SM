# Hybrid-3D hillslope hydrological model

The model represents subsurface flow and groundwater discharge along an
idealized hillslope within each grid cell. It solves the Dupuit–Boussinesq
form of the lateral flow equation for the saturated thickness $h(x,t)$,
using a finite‐difference implicit solver. The formulation includes
variable transmissivity, slope‐driven gradients, variable drainable
porosity, and recharge coupling.

>$$
\frac{\partial h}{\partial t} = \frac{1}{f_{\text{drain}}} \frac{\partial}{\partial x} \left[ T(x,h)\left(\frac{\partial h}{\partial }\cos\theta + \sin\theta\right)\right]+ \frac{R}{f_{\text{drain}}}
$$

