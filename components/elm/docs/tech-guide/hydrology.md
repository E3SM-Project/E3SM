# Hydrology

The model parameterizes interception, throughfall, canopy drip, snow
accumulation and melt, water transfer between snow layers, infiltration,
evaporation, surface runoff, sub-surface drainage, redistribution within
the soil column, and groundwater discharge and recharge to simulate
changes in canopy water $\Delta W_{can,\,liq}$, canopy snow water
$\Delta W_{can,\,sno}$ surface water $\Delta W_{sfc}$, snow water
$\Delta W_{sno}$, soil water $\Delta w_{liq,\, i}$, and soil ice
$\Delta w_{ice,\, i}$, and water in the unconfined aquifer
$\Delta W_{a}$ (all in $kg\ m^{-2}$ or $mm\ of\ H_2O$)
(`Figure Hydrologic processes`{.interpreted-text role="numref"}).

The total water balance of the system is

$$
\begin{aligned}
\Delta W_{can,\,liq} + \Delta W_{can,\,sno} + \Delta W_{sfc} + \Delta W_{sno} +
\sum_{i=1}^{N_{levsoi}} (\Delta w_{liq,\, i} + \Delta w_{ice,\, i}) + \Delta W_{a} \\
= \left( q_{rain} + q_{sno} - E_{v} - E_{g} - q_{over} - q_{h2osfc} - q_{drai} - q_{rgwl} - q_{snwcp,\, ice} \right) \Delta t
\end{aligned}
$$

where $q_{rain}$ is the liquid part of precipitation, $q_{sno}$ is the
solid part of precipitation, $E_{v}$ is ET from vegetation (Chapter
`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`{.interpreted-text
role="numref"}), $E_{g}$ is ground evaporation (Chapter
`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`{.interpreted-text
role="numref"}), $q_{over}$ is surface runoff (section
`Surface Runoff`{.interpreted-text role="numref"}), $q_{h2osfc}$ is
runoff from surface water storage (section
`Surface Runoff`{.interpreted-text role="numref"}), $q_{drai}$ is
sub-surface drainage (section
`Lateral Sub-surface Runoff`{.interpreted-text role="numref"}),
$q_{rgwl}$ and $q_{snwcp,ice}$ are liquid and solid runoff from glaciers
and lakes, and runoff from other surface types due to snow capping
(section
`Runoff from glaciers and snow-capped surfaces`{.interpreted-text
role="numref"}) (all in $kg\ m^{-2} s^{-1}$), $N_{levsoi}$ is the number of
soil layers (note that hydrology calculations are only done over soil
layers 1 to $N_{levsoi}$; ground levels $N_{levsoi} +1$ to $N_{levgrnd}$
are currently hydrologically inactive;
`(Lawrence et al. 2008) <Lawrenceetal2008>`{.interpreted-text
role="ref"} and $\Delta t$ is the time step (s).

::: {#Figure Hydrologic processes}
![Hydrologic processes represented in CLM.](hydrologic.processes.png)
:::

## Canopy Water

Liquid precipitation is either intercepted by the canopy, falls directly
to the snow/soil surface (throughfall), or drips off the vegetation
(canopy drip). Solid precipitation is treated similarly, with the
addition of unloading of previously intercepted snow. Interception by
vegetation is divided between liquid and solid phases $q_{intr,\,liq}$
and $q_{intr,\,ice}$ ($kg\ m^{-2} s^{-1}$)

$$q_{intr,\,liq} = f_{pi,\,liq} \ q_{rain}$$

$$q_{intr,\,ice} = f_{pi,\,ice} \ q_{sno}$$

where $f_{pi,\,liq}$ and $f_{pi,\,ice}$ are the fractions of intercepted
precipitation of rain and snow, respectively

$$f_{pi,\,liq} = \alpha_{liq} \ tanh \left(L+S\right)$$

$$f_{pi,\,ice} =\alpha_{sno} \ \left\(1-\exp \left[-0.5\left(L+S\right)\right]\right\) $$

and $L$ and $S$ are the exposed leaf and stem area index, respectively
(section `Phenology and vegetation burial by snow`{.interpreted-text
role="numref"}), and the $\alpha$\'s scale the fractional area of a leaf
that collects water
(`Lawrence et al. 2007 <Lawrenceetal2007>`{.interpreted-text
role="ref"}). Default values of $\alpha_{liq}$ and $\alpha_{sno}$ are
set to 1. Throughfall ($kg\ m^{-2} s^{-1}$) is also divided into liquid and
solid phases, reaching the ground (soil or snow surface) as

$$q_{thru,\, liq} = q_{rain} \left(1 - f_{pi,\,liq}\right)$$

$$q_{thru,\, ice} = q_{sno} \left(1 - f_{pi,\,ice}\right)$$

Similarly, the liquid and solid canopy drip fluxes are

$$q_{drip,\, liq} =\frac{W_{can,\,liq}^{intr} -W_{can,\,liq}^{max } }{\Delta t} \ge 0$$

$$q_{drip,\, ice} =\frac{W_{can,\,sno}^{intr} -W_{can,\,sno}^{max } }{\Delta t} \ge 0$$

where

$$W_{can,liq}^{intr} =W_{can,liq}^{n} +q_{intr,\, liq} \Delta t\ge 0$$

and

$$W_{can,sno}^{intr} =W_{can,sno}^{n} +q_{intr,\, ice} \Delta t\ge 0$$

are the the canopy liquid water and snow water equivalent after
accounting for interception, $W_{can,\,liq}^{n}$ and $W_{can,\,sno}^{n}$
are the canopy liquid and snow water from the previous time step, and
$W_{can,\,liq}^{max }$ and $W_{can,\,snow}^{max }$ ($kg\ m^{-2}$ or $mm\ of\
H_2O) are the maximum amounts of liquid water and snow the canopy can
hold. They are defined by

$$W_{can,\,liq}^{max } =p_{liq}\left(L+S\right)$$

$$W_{can,\,sno}^{max } =p_{sno}\left(L+S\right).$$

The maximum storage of liquid water is $p_{liq}=0.1$ $kg\ m^{-2}$
(`Dickinson et al. 1993 <Dickinsonetal1993>`{.interpreted-text
role="ref"}), and that of snow is $p_{sno}=6$ $kg\ m^{-2}$, consistent with
reported field measurements
(`Pomeroy et al. 1998 <Pomeroyetal1998>`{.interpreted-text role="ref"}).

Canopy snow unloading from wind speed $u$ and above-freezing
temperatures are modeled from linear fluxes and e-folding times similar
to `Roesch et al. (2001) <Roeschetal2001>`{.interpreted-text role="ref"}

$$q_{unl,\, wind} =\frac{u W_{can,sno}}{1.56\times 10^5 \text{ m}}$$

$$q_{unl,\, temp} =\frac{W_{can,sno}(T-270 \textrm{ K})}{1.87\times 10^5 \text{ K s}} > 0$$

$$q_{unl,\, tot} =\min \left( q_{unl,\, wind} +q_{unl,\, temp} ,W_{can,\, sno} \right)$$

The canopy liquid water and snow water equivalent are updated as

$$W_{can,\, liq}^{n+1} =W_{can,liq}^{n} + q_{intr,\, liq} - q_{drip,\, liq} \Delta t - E_{v}^{liq} \Delta t \ge 0$$

and

$$W_{can,\, sno}^{n+1} =W_{can,sno}^{n} + q_{intr,\, ice} - \left(q_{drip,\, ice}+q_{unl,\, tot} \right)\Delta t
                      - E_{v}^{ice} \Delta t \ge 0$$

where $E_{v}^{liq}$ and $E_{v}^{ice}$ are partitioned from the stem and
leaf surface evaporation $E_{v}^{w}$ (Chapter
`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`{.interpreted-text
role="numref"}) based on the vegetation temperature $T_{v}$ (K) (Chapter
`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`{.interpreted-text
role="numref"}) and its relation to the freezing temperature of water
$T_{f}$ (K) (`Table Physical Constants`{.interpreted-text
role="numref"})

$$\begin{aligned}
E_{v}^{liq} =
\left\(\begin{array}{lr}
E_{v}^{w} &  T_v > T_{f} \\
0         &  T_v \le T_f
\end{array}\right\)
\end{aligned}$$

$$\begin{aligned}
E_{v}^{ice} =
\left\(\begin{array}{lr}
0         & T_v > T_f \\
E_{v}^{w} & T_v \le T_f
\end{array}\right\).
\end{aligned}$$

The total rate of liquid and solid precipitation reaching the ground is
then

$$q_{grnd,liq} =q_{thru,\, liq} +q_{drip,\, liq}$$

$$q_{grnd,ice} =q_{thru,\, ice} +q_{drip,\, ice} +q_{unl,\, tot} .$$

Solid precipitation reaching the soil or snow surface,
$q_{grnd,\, ice} \Delta t$, is added immediately to the snow pack
(Chapter `rst_Snow Hydrology`{.interpreted-text role="numref"}). The
liquid part, $q_{grnd,\, liq} \Delta t$ is added after surface fluxes
(Chapter
`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`{.interpreted-text
role="numref"}) and snow/soil temperatures (Chapter
`rst_Soil and Snow Temperatures`{.interpreted-text role="numref"}) have
been determined.

The wetted fraction of the canopy (stems plus leaves), which is required
for surface flux (Chapter
`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`{.interpreted-text
role="numref"}) calculations, is
(`Dickinson et al.1993 <Dickinsonetal1993>`{.interpreted-text
role="ref"})

$$
f_{\text{wet}} =
\begin{cases}
\left[\dfrac{W_{\text{can}}}{\rho_{\text{liq}}(L + S)} \right]^{2/3}, & \text{if } L + S > 0 \text{ and } \left[\dfrac{W_{\text{can}}}{\rho_{\text{liq}}(L + S)} \right]^{2/3} \leq 1 \\
1, & \text{if } L + S > 0 \text{ and } \left[\dfrac{W_{\text{can}}}{\rho_{\text{liq}}(L + S)} \right]^{2/3} > 1 \\
0, & \text{if } L + S = 0
\end{cases}
$$

while the fraction of the canopy that is dry and transpiring is

$$
f_{\text{dry}} =
\begin{cases}
\dfrac{(1 - f_{\text{wet}}) L}{L + S}, & \text{if } L + S > 0 \\
0, & \text{if } L + S = 0
\end{cases}
$$

Similarly, the snow-covered fraction of the canopy is used for surface
alebdo when intercepted snow is present (Chapter
`rst_Surface Albedos`{.interpreted-text role="numref"})

$$
f_{\text{can,\,sno}} =
\begin{cases}
\min\left( \left[\dfrac{W_{\text{can,\,sno}}}{\rho_{\text{sno}}(L + S)} \right]^{3/20},\ 1 \right), & \text{if } L + S > 0 \\
0, & \text{if } L + S = 0
\end{cases}
$$

## Surface Runoff, Surface Water Storage, and Infiltration

The moisture input at the grid cell surface,$q_{liq,\, 0}$, is the sum
of liquid precipitation reaching the ground and melt water from snow ($kg\
m^{-2} s^{-1}$). The moisture flux is then partitioned between surface
runoff, surface water storage, and infiltration into the soil.

### Surface Runoff

The simple TOPMODEL-based
(`Beven and Kirkby 1979 <BevenKirkby1979>`{.interpreted-text
role="ref"}) runoff model (SIMTOP) described by
`Niu et al. (2005) <Niuetal2005>`{.interpreted-text role="ref"} is
implemented to parameterize runoff. A key concept underlying this
approach is that of fractional saturated area $f_{sat}$, which is
determined by the topographic characteristics and soil moisture state of
a grid cell. The saturated portion of a grid cell contributes to surface
runoff, $q_{over}$, by the saturation excess mechanism (Dunne runoff)

$$q_{over} =f_{sat} \ q_{liq,\, 0}$$

The fractional saturated area is a function of soil moisture

$$f_{sat} =f_{\max } \ \exp \left(-0.5f_{over} z_{\nabla } \right)$$

where $f_{\max }$ is the potential or maximum value of $f_{sat}$,
$f_{over}$ is a decay factor ($m^{-1}$), and $z_{\nabla}$ is the water
table depth (m) (section `Lateral Sub-surface Runoff`{.interpreted-text
role="numref"}). The maximum saturated fraction, $f_{\max }$, is defined
as the value of the discrete cumulative distribution function (CDF) of
the topographic index when the grid cell mean water table depth is zero.
Thus, $f_{\max }$ is the percent of pixels in a grid cell whose
topographic index is larger than or equal to the grid cell mean
topographic index. It should be calculated explicitly from the CDF at
each grid cell at the resolution that the model is run. However, because
this is a computationally intensive task for global applications,
$f_{\max }$ is calculated once at 0.125° resolution using the 1-km
compound topographic indices (CTIs) based on the HYDRO1K dataset
(`Verdin and Greenlee 1996 <VerdinGreenlee1996>`{.interpreted-text
role="ref"}) from USGS following the algorithm in
`Niu et al. (2005) <Niuetal2005>`{.interpreted-text role="ref"} and then
area-averaged to the desired model resolution (section
`Surface Data`{.interpreted-text role="numref"}). Pixels with CTIs
exceeding the 95 percentile threshold in each 0.125° grid cell are
excluded from the calculation to eliminate biased estimation of
statistics due to large CTI values at pixels on stream networks. For
grid cells over regions without CTIs such as Australia, the global mean
$f_{\max }$ is used to fill the gaps. See
`Li et al. (2013b) <Lietal2013b>`{.interpreted-text role="ref"} for
additional details. The decay factor $f_{over}$ for global simulations
was determined through sensitivity analysis and comparison with observed
runoff to be 0.5 $m^{-1}$.

### Surface Water Storage

A surface water store has been added to the model to represent wetlands
and small, sub-grid scale water bodies. As a result, the wetland land
unit has been removed as of CLM4.5. The state variables for surface
water are the mass of water $W_{sfc}$ ($kg\ m^{-2}$) and temperature
$T_{h2osfc}$ (Chapter `rst_Soil and Snow Temperatures`{.interpreted-text
role="numref"}). Surface water storage and outflow are functions of fine
spatial scale elevation variations called microtopography. The
microtopography is assumed to be distributed normally around the grid
cell mean elevation. Given the standard deviation of the
microtopographic distribution, $\sigma_{micro}$ (m), the fractional
area of the grid cell that is inundated can be calculated. Surface water
storage, $Wsfc$, is related to the height (relative to the grid cell
mean elevation) of the surface water, $d$, by

$$W_{sfc} =\frac{d}{2} \left(1+erf\left(\frac{d}{\sigma_{micro} \sqrt{2} } \right)\right)+\frac{\sigma_{micro} }{\sqrt{2\pi } } e^{\frac{-d^{2} }{2\sigma_{micro} ^{2} } }$$

where $erf$ is the error function. For a given value of $W_{sfc}$,
`7.66`{.interpreted-text role="eq"} can be solved for $d$ using the
Newton-Raphson method. Once $d$ is known, one can determine the fraction
of the area that is inundated as

$$f_{h2osfc} =\frac{1}{2} \left(1+erf\left(\frac{d}{\sigma_{micro} \sqrt{2} } \right)\right)$$

No global datasets exist for microtopography, so the default
parameterization is a simple function of slope

$$\sigma_{micro} =\left(\beta +\beta_{0} \right)^{\eta }$$

where $\beta$ is the topographic slope,
$\beta_{0} =\left(\sigma_{\max } \right)^{\frac{1}{\eta } }$ determines
the maximum value of $\sigma_{micro}$, and $\eta$ is an adjustable
parameter. Default values in the model are $\sigma_{\max } =0.4$ and
$\eta =-3$.

If the spatial scale of the microtopography is small relative to that of
the grid cell, one can assume that the inundated areas are distributed
randomly within the grid cell. With this assumption, a result from
percolation theory can be used to quantify the fraction of the inundated
portion of the grid cell that is interconnected

$$\begin{aligned}
\begin{array}{lr} f_{connected} =\left(f_{h2osfc} -f_{c} \right)^{\mu } & \qquad f_{h2osfc} >f_{c}  \\ f_{connected} =0 &\qquad  f_{h2osfc} \le f_{c}  \end{array}
\end{aligned}$$

where $f_{c}$ is a threshold below which no single connected inundated
area spans the grid cell and $\mu$ is a scaling exponent. Default values
of $f_{c}$ and $\mu$ are 0.4 and 0.14, respectively. When the inundated
fraction of the grid cell surpasses $f_{c}$, the surface water store
acts as a linear reservoir

$$q_{out,h2osfc}=k_{h2osfc} \ f_{connected} \ (Wsfc-Wc)\frac{1}{\Delta t}$$

where $q_{out,h2osfc}$ is the surface water runoff, $k_{h2osfc}$ is a
constant, $Wc$ is the amount of surface water present when
$f_{h2osfc} =f_{c}$, and $\Delta t$ is the model time step. The linear
storage coefficent $k_{h2osfc} = \sin \left(\beta \right)$ is a function
of grid cell mean topographic slope where $\beta$ is the slope in
radians.

### Infiltration

The surface moisture flux remaining after surface runoff has been
removed,

$$q_{in,surface} = (1-f_{sat}) \ q_{liq,\, 0}$$

is divided into inputs to surface water ($q_{in,\, h2osfc}$ ) and the
soil $q_{in,soil}$. If $q_{in,soil}$ exceeds the maximum soil
infiltration capacity ($kg\ m^{-2} s^{-1}$),

$$q_{infl,\, \max } =(1-f_{sat}) \ \Theta_{ice} k_{sat}$$

where $\Theta_{ice}$ is an ice impedance factor (section
`Hydraulic Properties`{.interpreted-text role="numref"}), infiltration
excess (Hortonian) runoff is generated

$$q_{infl,\, excess} =\max \left(q_{in,soil} -\left(1-f_{h2osfc} \right)q_{\inf l,\max } ,0\right)$$

and transferred from $q_{in,soil}$ to $q_{in,h2osfc}$. After evaporative
losses have been removed, these moisture fluxes are

$$q_{in,\, h2osfc} = f_{h2osfc} q_{in,surface} + q_{infl,excess} - q_{evap,h2osfc}$$

and

$$q_{in,soil} = (1-f_{h2osfc} ) \ q_{in,surface} - q_{\inf l,excess} - (1 - f_{sno} - f_{h2osfc} ) \ q_{evap,soil}.$$

The balance of surface water is then calculated as

$$\Delta W_{sfc} =\left(q_{in,h2osfc} - q_{out,h2osfc} - q_{drain,h2osfc} \right) \ \Delta t.$$

Bottom drainage from the surface water store

$$q_{drain,h2osfc} = \min \left(f_{h2osfc} q_{\inf l,\max } ,\frac{W_{sfc} }{\Delta t} \right)$$

is then added to $q_{in,soil}$ giving the total infiltration into the
surface soil layer

$$q_{infl} = q_{in,soil} + q_{drain,h2osfc}$$

Infiltration $q_{infl}$ and explicit surface runoff $q_{over}$ are not
allowed for glaciers.

## Soil Water

Soil water is predicted from a multi-layer model, in which the vertical
soil moisture transport is governed by infiltration, surface and
sub-surface runoff, gradient diffusion, gravity, and canopy
transpiration through root extraction
(`Figure Hydrologic processes`{.interpreted-text role="numref"}).

For one-dimensional vertical water flow in soils, the conservation of
mass is stated as

$$\frac{\partial \theta }{\partial t} =-\frac{\partial q}{\partial z} - e$$

where $\theta$ is the volumetric soil water content (mm^3^ of water /
mm^-3^ of soil), $t$ is time (s), $z$ is height above some datum in the
soil column (mm) (positive upwards), $q$ is soil water flux ($kg\ m^{-2}
s^{-1}$ or $mm\ s^{-1}$) (positive upwards), and $e$ is a soil moisture sink
term ($mm\ of\ water mm^{-1}\ of\ soil\ s^{-1}$) (ET loss). This equation is
solved numerically by dividing the soil column into multiple layers in
the vertical and integrating downward over each layer with an upper
boundary condition of the infiltration flux into the top soil layer
$q_{infl}$ and a zero-flux lower boundary condition at the bottom of the
soil column (sub-surface runoff is removed later in the timestep,
section `Lateral Sub-surface Runoff`{.interpreted-text role="numref"}).

The soil water flux $q$ in equation `7.79`{.interpreted-text role="eq"}
can be described by Darcy\'s law
`(Dingman 2002) <Dingman2002>`{.interpreted-text role="ref"}

$$q = -k \frac{\partial \psi_{h} }{\partial z}$$

where $k$ is the hydraulic conductivity ($mm\ s^{-1}$), and $\psi_{h}$ is
the hydraulic potential (mm). The hydraulic potential is

$$\psi_{h} =\psi_{m} +\psi_{z}$$

where $\psi_{m}$ is the soil matric potential (mm) (which is related to
the adsorptive and capillary forces within the soil matrix), and
$\psi_{z}$ is the gravitational potential (mm) (the vertical distance
from an arbitrary reference elevation to a point in the soil). If the
reference elevation is the soil surface, then $\psi_{z} =z$. Letting
$\psi =\psi_{m}$, Darcy\'s law becomes

$$q = -k \left[\frac{\partial \left(\psi +z\right)}{\partial z} \right].$$

Equation `7.82`{.interpreted-text role="eq"} can be further manipulated
to yield

$$q = -k \left[\frac{\partial \left(\psi +z\right)}{\partial z} \right]
= -k \left(\frac{\partial \psi }{\partial z} + 1 \right) \ .$$

Substitution of this equation into equation `7.79`{.interpreted-text
role="eq"}, with $e = 0$, yields the Richards equation
`(Dingman 2002) <Dingman2002>`{.interpreted-text role="ref"}

$$\frac{\partial \theta }{\partial t} =
\frac{\partial }{\partial z} \left[k\left(\frac{\partial \psi }{\partial z} + 1
\right)\right].$$

In practice (Section `Numerical Solution Hydrology`{.interpreted-text
role="numref"}), changes in soil water content are predicted from
`7.79`{.interpreted-text role="eq"} using finite-difference
approximations for `7.84`{.interpreted-text role="eq"}.

### Hydraulic Properties

The hydraulic conductivity $k_{i}$ ($mm\ s^{-1}$) and the soil matric
potential $\psi_{i}$ (mm) for layer $i$ vary with volumetric soil water
$\theta_{i}$ and soil texture. As with the soil thermal properties
(section `Soil And Snow Thermal Properties`{.interpreted-text
role="numref"}) the hydraulic properties of the soil are assumed to be a
weighted combination of the mineral properties, which are determined
according to sand and clay contents based on work by
`Clapp and Hornberger (1978) <ClappHornberger1978>`{.interpreted-text
role="ref"} and `Cosby et al. (1984) <Cosbyetal1984>`{.interpreted-text
role="ref"}, and organic properties of the soil
(`Lawrence and Slater 2008 <LawrenceSlater2008>`{.interpreted-text
role="ref"}).

The hydraulic conductivity is defined at the depth of the interface of
two adjacent layers $z_{h,\, i}$
(`Figure Water flux schematic`{.interpreted-text role="numref"}) and is
a function of the saturated hydraulic conductivity
$k_{sat} \left[z_{h,\, i} \right]$, the liquid volumetric soil moisture
of the two layers $\theta_{i}$ and $\theta_{i+1}$ and an ice impedance
factor $\Theta_{ice}$

$$
k\left[z_{h,\, i} \right] =
\begin{cases}
\Theta_{\text{ice}} \, k_{\text{sat}}\left[z_{h,\, i} \right] \left[\dfrac{0.5\left(\theta_{i} + \theta_{i+1} \right)}{0.5\left(\theta_{\text{sat},\, i} + \theta_{\text{sat},\, i+1} \right)} \right]^{2B_{i} + 3}, & \text{for } 1 \le i \le N_{\text{levsoi}} - 1 \\
\Theta_{\text{ice}} \, k_{\text{sat}}\left[z_{h,\, i} \right] \left(\dfrac{\theta_{i}}{\theta_{\text{sat},\, i}} \right)^{2B_{i} + 3}, & \text{for } i = N_{\text{levsoi}}
\end{cases}
$$

The ice impedance factor is a function of ice content, and is meant to
quantify the increased tortuosity of the water flow when part of the
pore space is filled with ice.
`Swenson et al. (2012) <Swensonetal2012>`{.interpreted-text role="ref"}
used a power law form

$$\Theta_{ice} = 10^{-\Omega F_{ice} }$$

where $\Omega = 6$ and $F_{ice} = \frac{\theta_{ice} }{\theta_{sat} }$ is
the ice-filled fraction of the pore space.

Because the hydraulic properties of mineral and organic soil may differ
significantly, the bulk hydraulic properties of each soil layer are
computed as weighted averages of the properties of the mineral and
organic components. The water content at saturation (i.e. porosity) is

$$\theta_{sat,i} =(1-f_{om,i} )\theta_{sat,\min ,i} +f_{om,i} \theta_{sat,om}$$

where $f_{om,i}$ is the soil organic matter fraction, $\theta_{sat,om}$
is the porosity of organic matter, and the porosity of the mineral soil
$\theta_{sat,\min,i}$ is

$$\theta_{sat, min, i} = 0.489 - 0.00126(\text{% sand})_{i}$$

The exponent $B_{i}$ is

$$B_{i} =(1-f_{om,i} )B_{\min ,i} +f_{om,i} B_{om}$$

where $B_{om}$ is for organic matter and

$$B_{\min ,i} =2.91+0.159(\text{% clay})_{i} $$

The soil matric potential (mm) is defined at the node depth $z_{i}$ of
each layer $i$ (`Figure Water flux schematic`{.interpreted-text
role="numref"})

$$\psi_{i} =\psi_{sat,\, i} \left(\frac{\theta_{\, i} }{\theta_{sat,\, i} } \right)^{-B_{i} } \ge -1\times 10^{8} \qquad 0.01\le \frac{\theta_{i} }{\theta_{sat,\, i}} \le 1$$

where the saturated soil matric potential (mm) is

$$\psi_{sat,i} =(1-f_{om,i} )\psi_{sat,\min ,i} +f_{om,i} \psi_{sat,om}$$

where $\psi_{sat,om}$ is the saturated organic matter matric potential
and the saturated mineral soil matric potential $\psi_{sat,\min,i}$ is

$$\psi_{sat,\, \min ,\, i} =-10.0\times 10^{1.88-0.0131(\text{% sand})_{i}} $$

The saturated hydraulic conductivity, $k_{sat} \left[z_{h,\, i} \right]$
($mm\ s^{-1}$), for organic soils ($k_{sat,\, om}$ ) may be two to three
orders of magnitude larger than that of mineral soils
($k_{sat,\, \min }$ ). Bulk soil layer values of $k_{sat}$ calculated as
weighted averages based on $f_{om}$ may therefore be determined
primarily by the organic soil properties even for values of $f_{om}$ as
low as 1 %. To better represent the influence of organic soil material
on the grid cell average saturated hydraulic conductivity, the soil
organic matter fraction is further subdivided into \"connected\" and
\"unconnected\" fractions using a result from percolation theory
(`Stauffer and Aharony 1994 <StaufferAharony1994>`{.interpreted-text
role="ref"},
`Berkowitz and Balberg 1992 <BerkowitzBalberg1992>`{.interpreted-text
role="ref"}). Assuming that the organic and mineral fractions are
randomly distributed throughout a soil layer, percolation theory
predicts that above a threshold value $f_{om} = f_{threshold}$,
connected flow pathways consisting of organic material only exist and
span the soil space. Flow through these pathways interacts only with
organic material, and thus can be described by $k_{sat,\, om}$. This
fraction of the grid cell is given by

$$
f_{\text{perc}} =
\begin{cases}
N_{\text{perc}} \left(f_{\text{om}} - f_{\text{threshold}} \right)^{\beta_{\text{perc}}} f_{\text{om}}, & \text{if } f_{\text{om}} \ge f_{\text{threshold}} \\
0, & \text{if } f_{\text{om}} < f_{\text{threshold}}
\end{cases}
$$

where $\beta ^{perc} =0.139$, $f_{threshold} =0.5$, and
$N_{perc} =\left(1-f_{threshold} \right)^{-\beta_{perc} }$. In the
unconnected portion of the grid cell,
$f_{uncon} =\; \left(1-f_{perc} {\rm \; }\right)$, the saturated
hydraulic conductivity is assumed to correspond to flow pathways that
pass through the mineral and organic components in series

$$k_{sat,\, uncon} =f_{uncon} \left(\frac{\left(1-f_{om} \right)}{k_{sat,\, \min } } +\frac{\left(f_{om} -f_{perc} \right)}{k_{sat,\, om} } \right)^{-1} .$$

where saturated hydraulic conductivity for mineral soil depends on soil
texture (`Cosby et al. 1984 <Cosbyetal1984>`{.interpreted-text
role="ref"}) as

$$k_{sat,\, \min } \left[z_{h,\, i} \right]=0.0070556\times 10^{-0.884+0.0153\left(\text{% sand}\right)_{i} } .$$

The bulk soil layer saturated hydraulic conductivity is then computed as

$$k_{sat} \left[z_{h,\, i} \right]=f_{uncon,\, i} k_{sat,\, uncon} \left[z_{h,\, i} \right]+(1-f_{uncon,\, i} )k_{sat,\, om} \left[z_{h,\, i} \right].$$

The soil organic matter properties implicitly account for the standard
observed profile of organic matter properties as

$$\theta_{sat,om} = max(0.93 - 0.1\times z_{i} / zsapric, 0.83).$$

$$B_{om} = min(2.7 + 9.3\times z_{i} / zsapric, 12.0).$$

$$\psi_{sat,om} = min(10.3 - 0.2\times z_{i} / zsapric, 10.1).$$

$$k_{sat,om} = max(0.28 - 0.2799\times z_{i} / zsapric, k_{sat,\, \min } \left[z_{h,\, i} \right]).$$

where $zsapric =0.5$ m is the depth that organic matter takes on the
characteristics of sapric peat.

### Numerical Solution

With reference to `Figure Water flux schematic`{.interpreted-text
role="numref"}, the equation for conservation of mass (equation
`7.79`{.interpreted-text role="eq"}) can be integrated over each layer
as

$$\int_{-z_{h,\, i} }^{-z_{h,\, i-1} }\frac{\partial \theta }{\partial t} \,  dz=-\int_{-z_{h,\, i} }^{-z_{h,\, i-1} }\frac{\partial q}{\partial z}  \, dz-\int_{-z_{h,\, i} }^{-z_{h,\, i-1} } e\, dz .$$

Note that the integration limits are negative since $z$ is defined as
positive upward from the soil surface. This equation can be written as

$$\Delta z_{i} \frac{\partial \theta_{liq,\, i} }{\partial t} =-q_{i-1} +q_{i} -e_{i}$$

where $q_{i}$ is the flux of water across interface $z_{h,\, i}$,
$q_{i-1}$ is the flux of water across interface $z_{h,\, i-1}$, and
$e_{i}$ is a layer-averaged soil moisture sink term (ET loss) defined as
positive for flow out of the layer ($mm\ s^{-1}$). Taking the finite
difference with time and evaluating the fluxes implicitly at time $n+1$
yields

$$\frac{\Delta z_{i} \Delta \theta_{liq,\, i} }{\Delta t} =-q_{i-1}^{n+1} +q_{i}^{n+1} -e_{i}$$

where
$\Delta \theta_{liq,\, i} =\theta_{liq,\, i}^{n+1} -\theta_{liq,\, i}^{n}$
is the change in volumetric soil liquid water of layer $i$ in time
$\Delta t$and $\Delta z_{i}$ is the thickness of layer $i$ (mm).

The water removed by transpiration in each layer $e_{i}$ is a function
of the total transpiration $E_{v}^{t}$ (Chapter
`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`{.interpreted-text
role="numref"}) and the effective root fraction $r_{e,\, i}$

$$e_{i} =r_{e,\, i} E_{v}^{t} .$$

::: {#Figure Water flux schematic}
![Schematic diagram of numerical scheme used to solve for soil water
fluxes.](image2.png)
:::

Shown are three soil layers, $i-1$, $i$, and $i+1$. The soil matric
potential $\psi$ and volumetric soil water $\theta_{liq}$ are defined at
the layer node depth $z$. The hydraulic conductivity
$k\left[z_{h} \right]$ is defined at the interface of two layers
$z_{h}$. The layer thickness is $\Delta z$. The soil water fluxes
$q_{i-1}$ and $q_{i}$ are defined as positive upwards. The soil moisture
sink term $e$ (ET loss) is defined as positive for flow out of the
layer.

Note that because more than one plant functional type (PFT) may share a
soil column, the transpiration $E_{v}^{t}$ is a weighted sum of
transpiration from all PFTs whose weighting depends on PFT area as

$$
E_{v}^{t} = \sum_{j=1}^{n_{\text{pft}}} \left( E_{v,j}^{t} \cdot \text{wt}_j \right)
$$

where $npft$ is the number of PFTs sharing a soil column,
$\left(E_{v}^{t} \right)\_{j}$ is the transpiration from the $j^{th}$ PFT
on the column, and $\left(wt\right)\_{j}$ is the relative area of the
$j^{th}$ PFT with respect to the column. The effective root fraction
$r_{e,\, i}$ is also a column-level quantity that is a weighted sum over
all PFTs. The weighting depends on the per unit area transpiration of
each PFT and its relative area as

$$r_{e,i} = \frac{\sum_{j=1}^{npft} (r_{e,i})\_j (E_v^t)\_j (wt)\_j}{\sum_{j=1}^{npft} (E_v^t)\_j (wt)\_j}$$

where $\left(r_{e,i} \right)_{j}$ is the effective root fraction for
the $j^{th}$ PFT

$$
\begin{aligned}
\begin{array}{lr}
\left(r_{e,i} \right)\_{j} = \frac{\left(r_{i} \right)\_{j} \left(w_{i} \right)\_{j}}{\left(\beta_{t} \right)\_{j}} & \qquad \left(\beta_{t} \right)\_{j} > 0 \\
\left(r_{e,i} \right)\_{j} = 0 & \qquad \left(\beta_{t} \right)_{j} = 0
\end{array}
\end{aligned}
$$

and $\left(r_{i} \right)_{j}$ is the fraction of roots in layer $i$
(Chapter `rst_Stomatal Resistance and Photosynthesis`{.interpreted-text
role="numref"}), $\left(w_{i} \right)_{j}$ is a soil dryness or plant
wilting factor for layer $i$ (Chapter
`rst_Stomatal Resistance and Photosynthesis`{.interpreted-text
role="numref"}), and $\left(\beta_{t} \right)_{j}$ is a wetness factor
for the total soil column for the $j^{th}$ PFT (Chapter
`rst_Stomatal Resistance and Photosynthesis`{.interpreted-text
role="numref"}).

The soil water fluxes in `7.103`{.interpreted-text role="eq"},, which
are a function of $\theta_{liq,\, i}$ and $\theta_{liq,\, i+1}$ because
of their dependence on hydraulic conductivity and soil matric potential,
can be linearized about $\theta$ using a Taylor series expansion as

$$q_{i}^{n+1} =q_{i}^{n} +\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } \Delta \theta_{liq,\, i} +\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} } \Delta \theta_{liq,\, i+1}$$

$$q_{i-1}^{n+1} =q_{i-1}^{n} +\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} } \Delta \theta_{liq,\, i-1} +\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } \Delta \theta_{liq,\, i} .$$

Substitution of these expressions for $q_{i}^{n+1}$ and $q_{i-1}^{n+1}$
into `7.103`{.interpreted-text role="eq"} results in a general
tridiagonal equation set of the form

$$r_{i} =a_{i} \Delta \theta_{liq,\, i-1} +b_{i} \Delta \theta_{liq,\, i} +c_{i} \Delta \theta_{liq,\, i+1}$$

where

$$a_{i} =-\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} }$$

$$b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}$$

$$c_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} }$$

$$r_{i} =q_{i-1}^{n} -q_{i}^{n} +e_{i} .$$

The tridiagonal equation set is solved over $i=1,\ldots,N_{levsoi}$.

The finite-difference forms of the fluxes and partial derivatives in
equations `7.111`{.interpreted-text role="eq"} -
`7.114`{.interpreted-text role="eq"} can be obtained from equation
`7.82`{.interpreted-text role="eq"} as

$$q_{i-1}^{n} =-k\left[z_{h,\, i-1} \right]\left[\frac{\left(\psi_{i-1} -\psi_{i} \right)+\left(z_{i} - z_{i-1} \right)}{z_{i} -z_{i-1} } \right]$$

$$q_{i}^{n} =-k\left[z_{h,\, i} \right]\left[\frac{\left(\psi_{i} -\psi_{i+1} \right)+\left(z_{i+1} - z_{i} \right)}{z_{i+1} -z_{i} } \right]$$

$$\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} } =-\left[\frac{k\left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } \frac{\partial \psi_{i-1} }{\partial \theta_{liq,\, i-1} } \right]-\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta_{liq,\, i-1} } \left[\frac{\left(\psi_{i-1} -\psi_{i} \right)+\left(z_{i} - z_{i-1} \right)}{z_{i} - z_{i-1} } \right]$$

$$\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } =\left[\frac{k\left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } \frac{\partial \psi_{i} }{\partial \theta_{liq,\, i} } \right]-\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta_{liq,\, i} } \left[\frac{\left(\psi_{i-1} -\psi_{i} \right)+\left(z_{i} - z_{i-1} \right)}{z_{i} - z_{i-1} } \right]$$

$$\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } =-\left[\frac{k\left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \frac{\partial \psi_{i} }{\partial \theta_{liq,\, i} } \right]-\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta_{liq,\, i} } \left[\frac{\left(\psi_{i} -\psi_{i+1} \right)+\left(z_{i+1} - z_{i} \right)}{z_{i+1} - z_{i} } \right]$$

$$\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} } =\left[\frac{k\left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \frac{\partial \psi_{i+1} }{\partial \theta_{liq,\, i+1} } \right]-\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta_{liq,\, i+1} } \left[\frac{\left(\psi_{i} -\psi_{i+1} \right)+\left(z_{i+1} - z_{i} \right)}{z_{i+1} - z_{i} } \right].$$

The derivatives of the soil matric potential at the node depth are
derived from `7.94`{.interpreted-text role="eq"}

$$\frac{\partial \psi_{i-1} }{\partial \theta_{liq,\, \, i-1} } =-B_{i-1} \frac{\psi_{i-1} }{\theta_{\, \, i-1} }$$

$$\frac{\partial \psi_{i} }{\partial \theta_{\, liq,\, i} } =-B_{i} \frac{\psi_{i} }{\theta_{i} }$$

$$\frac{\partial \psi_{i+1} }{\partial \theta_{liq,\, i+1} } =-B_{i+1} \frac{\psi_{i+1} }{\theta_{\, i+1} }$$

with the constraint
$0.01\, \theta_{sat,\, i} \le \theta_{\, i} \le \theta_{sat,\, i}$.

The derivatives of the hydraulic conductivity at the layer interface are
derived from `7.85`{.interpreted-text role="eq"}

$$
\begin{array}{l}
\frac{\partial k\left[z_{h,i-1} \right]}{\partial \theta_{liq,i-1}} 
= \frac{\partial k\left[z_{h,i-1} \right]}{\partial \theta_{liq,i}} 
= \left(2B_{i-1} + 3\right) \, \overline{\Theta}\_{ice} \, k_{sat}\left[z_{h,i-1} \right] 
\left( \frac{\overline{\theta}\_{liq}}{\overline{\theta}\_{sat}} \right)^{2B_{i-1} + 2} 
\left( \frac{0.5}{\overline{\theta}_{sat}} \right)
\end{array}
$$

where $\overline{\Theta}\_{ice} = \Theta(\overline{\theta}\_{ice})$
`7.86`{.interpreted-text role="eq"},
$\overline{\theta}\_{ice} = 0.5\left(\theta_{ice\, i-1} +\theta_{ice\, i} \right)$,
$\overline{\theta}\_{liq} = 0.5\left(\theta_{liq\, i-1} +\theta_{liq\, i} \right)$,
and
$\overline{\theta}\_{sat} = 0.5\left(\theta_{sat,\, i-1} +\theta_{sat,\, i} \right)$

and

$$
\begin{array}{l}
{\frac{\partial k\left[z_{h,i} \right]}{\partial \theta_{liq,i}}
= \frac{\partial k\left[z_{h,i} \right]}{\partial \theta_{liq,i+1}}
= \left(2B_{i} +3\right) \, \overline{\Theta}\_{ice} \ k_{sat} \left[z_{h,i} \right] 
\ \left(\frac{\overline{\theta}\_{liq}}{\overline{\theta}\_{sat}} \right)^{2B_{i} +2} \left(\frac{0.5}{\overline{\theta}_{sat}} \right)} \end{array}.$$

where
$\overline{\theta}\_{liq} = 0.5\left(\theta_{\, i} +\theta_{\, i+1} \right)$,
$\overline{\theta}\_{sat} = 0.5\left(\theta_{sat,\, i} +\theta_{sat,\, i+1} \right)$.

#### Equation set for layer $i=1$

For the top soil layer ($i=1$), the boundary condition is the
infiltration rate (section `Surface Runoff`{.interpreted-text
role="numref"}), $q_{i-1}^{n+1} =-q_{infl}^{n+1}$, and the water balance
equation is

$$\frac{\Delta z_{i} \Delta \theta_{liq,\, i} }{\Delta t} =q_{infl}^{n+1} +q_{i}^{n+1} -e_{i} .$$

After grouping like terms, the coefficients of the tridiagonal set of
equations for $i=1$ are

$$a_{i} =0$$

$$b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}$$

$$c_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} }$$

$$r_{i} =q_{infl}^{n+1} -q_{i}^{n} +e_{i} .$$

#### Equation set for layers $i=2,\ldots ,N_{levsoi} -1$

The coefficients of the tridiagonal set of equations for
$i=2,\ldots,N_{levsoi} -1$ are

$$a_{i} =-\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} }$$

$$b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}$$

$$c_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} }$$

$$r_{i} =q_{i-1}^{n} -q_{i}^{n} +e_{i} .$$

#### Equation set for layer $i=N_{levsoi}$

For the lowest soil layer ($i=N_{levsoi}$ ), a zero-flux bottom boundary
condition is applied ($q_{i}^{n} =0$) and the coefficients of the
tridiagonal set of equations for $i=N_{levsoi}$ are

$$a_{i} =-\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} }$$

$$b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}$$

$$c_{i} =0$$

$$r_{i} =q_{i-1}^{n} +e_{i} .$$

#### Adaptive Time Stepping

The length of the time step is adjusted in order to improve the accuracy
and stability of the numerical solutions. The difference between two
numerical approximations is used to estimate the temporal truncation
error, and then the step size $\Delta t_{sub}$ is adjusted to meet a
user-prescribed error tolerance
`[Kavetski et al., 2002]<Kavetskietal2002>`{.interpreted-text
role="ref"}. The temporal truncation error is estimated by comparing the
flux obtained from the first-order Taylor series expansion
($q_{i-1}^{n+1}$ and $q_{i}^{n+1}$, equations `7.108`{.interpreted-text
role="eq"} and `7.109`{.interpreted-text role="eq"}) against the flux at
the start of the time step ($q_{i-1}^{n}$ and $q_{i}^{n}$). Since the
tridiagonal solution already provides an estimate of
$\Delta \theta_{liq,i}$, it is convenient to compute the error for each
of the $i$ layers from equation `7.103`{.interpreted-text role="eq"} as

$$\epsilon_{i} = \left[ \frac{\Delta \theta_{liq,\, i} \Delta z_{i}}{\Delta t_{sub}} -
\left( q_{i-1}^{n} - q_{i}^{n} + e_{i}\right) \right] \ \frac{\Delta t_{sub}}{2}$$

and the maximum absolute error across all layers as

$$\begin{array}{lr}
\epsilon_{crit} = {\rm max} \left( \left| \epsilon_{i} \right| \right) & \qquad 1 \le i \le nlevsoi
\end{array} \ .$$

The adaptive step size selection is based on specified upper and lower
error tolerances, $\tau_{U}$ and $\tau_{L}$. The solution is accepted if
$\epsilon_{crit} \le \tau_{U}$ and the procedure repeats until the
adaptive sub-stepping spans the full model time step (the sub-steps are
doubled if $\epsilon_{crit} \le \tau_{L}$, i.e., if the solution is very
accurate). Conversely, the solution is rejected if
$\epsilon_{crit} > \tau_{U}$. In this case the length of the sub-steps
is halved and a new solution is obtained. The halving of substeps
continues until either $\epsilon_{crit} \le \tau_{U}$ or the specified
minimum time step length is reached.

Upon solution of the tridiagonal equation set, the liquid water contents
are updated as follows

$$w_{liq,\, i}^{n+1} =w_{liq,\, i}^{n} +\Delta \theta_{liq,\, i} \Delta z_{i} \qquad i=1,\ldots ,N_{levsoi} .$$

The volumetric water content is

$$\theta_{i} =\frac{w_{liq,\, i} }{\Delta z_{i} \rho_{liq} } +\frac{w_{ice,\, i} }{\Delta z_{i} \rho_{ice} } .$$

## Frozen Soils and Perched Water Table

When soils freeze, the power-law form of the ice impedance factor
(section `Hydraulic Properties`{.interpreted-text role="numref"}) can
greatly decrease the hydraulic conductivity of the soil, leading to
nearly impermeable soil layers. When unfrozen soil layers are present
above relatively ice-rich frozen layers, the possibility exists for
perched saturated zones. Lateral drainage from perched saturated regions
is parameterized as a function of the thickness of the saturated zone

$$q_{drai,perch} =k_{drai,\, perch} \left(z_{frost} -z_{\nabla ,perch} \right)$$

where $k_{drai,\, perch}$ depends on topographic slope and soil
hydraulic conductivity,

$$k_{drai,\, perch} =10^{-5} \sin (\beta )\left(\frac{\sum_{i=N_{perch} }^{i=N_{frost} }\Theta_{ice,i} k_{sat} \left[z_{i} \right]\Delta z_{i}  }{\sum_{i=N_{perch} }^{i=N_{frost} }\Delta z_{i}  } \right)$$

where $\Theta_{ice}$ is an ice impedance factor, $\beta$ is the mean
grid cell topographic slope in radians, $z_{frost}$ is the depth to the
frost table, and $z_{\nabla,perch}$ is the depth to the perched
saturated zone. The frost table $z_{frost}$ is defined as the shallowest
frozen layer having an unfrozen layer above it, while the perched water
table $z_{\nabla,perch}$ is defined as the depth at which the volumetric
water content drops below a specified threshold. The default threshold
is set to 0.9. Drainage from the perched saturated zone $q_{drai,perch}$
is removed from layers $N_{perch}$ through $N_{frost}$, which are the
layers containing $z_{\nabla,perch}$ and, $z_{frost}$ respectively.

## Lateral Sub-surface Runoff

Lateral sub-surface runoff occurs when saturated soil moisture
conditions exist within the soil column. Sub-surface runoff is

$$q_{drai} = \Theta_{ice} K_{baseflow} tan \left( \beta \right)
\Delta z_{sat}^{N_{baseflow}} \ ,$$

where $K_{baseflow}$ is a calibration parameter, $\beta$ is the
topographic slope, the exponent $N_{baseflow}$ = 1, and $\Delta z_{sat}$
is the thickness of the saturated portion of the soil column.

The saturated thickness is

$$\Delta z_{sat} = z_{bedrock} - z_{\nabla},$$

where the water table $z_{\nabla}$ is determined by finding the first
soil layer above the bedrock depth (section
`Depth to Bedrock`{.interpreted-text role="numref"}) in which the
volumetric water content drops below a specified threshold. The default
threshold is set to 0.9.

The specific yield, $S_{y}$, which depends on the soil properties and
the water table location, is derived by taking the difference between
two equilibrium soil moisture profiles whose water tables differ by an
infinitesimal amount

$$S_{y} =\theta_{sat} \left(1-\left(1+\frac{z_{\nabla } }{\Psi_{sat} } \right)^{\frac{-1}{B} } \right)$$

where B is the Clapp-Hornberger exponent. Because $S_{y}$ is a function
of the soil properties, it results in water table dynamics that are
consistent with the soil water fluxes described in section
`Soil Water`{.interpreted-text role="numref"}.

After the above calculations, two numerical adjustments are implemented
to keep the liquid water content of each soil layer ($w_{liq,\, i}$ )
within physical constraints of
$w_{liq}^{\min } \le w_{liq,\, i} \le \left(\theta_{sat,\, i} -\theta_{ice,\, i} \right)\Delta z_{i}$
where $w_{liq}^{\min } =0.01$ (mm). First, beginning with the bottom
soil layer $i=N_{levsoi}$, any excess liquid water in each soil layer
($w_{liq,\, i}^{excess} =w_{liq,\, i} -\left(\theta_{sat,\, i} -\theta_{ice,\, i} \right)\Delta z_{i} \ge 0$)
is successively added to the layer above. Any excess liquid water that
remains after saturating the entire soil column is added to drainage
$q_{drai}$. Second, to prevent negative $w_{liq,\, i}$, each layer is
successively brought up to $w_{liq,\, i} =w_{liq}^{\min }$ by taking the
required amount of water from the layer below. If this results in
$w_{liq,\, N_{levsoi} } <w_{liq}^{\min }$, then the layers above are
searched in succession for the required amount of water
($w_{liq}^{\min } -w_{liq,\, N_{levsoi} }$ ) and removed from those
layers subject to the constraint $w_{liq,\, i} \ge w_{liq}^{\min }$. If
sufficient water is not found, then the water is removed from $W_{t}$
and $q_{drai}$.

The soil surface layer liquid water and ice contents are then updated
for dew $q_{sdew}$, frost $q_{frost}$, or sublimation $q_{subl}$
(section
`Update of Ground Sensible and Latent Heat Fluxes`{.interpreted-text
role="numref"}) as

$$w_{liq,\, 1}^{n+1} =w_{liq,\, 1}^{n} +q_{sdew} \Delta t$$

$$w_{ice,\, 1}^{n+1} =w_{ice,\, 1}^{n} +q_{frost} \Delta t$$

$$w_{ice,\, 1}^{n+1} =w_{ice,\, 1}^{n} -q_{subl} \Delta t.$$

Sublimation of ice is limited to the amount of ice available.

## Runoff from glaciers and snow-capped surfaces

All surfaces are constrained to have a snow water equivalent
$W_{sno} \le W_{cap} = 10,000$ $kg\ m^{-2}$. For snow-capped columns, any
addition of mass at the top (precipitation, dew/riping) is balanced by
an equally large mass flux at the bottom of the snow column. This
so-called capping flux is separated into solid $q_{snwcp,ice}$ and
liquid $q_{snwcp,liq}$ runoff terms. The partitioning of these phases is
based on the phase ratio in the bottom snow layer at the time of the
capping, such that phase ratio in this layer is unaltered.

The $q_{snwcp,ice}$ runoff is sent to MOSART (Chapter
`rst_MOSART`{.interpreted-text role="numref"}) where it is routed to the
ocean as an ice stream and, if applicable, the ice is melted there.

For snow-capped surfaces other than glaciers and lakes the
$q_{snwcp,liq}$ runoff is assigned to the glaciers and lakes runoff term
$q_{rgwl}$ (e.g. $q_{rgwl} =q_{snwcp,liq}$ ). For glacier surfaces the
runoff term $q_{rgwl}$ is calculated from the residual of the water
balance

$$q_{rgwl} =q_{grnd,ice} +q_{grnd,liq} -E_{g} -E_{v} -\frac{\left(W_{b}^{n+1} -W_{b}^{n} \right)}{\Delta t} -q_{snwcp,ice}$$

where $W_{b}^{n}$ and $W_{b}^{n+1}$ are the water balances at the
beginning and ending of the time step defined as

$$W_{b} =W_{can} +W_{sno} +\sum_{i=1}^{N}\left(w_{ice,i} +w_{liq,i} \right) .$$

Currently, glaciers are non-vegetated and $E_{v} =W_{can} =0$. The
contribution of lake runoff to $q_{rgwl}$ is described in section
`Precipitation, Evaporation, and Runoff Lake`{.interpreted-text
role="numref"}. The runoff term $q_{rgwl}$ may be negative for glaciers
and lakes, which reduces the total amount of runoff available to the
river routing model (Chapter `rst_MOSART`{.interpreted-text
role="numref"}).
