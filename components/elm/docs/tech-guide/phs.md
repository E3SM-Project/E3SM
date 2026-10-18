# Plant Hydraulics

## Roots {#Roots}

### Vertical Root Distribution {#Vertical Root Distribution}

The root fraction $r_{i}$ in each soil layer depends on the plant
functional type

$$r_{i} =
\begin{array}{lr}
\left(\beta^{z_{h,\, i-1} \cdot 100} - \beta^{z_{h,\, i} \cdot 100} \right) & \qquad {\rm for\; }1 \le i \le N_{levsoi}
\end{array}$$

where $z_{h,\, i}$ (m) is the depth from the soil surface to the
interface between layers $i$ and $i+1$ ($z_{h,\, 0}$, the soil surface)
(section `Vertical Discretization`{.interpreted-text role="numref"}),
the factor of 100 converts from m to cm, and $\beta$ is a
plant-dependent root distribution parameter adopted from
`Jackson et al. (1996)<Jacksonetal1996>`{.interpreted-text role="ref"}
(`Table Plant functional type root distribution parameters`{.interpreted-text
role="numref"}).

::: {#Table Plant functional type root distribution parameters}
  -----------------------------------------------------
  Plant Functional Type              $\beta$
  ---------------------------------- ------------------
  NET Temperate                      0.976

  NET Boreal                         0.943

  NDT Boreal                         0.943

  BET Tropical                       0.993

  BET temperate                      0.966

  BDT tropical                       0.993

  BDT temperate                      0.966

  BDT boreal                         0.943

  BES temperate                      0.964

  BDS temperate                      0.964

  BDS boreal                         0.914

  C~3~ grass arctic                  0.914

  C~3~ grass                         0.943

  C~4~ grass                         0.943

  Crop R                             0.943

  Crop I                             0.943

  Corn R                             0.943

  Corn I                             0.943

  Temp Cereal R                      0.943

  Temp Cereal I                      0.943

  Winter Cereal R                    0.943

  Winter Cereal I                    0.943

  Soybean R                          0.943

  Soybean I                          0.943

  Miscanthus R                       0.943

  Miscanthus I                       0.943

  Switchgrass R                      0.943

  Switchgrass I                      0.943
  -----------------------------------------------------

  : Plant functional type root distribution parameters
:::

### Root Spacing {#Root Spacing}

To determine the conductance along the soil to root pathway (section
`Soil-to-root`{.interpreted-text role="numref"}) an estimate of the
spacing between the roots within a soil layer is required. The distance
between roots $dx_{root,i}$ (m) is calculated by assuming that roots are
distributed uniformly throughout the soil
(`Gardner 1960<Gardner1960>`{.interpreted-text role="ref"})

$$dx_{root,i} = \left(\pi \cdot L_i\right)^{- \frac{1}{2}}$$

where $L_{i}$ is the root length density (m m ^-3^)

$$L_{i} = \frac{B_{root,i}}{\rho_{root} {CA}_{root}} \ ,$$

$B_{root,i}$ is the root biomass density (kg m ^-3^)

$$B_{root,i} = \frac{c\_to\_b \cdot C_{fineroot} \cdot r_{i}}{dz_{i}}$$

where $c\_to\_b = 2$ (kg biomass kg carbon ^-1^) and $C_{fineroot}$ is
the amount of fine root carbon (kg m ^-2^).

$\rho_{root}$ is the root density (kg m ^-3^), and ${CA}_{root}$ is the
fine root cross sectional area (m ^2^)

$$CA_{root} = \pi r_{root}^{2}$$

where $r_{root}$ is the root radius (m).

## Plant Hydraulic Stress {#Plant Hydraulic Stress}

The Plant Hydraulic Stress (PHS) routine explicitly models water
transport through the vegetation according to a simple hydraulic
framework following Darcy\'s Law for porous media flow equations
influenced by `Bonan et al. (2014) <Bonanetal2014>`{.interpreted-text
role="ref"}, `Chuang et al. (2006) <Chuangetal2006>`{.interpreted-text
role="ref"}, `Sperry et al. (1998) <Sperryetal1998>`{.interpreted-text
role="ref"},
`Sperry and Love (2015) <SperryandLove2015>`{.interpreted-text
role="ref"},
`Williams et al (1996) <Williamsetal1996>`{.interpreted-text
role="ref"}.

PHS solves for the vegetation water potential that matches water supply
with transpiration demand. Water supply is modeled according to the
circuit analog in `Figure Plant hydraulic circuit`{.interpreted-text
role="numref"}. Transpiration demand is modeled relative to maximum
transpiration by a transpiration loss function dependent on leaf water
potential.

::: {#Figure Plant hydraulic circuit}
![Circuit diagram of plant hydraulics scheme](circuit.jpg)
:::

### Plant Water Supply {#Plant Water Supply}

The supply equations are used to solve for vegetation water potential
forced by transpiration demand and the set of layer-by-layer soil water
potentials. The water supply is discretized into segments: soil-to-root,
root-to-stem, and stem-to-leaf. There are typically several (1-49)
soil-to-root flows operating in parallel, one per soil layer. There are
two stem-to-leaf flows operating in parallel, corresponding to the
sunlit and shaded \"leaves\".

In general the water fluxes (e.g. soil-to-root, root-to-stem, etc.) are
modeled according to Darcy\'s Law for porous media flow as:

$$q = kA\left( \psi_1 - \psi_2 \right)$$

$q$ is the flux of water (mmH~2~O/s) spanning the segment between
$\psi_1$ and $\psi_2$

$k$ is the hydraulic conductance (s^-1^)

$A$ is the area basis (m^2^/m^2^) relating the conducting area basis to
ground area $\psi_1 - \psi_2$ is the gradient in water potential
(mmH~2~O) across the segment The segments in
`Figure Plant hydraulic circuit`{.interpreted-text role="numref"} have
variable resistance, as the water potentials become lower, hydraulic
conductance decreases. This is captured by multiplying the maximum
segment conductance by a sigmoidal function capturing the percent loss
of conductivity. The function uses two parameters to fit experimental
vulnerability curves: the water potential at 50% loss of conductivity
($p50$) and a shape fitting parameter ($c_k$).

$$k=k_{max}\cdot 2^{-\left(\dfrac{\psi_1}{p50}\right)^{c_k}}$$

$k_{max}$ is the maximum segment conductance (s^-1^)

$p50$ is the water potential at 50% loss of conductivity (mmH~2~O)

$\psi_1$ is the water potential of the lower segment terminus (mmH~2~O)

#### Stem-to-leaf {#Stem-to-leaf}

The area basis and conductance parameterization varies by segment. There
are two stem-to-leaf fluxes in parallel, from stem to sunlit leaf and
from stem to shaded leaf ($q_{1a}$ and $q_{1a}$). The water flux from
stem-to-leaf is the product of the segment conductance, the conducting
area basis, and the water potential gradient from stem to leaf.
Stem-to-leaf conductance is defined as the maximum conductance
multiplied by the percent of maximum conductance, as calculated by the
sigmoidal vulnerability curve. The maximum conductance is a PFT
parameter representing the maximum conductance of water from stem to
leaf per unit leaf area. This parameter can be defined separately for
sunlit and shaded segments and should already include the appropriate
length scaling (in other words this is a conductance, not conductivity).
The water potential gradient is the difference between leaf water
potential and stem water potential. There is no gravity term, assuming a
negligible difference in height across the segment. The area basis is
the leaf area index (either sunlit or shaded).

$$q_{1a}=k_{1a}\cdot\mbox{LAI}_{sun}\cdot\left(\psi_{stem}-\psi_{sunleaf} \right)$$

$$q_{1b}=k_{1b}\cdot\mbox{LAI}_{shade}\cdot\left(\psi_{stem}-\psi_{shadeleaf} \right)$$

$$k_{1a}=k_{1a,max}\cdot 2^{-\left(\dfrac{\psi_{stem}}{p50_1}\right)^{c_k}}$$

$$k_{1b}=k_{1b,max}\cdot 2^{-\left(\dfrac{\psi_{stem}}{p50_1}\right)^{c_k}}$$

Variables:

$q_{1a}$ = flux of water (mmH~2~O/s) from stem to sunlit leaf

$q_{1b}$ = flux of water (mmH~2~O/s) from stem to shaded leaf

$LAI_{sun}$ = sunlit leaf area index (m2/m2)

$LAI_{shade}$ = shaded leaf area index (m2/m2)

$\psi_{stem}$ = stem water potential (mmH~2~O)

$\psi_{sunleaf}$ = sunlit leaf water potential (mmH~2~O)

$\psi_{shadeleaf}$ = shaded leaf water potential (mmH~2~O)

Parameters:

$k_{1a,max}$ = maximum leaf conductance (s^-1^)

$k_{1b,max}$ = maximum leaf conductance (s^-1^)

$p50_{1}$ = water potential at 50% loss of conductance (mmH~2~O)

$c_{k}$ = vulnerability curve shape-fitting parameter (-)

#### Root-to-stem {#Root-to-stem}

There is one root-to-stem flux. This represents a flux from the root
collar to the upper branch reaches. The water flux from root-to-stem is
the product of the segment conductance, the conducting area basis, and
the water potential gradient from root to stem. Root-to-stem conductance
is defined as the maximum conductance multiplied by the percent of
maximum conductance, as calculated by the sigmoidal vulnerability curve
(two parameters). The maximum conductance is defined as the maximum
root-to-stem conductivity per unit stem area (PFT parameter) divided by
the length of the conducting path, which is taken to be the vegetation
height. The area basis is the stem area index. The gradient in water
potential is the difference between the root water potential and the
stem water potential less the difference in gravitational potential.

$$q_2=k_2 \cdot SAI \cdot \left( \psi_{root} - \psi_{stem} - \Delta \psi_z  \right)$$

$$k_2=\dfrac{k_{2,max}}{z_2} \cdot 2^{-\left(\dfrac{\psi_{root}}{p50_2}\right)^{c_k}}$$

Variables:

$q_2$ = flux of water (mmH~2~O/s) from root to stem

$SAI$ = stem area index (m2/m2)

$\Delta\psi_z$ = gravitational potential (mmH~2~O)

$\psi_{root}$ = root water potential (mmH~2~O)

$\psi_{stem}$ = stem water potential (mmH~2~O)

Parameters:

$k_{2,max}$ = maximum stem conductivity (m/s)

$p50_2$ = water potential at 50% loss of conductivity (mmH~2~O)

$z_2$ = vegetation height (m)

#### Soil-to-root {#Soil-to-root}

There are several soil-to-root fluxes operating in parallel (one for
each root-containing soil layer). Each represents a flux from the given
soil layer to the root collar. The water flux from soil-to-root is the
product of the segment conductance, the conducting area basis, and the
water potential gradient from soil to root. The area basis is a proxy
for root area index, defined as the summed leaf and stem area index
multiplied by the root-to-shoot ratio (PFT parameter) multiplied by the
layer root fraction. The root fraction comes from an empirical root
profile (section `Vertical Root Distribution`{.interpreted-text
role="numref"}).

The gradient in water potential is the difference between the soil water
potential and the root water potential less the difference in
gravitational potential. There is only one root water potential to which
all soil layers are connected in parallel. A soil-to-root flux can be
either positive (vegetation water uptake) or negative (water
deposition), depending on the relative values of the root and soil water
potentials. This allows for the occurrence of hydraulic redistribution
where water moves through vegetation tissue from one soil layer to
another.

Soil-to-root conductance is the result of two resistances in series,
first across the soil-root interface and then through the root tissue.
The root tissue conductance is defined as the maximum conductance
multiplied by the percent of maximum conductance, as calculated by the
sigmoidal vulnerability curve. The maximum conductance is defined as the
maximum root-tissue conductivity (PFT parameter) divided by the length
of the conducting path, which is taken to be the soil layer depth plus
lateral root length.

The soil-root interface conductance is defined as the soil conductivity
divided by the conducting length from soil to root. The soil
conductivity varies by soil layer and is calculated based on soil
potential and soil properties, via the Brooks-Corey theory. The
conducting length is determined from the characteristic root spacing
(section `Root Spacing`{.interpreted-text role="numref"}).

$$q_{3,i}=k_{3,i} \cdot \left(\psi_{soil,i}-\psi_{root} + \Delta\psi_{z,i} \right)$$

$$k_{3,i}=\dfrac{k_{r,i} \cdot k_{s,i}}{k_{r,i}+k_{s,i}}$$

$$k_{r,i}=\dfrac{k_{3,max}}{z_{3,i}} \cdot RAI \cdot 2^{-\left(\dfrac{\psi_{soil,i}}{p50_3}\right)^{c_k}}$$

$$RAI=\left(LAI+SAI \right) \cdot r_i \cdot f_{root-leaf}$$

$$k_{s,i} = \dfrac{k_{soil,i}}{dx_{root,i}}$$

Variables:

$q_{3,i}$ = flux of water (mmH~2~O/s) from soil layer $i$ to root

$\Delta\psi_{z,i}$ = change in gravitational potential from soil layer
$i$ to surface (mmH~2~O)

$LAI$ = total leaf area index (m2/m2)

$SAI$ = stem area index (m2/m2)

$\psi_{soil,i}$ = water potential in soil layer $i$ (mmH~2~O)

$\psi_{root}$ = root water potential (mmH~2~O)

$z_{3,i}$ = length of root tissue conducting path = soil layer depth +
root lateral length (m)

$r_i$ = root fraction in soil layer $i$ (-)

$k_{soil,i}$ = Brooks-Corey soil conductivity in soil layer $i$ (m/s)

Parameters:

$f_{root-leaf}$ = root-to-shoot ratio (-)

$p50_3$ = water potential at 50% loss of root tissue conductance
(mmH~2~O)

$ck$ = shape-fitting parameter for vulnerability curve (-)

### Plant Water Demand {#Plant Water Demand}

Plant water demand depends on stomatal conductance, which is described
in section `Stomatal resistance`{.interpreted-text role="numref"}. Here
we describe the influence of PHS and the coupling of vegetation water
demand and supply. PHS models vegetation water demand as transpiration
attenuated by a transpiration loss function based on leaf water
potential. Sunlit leaf transpiration is modeled as the maximum sunlit
leaf transpiration multiplied by the percent of maximum transpiration as
modeled by the sigmoidal loss function. The same follows for shaded leaf
transpiration. Maximum stomatal conductance is calculated from the
Medlyn model `(Medlyn et al. 2011) <Medlynetal2011>`{.interpreted-text
role="ref"} absent water stress and used to calculate the maximum
transpiration (see section
`Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces`{.interpreted-text
role="numref"}). Water stress is calculated as the ratio of attenuated
stomatal conductance to maximum stomatal conductance. Water stress is
calculated with distinct values for sunlit and shaded leaves. Vegetation
water stress is calculated based on leaf water potential and is used to
attenuate photosynthesis (see section `Photosynthesis`{.interpreted-text
role="numref"})

$$E_{sun} = E_{sun,max} \cdot 2^{-\left(\dfrac{\psi_{sunleaf}}{p50_e}\right)^{c_k}}$$

$$E_{shade} = E_{shade,max} \cdot 2^{-\left(\dfrac{\psi_{shadeleaf}}{p50_e}\right)^{c_k}}$$

$$\beta_{t,sun} = \dfrac{g_{s,sun}}{g_{s,sun,\beta_t=1}}$$

$$\beta_{t,shade} = \dfrac{g_{s,shade}}{g_{s,shade,\beta_t=1}}$$

$E_{sun}$ = sunlit leaf transpiration (mm/s)

$E_{shade}$ = shaded leaf transpiration (mm/s)

$E_{sun,max}$ = sunlit leaf transpiration absent water stress (mm/s)

$E_{shade,max}$ = shaded leaf transpiration absent water stress (mm/s)

$\psi_{sunleaf}$ = sunlit leaf water potential (mmH~2~O)

$\psi_{shadeleaf}$ = shaded leaf water potential (mmH~2~O)

$\beta_{t,sun}$ = sunlit transpiration water stress (-)

$\beta_{t,shade}$ = shaded transpiration water stress (-)

$g_{s,sun}$ = stomatal conductance of water corresponding to $E_{sun}$

$g_{s,shade}$ = stomatal conductance of water corresponding to
$E_{shade}$

$g_{s,sun,max}$ = stomatal conductance of water corresponding to
$E_{sun,max}$

$g_{s,shade,max}$ = stomatal conductance of water corresponding to
$E_{shade,max}$

### Vegetation Water Potential {#Vegetation Water Potential}

Both plant water supply and demand are functions of vegetation water
potential. PHS explicitly models root, stem, shaded leaf, and sunlit
leaf water potential at each timestep. PHS iterates to find the
vegetation water potential $\psi$ (vector) that satisfies continuity
between the non-linear vegetation water supply and demand (equations
`11.103`{.interpreted-text role="eq"}, `11.104`{.interpreted-text
role="eq"}, `11.107`{.interpreted-text role="eq"},
`11.109`{.interpreted-text role="eq"}, `11.201`{.interpreted-text
role="eq"}, `11.202`{.interpreted-text role="eq"}).

$$\psi=\left[\psi_{sunleaf},\psi_{shadeleaf},\psi_{stem},\psi_{root}\right]$$

$$\begin{aligned}
\begin{aligned}
E_{sun}&=q_{1a}\\
E_{shade}&=q_{1b}\\
E_{sun}+E_{shade}&=q_{1a}+q_{1b}\\
&=q_2\\
&=\sum_{i=1}^{nlevsoi}{q_{3,i}}
\end{aligned}
\end{aligned}$$

PHS finds the water potentials that match supply and demand. In the
plant water transport equations `11.302`{.interpreted-text role="eq"},
the demand terms (left-hand side) are decreasing functions of absolute
leaf water potential. As absolute leaf water potential becomes larger,
water stress increases, causing a decrease in transpiration demand. The
supply terms (right-hand side) are increasing functions of absolute leaf
water potential. As absolute leaf water potential becomes larger, the
gradients in water potential increase, causing an increase in vegetation
water supply. PHS takes a Newton\'s method approach to iteratively solve
for the vegetation water potentials that satisfy continuity
`11.302`{.interpreted-text role="eq"}.

### Numerical Implementation {#PHS Numerical Implementation}

The four plant water potential nodes are ( $\psi_{root}$,
$\psi_{xylem}$, $\psi_{shadeleaf}$, $\psi_{sunleaf}$). The fluxes
between each pair of nodes are labeled in Figure 1. $E_{sun}$ and
$E_{sha}$ are the transpiration from sunlit and shaded leaves,
respectively. We use the circuit-analog model to calculate the
vegetation water potential ( $\psi$) for the four plant nodes, forced by
soil matric potential and unstressed transpiration. The unstressed
transpiration is acquired by running the photosynthesis model with
$\beta_t=1$. The unstressed transpiration flux is attenuated based on
the leaf-level vegetation water potential. Using the attenuated
transpiration, we solve for $g_{s,stressed}$ and output
$\beta_t=\dfrac{g_{s,stressed}}{g_{s,unstressed}}$.

The continuity of water flow through the system yields four equations

$$\begin{aligned}
\begin{aligned}
E_{sun}&=q_{1a}\\
E_{shade}&=q_{1b}\\
q_{1a}+q_{1b}&=q_2\\
q_2&=\sum_{i=1}^{nlevsoi}{q_{3,i}}
\end{aligned}
\end{aligned}$$

We seek the set of vegetation water potential values,

$$\psi=\left[ \begin {array}{c}
\psi_{sunleaf}\cr\psi_{shadeleaf}\cr\psi_{stem}\cr\psi_{root}
\end {array} \right]$$

that satisfies these equations, as forced by the soil moisture and
atmospheric state. Each flux on the schematic can be represented in
terms of the relevant water potentials. Defining the transpiration
fluxes:

$$\begin{aligned}
\begin{aligned}
E_{sun} &= E_{sun,max} \cdot 2^{-\left(\dfrac{\psi_{sunleaf}}{p50_e}\right)^{c_k}} \\
E_{shade} &= E_{shade,max} \cdot 2^{-\left(\dfrac{\psi_{shadeleaf}}{p50_e}\right)^{c_k}}
\end{aligned}
\end{aligned}$$

Defining the water supply fluxes:

$$\begin{aligned}
\begin{aligned}
q_{1a}&=k_{1a,max}\cdot 2^{-\left(\dfrac{\psi_{stem}}{p50_1}\right)^{c_k}} \cdot\mbox{LAI}_{sun}\cdot\left(\psi_{stem}-\psi_{sunleaf} \right) \\
q_{1b}&=k_{1b,max}\cdot 2^{-\left(\dfrac{\psi_{stem}}{p50_1}\right)^{c_k}}\cdot\mbox{LAI}_{shade}\cdot\left(\psi_{stem}-\psi_{shadeleaf} \right) \\
q_2&=\dfrac{k_{2,max}}{z_2} \cdot 2^{-\left(\dfrac{\psi_{root}}{p50_2}\right)^{c_k}} \cdot SAI \cdot \left( \psi_{root} - \psi_{stem} - \Delta \psi_z  \right) \\
q_{soil}&=\sum_{i=1}^{nlevsoi}{q_{3,i}}=\sum_{i=1}^{nlevsoi}{k_{3,i}\cdot RAI\cdot\left(\psi_{soil,i}-\psi_{root} + \Delta\psi_{z,i} \right)}
\end{aligned}
\end{aligned}$$

We\'re looking to find the vector $\psi$ that fits with soil and
atmospheric forcings while satisfying water flow continuity. Due to the
model non-linearity, we use a linearized explicit approach, iterating
with Newton\'s method. The initial guess is the solution for $\psi$
(vector) from the previous time step. The general framework, from
iteration [m]{.title-ref} to [m+1]{.title-ref} is:

$$\begin{aligned}
q^{m+1}=q^m+\dfrac{\delta q}{\delta\psi}\Delta\psi \\
\psi^{m+1}=\psi^{m}+\Delta\psi
\end{aligned}$$

So for our first flux balance equation, at iteration [m+1]{.title-ref},
we have:

$$E_{sun}^{m+1}=q_{1a}^{m+1}$$

Which can be linearized to:

$$E_{sun}^{m}+\dfrac{\delta E_{sun}}{\delta\psi}\Delta\psi=q_{1a}^{m}+\dfrac{\delta q_{1a}}{\delta\psi}\Delta\psi$$

And rearranged to be:

$$\dfrac{\delta q_{1a}}{\delta\psi}\Delta\psi-\dfrac{\delta E_{sun}}{\delta\psi}\Delta\psi=E_{sun}^{m}-q_{1a}^{m}$$

And for the other 3 flux balance equations:

$$\begin{aligned}
\begin{aligned}
\dfrac{\delta q_{1b}}{\delta\psi}\Delta\psi-\dfrac{\delta E_{sha}}{\delta\psi}\Delta\psi&=E_{sha}^{m}-q_{1b}^{m} \\
\dfrac{\delta q_2}{\delta\psi}\Delta\psi-\dfrac{\delta q_{1a}}{\delta\psi}\Delta\psi-\dfrac{\delta q_{1b}}{\delta\psi}\Delta\psi&=q_{1a}^{m}+q_{1b}^{m}-q_2^{m} \\
\dfrac{\delta q_{soil}}{\delta\psi}\Delta\psi-\dfrac{\delta q_2}{\delta\psi}\Delta\psi&=q_2^{m}-q_{soil}^{m}
\end{aligned}
\end{aligned}$$

Putting all four together in matrix form:

$$\left[ \begin {array}{c}
\dfrac{\delta q_{1a}}{\delta\psi}-\dfrac{\delta E_{sun}}{\delta\psi} \cr
\dfrac{\delta q_{1b}}{\delta\psi}-\dfrac{\delta E_{sha}}{\delta\psi} \cr
\dfrac{\delta q_2}{\delta\psi}-\dfrac{\delta q_{1a}}{\delta\psi}-\dfrac{\delta q_{1b}}{\delta\psi} \cr
\dfrac{\delta q_{soil}}{\delta\psi}-\dfrac{\delta q_2}{\delta\psi}
\end {array} \right]
\Delta\psi=
\left[ \begin {array}{c}
E_{sun}^{m}-q_{1a}^{m} \cr
E_{sha}^{m}-q_{1b}^{m} \cr
q_{1a}^{m}+q_{1b}^{m}-q_2^{m} \cr
q_2^{m}-q_{soil}^{m}
\end {array} \right]$$

Now to expand the left-hand side, from generic $\psi$ to all four plant
water potential nodes, noting that many derivatives are zero (e.g.
$\dfrac{\delta E_{sun}}{\delta\psi_{sha}}=0$)

Introducing the notation: $A\Delta\psi=b$

$$\Delta\psi=\left[ \begin {array}{c}
\Delta\psi_{sunleaf} \cr
\Delta\psi_{shadeleaf} \cr
\Delta\psi_{stem} \cr
\Delta\psi_{root}
\end {array} \right]$$

$$A=
\left[ \begin {array}{cccc}
\dfrac{\delta q_{1a}}{\delta \psi_{sun}}-\dfrac{\delta E_{sun}}{\delta \psi_{sun}}&0&\dfrac{\delta q_{1a}}{\delta \psi_{stem}}&0\cr
0&\dfrac{\delta q_{1b}}{\delta \psi_{sha}}-\dfrac{\delta E_{sha}}{\delta \psi_{sha}}&\dfrac{\delta q_{1b}}{\delta \psi_{stem}}&0\cr
-\dfrac{\delta q_{1a}}{\delta \psi_{sun}}&
-\dfrac{\delta q_{1b}}{\delta \psi_{sha}}&
\dfrac{\delta q_2}{\delta \psi_{stem}}-\dfrac{\delta q_{1a}}{\delta \psi_{stem}}-\dfrac{\delta q_{1b}}{\delta \psi_{stem}}&
\dfrac{\delta q_2}{\delta \psi_{root}}\cr
0&0&-\dfrac{\delta q_2}{\delta \psi_{stem}}&\dfrac{\delta q_{soil}}{\delta \psi_{root}}-\dfrac{\delta q_2}{\delta \psi_{root}}
\end {array} \right]$$

$$b=
\left[ \begin {array}{c}
E_{sun}^{m}-q_{b1}^{m} \cr
E_{sha}^{m}-q_{b2}^{m} \cr
q_{b1}^{m}+q_{b2}^{m}-q_{stem}^{m} \cr
q_{stem}^{m}-q_{soil}^{m}
\end {array} \right]$$

Now we compute all the entries for $A$ and $b$ based on the soil
moisture and maximum transpiration forcings and can solve to find:

$$\Delta\psi=A^{-1}b$$

$$\psi_{m+1}=\psi_m+\Delta\psi$$

We iterate until $b\to 0$, signifying water flux balance through the
system. The result is a final set of water potentials ( $\psi_{root}$,
$\psi_{xylem}$, $\psi_{shadeleaf}$, $\psi_{sunleaf}$) satisfying
non-divergent water flux through the system. The magnitude of the water
flux is driven by soil matric potential and unstressed ( $\beta_t=1$)
transpiration.

We use the transpiration solution (corresponding to the final solution
for $\psi$) to compute stomatal conductance. The stomatal conductance is
then used to compute $\beta_t$.

$$\beta_{t,sun} = \dfrac{g_{s,sun}}{g_{s,sun,\beta_t=1}}$$

$$\beta_{t,shade} = \dfrac{g_{s,shade}}{g_{s,shade,\beta_t=1}}$$

The $\beta_t$ values are used in the Photosynthesis module (see section
`Photosynthesis`{.interpreted-text role="numref"}) to apply water
stress. The solution for $\psi$ is saved as a new variable (vegetation
water potential) and is indicative of plant water status. The
soil-to-root fluxes $\left( q_{3,1},q_{3,2},\mbox{...},q_{3,n}\right)$
are used as the soil transpiration sink in the Richards\' equation
subsurface flow equations (see section `Soil Water`{.interpreted-text
role="numref"}).

### Flow Diagram of Leaf Flux Calculations: {#Flow Diagram of Leaf Flux Calculations}

PHS runs nested in the loop that solves for sensible and latent heat
fluxes and temperature for vegetated surfaces (see section
`Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces`{.interpreted-text
role="numref"}). The scheme iterates for convergence of leaf temperature
($T_l$), transpiration water stress ($\beta_t$), and intercellular CO2
concentration ($c_i$). PHS is forced by maximum transpiration (absent
water stress, $\beta_t=1$), whereby we first solve for assimilation,
stomatal conductance, and intercellular CO2 with $\beta_{t,sun}$ and
$\beta_{t,shade}$ both set to 1. This involves iterating to convergence
of $c_i$ (see section `Photosynthesis`{.interpreted-text
role="numref"}).

Next, using the solutions for $E_{sun,max}$ and $E_{shade,max}$, PHS
solves for $\psi$, $\beta_{t,sun}$, and $\beta_{t,shade}$. The values
for $\beta_{t,sun}$, and $\beta_{t,shade}$ are inputs to the
photosynthesis routine, which now solves for attenuated photosynthesis
and stomatal conductance (reflecting water stress). Again this involves
iterating to convergence of $c_i$. Non-linearities between $\beta_t$ and
transpiration require also iterating to convergence of $\beta_t$. The
outermost level of iteration works towards convergence of leaf
temperature, reflecting leaf surface energy balance.

::: {#Figure PHS Flow Diagram}
![Flow diagram of leaf flux calculations](phs_iteration_schematic.*)
:::
