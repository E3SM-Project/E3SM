# Sediment Transport

Our treatment of the transport of sediment is based on a model developed by
Hairsine and Rose that accounts for size-selective sediment transport using a
particle size distribution.

## 2-D Hairsine-Rose (H-R) Equations

The H-R equations use a particle size distribution consisting of a set of
$P$ discrete particle/sediment size classes $p = 1, 2, ..., P$. Each size class
is represented by a sediment concentration $c_p$ and the mass $M_p$ of the layer
deposited by size class $p$ on the bed floor.

Each sediment concentration $c_p$ evolves in time according to its own transport
equation

\begin{equation}
\frac{\partial (hc)_p}{\partial t} + \nabla\cdot(h c_p \vec{u}) = e_p + e_{rp} + r_p + r_{rp} - d_p \tag{1}\label{1}
\end{equation}

where

* $h$ is the water height, as in the [shallow water equations](swe.md)
* $\vec{u} = (u, v)$ is the water flow velocity, along with which sediments are carried
* $e_p$ and $e_{rp}$ are the _rainfall-driven detachment_ and _re-detachment rates_
* $r_p$ and $r_{rp}$ are the _flow-induced entrainment and re-entrainment rates_
* $d_p$ is a _deposition rate_, expressed as mass per unit area per unit time

and $\nabla\cdot\vec{F} = \partial F_x/\partial x + \partial F_y/\partial y$ is
the 2D divergence of the spatial vector $\vec{F}$.

The deposited layer mass $M_p$ for each size class accumulates according to an
ordinary differential equation involving its deposition, re-detachment, and
re-entrainment rates:

\begin{equation}
\frac{\partial M_p}{\partial t} = d_p - e_{rp} - r_{rp}\tag{2}\label{2}
\end{equation}

All size classes deposit their layers to the bed floor, changing the bed
elevation according to the ordinary differential equation

\begin{equation}
(1-\beta)\rho_{s}\frac{\partial z}{\partial t} = \sum_{p=1}^{P}(d_p - e_p - e_{rp} - r_p - r_{rp})\tag{3}\label{3}
\end{equation}

where 

* $\beta$ is the porosity of the soil in its original state
* $\rho_s$ is the density of solid sediment, assumed to be the same for all size classes.

### Source terms

[Hairsine and Rose, 1992] specify forms for each of the source terms appearing in
the H-R equations above.

#### Rainfall-driven detachment and re-detachment rates

\begin{eqnarray}
e_p &=& F_w (1 - H) f_p a_0 R \\
e_{rp} &=& F_w H \frac{M_p}{M_t} a_d R \tag{4}\label{4}
\end{eqnarray}

where

* $f_p$ is the time-dependent ratio of the fraction of sediment in size class $p$
  to its proportion in the soil's original state (i.e. $f_p(0) = 1$)
* $a_0$ and $a_d$ are the detachability of uneroded and deposited soil, expressed
  in mass per unit volume
* R is the intensity of rainfall, expressed as the change in water height
  per unit time
* $M_t = \sum M_p$ is the total sediment mass in the deposited layer,
  expressed in mass per unit area
* $F_w$ is a _shield factor_ that attenuates the detachment and re-detachment
  rates under conditions where the water height is more than 3 times the
  diameter of a "typical" raindrop.
* $H = \min(M_t/(F_w M_t^*),1)$ is the proportion of shielding of the deposited
  layer, given in mass per unit area; here, $M_t^*$ is calibrated to the mass of
  deposited sediment needed to completely shield the soil in its original state.

The shield factor $F_w$ can be computed using a power law relation by [Proffitt et al. 1991]:

\begin{equation}
F_{w}=
\begin{cases}
  1             \quad & h \le h_{0} \\
  (h_{0}/h)^{b} \quad & h > h_{0}   \\
\end{cases} \tag{5}\label{5}
\end{equation}

where $h_0$ is a threshold height (typically $0.33 D_R$, with $D_R$ the mean
raindrop size).

The exponent $b$ in $\eqref{5}$ varies depending on the type of soil, and can be
obtained with a best fit using experimental data. For example, $b$ is 0.66 for
clay and 1.13 for loam.

#### Overland flow-driven entrainment and re-entrainment rates

\begin{eqnarray}
r_p &=& (1-H)f_p}\frac{F(\Omega-\Omega_{cr})}{J} \\
r_{rp} &=& H\frac{M_p}{M_{t}}\frac{F(\Omega - \Omega_{cr})}{(\rho_{s}-\rho_{w})gh/\rho_{s}} \tag{6}\label{6}
\end{eqnarray}

where

* $\Omega = \rho_{w}gh S_f \sqrt{u^2+v^2}$ is the _stream power_ in mass per cubic unit time,
  with $S_f = n^2 (u^2 + v^2) h^{-4/3}$
* $\Omega_{cr}$ is the _critical stream power_, below which neither soil
  entrainment or re-entrainment occur
* $F$ is the _effective fraction of excess stream power_ in entrainment or
  re-entrainment, which accounts for thermal energy dissipation
* $J$ is the _specific energy of entrainment_ in energy per unit mass, which
  indicates e.g. the energy required for soil of a given mass to be entrained
* $\rho_{w}$ is the density of water.

#### Size class deposition rate

\begin{equation}
d_p = v_p c_p \tag{7}\label{7}
\end{equation}

where $v_p$ is the _settling velocity_ of each size class with concentration
$c_p$, given as mass per unit volume. This model assumes that

* the suspended load in the water column is completely mixed in the vertical direction
* the infiltration rate does not affect size class settling velocities.

## Coupling the H-R equations with the Shallow Water Equations

Equations $\eqref{1}$ can be coupled with the [shallow water equations](swe.md)
by augment the solution vector $\mathbf{U}$ with water-height-weighted
sediment size-class concentrations:

\begin{align}
\mathbf{U} =
  \begin{bmatrix}
    h      \\[.5em]
    uh     \\[.5em]
    vh     \\
    c_1 h \\
    \vdots \\
    c_P h
  \end{bmatrix}.
\end{align}

We also augment the flux vectors $\mathbf{E}$ and $\mathbf{G}$ from the shallow
water equations with the flux terms for the sediment size class transport
equations:

\begin{align}
\mathbf{E} =
  \begin{bmatrix}
    u h                       \\[.5em]
    u^2 h + \frac{1}{2} g h^2 \\[.5em]
    u v h                     \\
    c_1 u h                   \\
    \vdots                    \\
    c_P u h
  \end{bmatrix},
\end{align}

\begin{align}
\mathbf{G} =
  \begin{bmatrix}
    v h                       \\[.5em]
    u v h                     \\[.5em]
    v^2 h + \frac{1}{2} g h^2 \\
    c_1 v h                   \\
    \vdots                    \\
    c_P v h
  \end{bmatrix}.
\end{align}

Additionally, we augment the shallow water equation source vector $\mathbf{S}$
with the (re)attachment, (re)entrainment, and deposition terms:

\begin{align}
\mathbf{S} =
  \begin{bmatrix}
    Q
    -g h\frac{\partial z}{\partial x} - C_D u\sqrt{u^2 + v^2} \\[.5em]
    -g h\frac{\partial z}{\partial y} - C_D v\sqrt{u^2 + v^2} \\
    e_1 + e_{r1} + r_1 + r_{r1} - d_1                                 \\
    \vdots                                                          \\
    e_P + e_{rP} + r_P + r_{rP} - d_P
  \end{bmatrix}.
\end{align}


Finally, to represent the deposition of mass on the bed floor, we define a
_deposited mass vector_ $\mathbf{M}$ and a _net deposition vector_ $\mathbf{D}$:

\begin{align}
\mathbf{M} =
  \begin{bmatrix}
    M_1    \\[.5em]
    \vdots   \\
    M_P
  \end{bmatrix},
\end{align}

\begin{align}
\mathbf{D} =
  \begin{bmatrix}
    d_1-e_{r1}-r_{r1} \\[.5em]
    \vdots              \\
    d_P-e_{rI}-r_{rP}
  \end{bmatrix}.
\end{align}

With these augmented and additional quantities, we can merge the H-R equations
with the shallow water equations:

\begin{eqnarray}
\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{E}}{\partial x} + \frac{\partial \mathbf{G}}{\partial y} &=& \mathbf{S}\\
\frac{\partial \mathbf{M}}{\partial t} &=& \mathbf{D} \tag{8}\label{8}
\end{eqnarray}

As in the case of the shallow water equations by themselves, we can form a
multicomponent spatial flux vector $\mathbf{\vec{F}} = (\mathbf{E}, \mathbf{G})$
to better accommodate our numerical treatment.

## TELEMAC/GAIA source terms

The TELEMAC/GAIA sediment transport model solves the coupled H-R/shallow water
equations, but uses simplified source terms in the size-class specific transport
equations:

\begin{equation}
\frac{\partial (h c_p)}{\partial t} + \nabla\cdot\left(h c_p \vec{u}\right) = E_p - D_p \tag{9} \label{9}
\end{equation}

with source terms $E_p$ and $D_p$ representing size-class-specific erosion and
deposition rates, each expressed as mass per unit area per unit time.

The GAIA model calculates these erosion and deposition rates from the following
expressions for each size class $p$:

\begin{eqnarray}
E_p &=& \mathcal{M} \left( \frac{\tau_b - \tau_{ce}}{\tau_{ce}} \right) \\
D_p &=& w_p c_p \left[ 1 - \left( \frac{\tau_b}{\tau_{cd}} \right) \right] \tag{10} \label{10}
\end{eqnarray}

where

* $\mathcal{M}$ is the Krone-Partheniades erosion law constant, sometimes called
  the "erodibility coefficient"
* $w_p$ is the settling velocity for sediment class $p$
* $\tau_{ce}$ is the critical shear stress for erosion
* $\tau_{cd}$ is the critical shear stress for deposition, and
* $\tau_b = \rho_s C_D u \sqrt{u^2 + v^2}$ is the bed bottom shear stress.

## Spatial discretization

The spatial discretization for the coupled H-R/shallow water equations is very
similar to the treatment described for [the shallow water equations](swe.md),
but uses the augmented solution vector $\mathbf{U}$, flux $\mathbf{\vec{F}}$,
and source term $\mathbf{S}$, which have analogous eigenvectors that can be
used to solve the Riemann problem with the Roe method.

Defining quantities normal to the face separating two cells (or on the boundary)
with a $\parallel$ subscript and the angle $\phi$ separating the face normal
from the $x$ axis, the normal flux is

\begin{align}
\mathbf{\vec{F}} \cdot \vec{n} =
  \begin{bmatrix}
    hu_{\parallel}                                                                  \\[.5em]
    huu_{\parallel} + \frac{1}{2}gh^{2}cos \phi + \frac{1}{24}g\Delta h^{2}cos \phi \\
    hvu_{\parallel} + \frac{1}{2}gh^{2}sin \phi + \frac{1}{24}g\Delta h^{2}sin \phi \\
    hc_1 u_{\parallel}                                                             \\
    \vdots                                                                          \\
    hc_P u_{\parallel}                                                             \\
  \end{bmatrix}\tag{11}\label{11}.
\end{align}

The last terms in the second and third rows of $\eqref{11}$ are _hydrostatic
thrust correction_ terms suggested by Bradford and Sanders [2002]. These terms
balance the bed slope terms for the still water condition.

The fluxes at the interface between cells can be approximated with Roe's method:

\begin{equation}
\mathbf{F} \cdot \mathbf{n} \approx \mathbf{F}_{\parallel,f} =
\frac{1}{2} \left(\mathbf{F}_{\parallel,L} + \mathbf{F}_{\parallel,R}-\mathbf{\hat{R}} |\mathbf{\hat{\Lambda}}| \mathbf{\Delta\hat{V}} \right)
\end{equation}

where the subscript $f$ annotates the interface between two adjacent cells,
$L$ and $R$ indicate the "left" and "right" states for the interface, and
$\Delta$ denotes the difference in quantities across the interface. The terms
$\mathbf{\hat{R}}$ and $\mathbf{\hat{\Lambda}}$ are the right eigenvector and the
eigenvalue of the Jacobian of $\mathbf{F}_{\parallel}$, and
$\mathbf{\Delta}\mathbf{\hat{V}}=\hat{L}\Delta U$ denotes the wave strength,
with $\hat{L}$ the left eigenvector of the Jacobian of $\mathbf{F}_{\parallel}$.

\begin{align}
\mathbf{\hat{R}} =
  \begin{bmatrix}
    1                         & 0           & 1                         & 0      & \ldots & 0      \\
    \hat{u}-\hat{a}cos \theta & -sin \theta & \hat{u}+\hat{a}cos \theta & 0      & \ldots & 0      \\
    \hat{v}-\hat{a}sin \theta &  cos \theta & \hat{v}+\hat{a}sin \theta & 0      & \ldots & 0      \\
    \hat{c_1}               & 0           & \hat{c_1}               & 1      & \ldots & 0      \\
    \vdots                    & \vdots      & \vdots                    & \vdots & \ddots & \vdots \\ 
    \hat{c_P}               & 0           & \hat{c_P}               & 0      & \ldots & 1      \\
  \end{bmatrix}
\end{align}

\begin{align}
\mathbf{\hat{\Lambda}} =
  \begin{bmatrix}
    |\hat{u_{\parallel}}-\hat{a}|^{*} &                   &                               &                   &        &                   \\
                                  & |\hat{u_{\parallel}}| &                               &                   &        &                   \\
                                  &                   & |\hat{u_{\parallel}}+\hat{a}|^{*} &                   &        &                   \\
                                  &                   &                               & |\hat{u_{\parallel}}| &        &                   \\
                                  &                   &                               &                   & \ddots &                   \\
                                  &                   &                               &                   &        & |\hat{u_{\parallel}}| \\
  \end{bmatrix}
\end{align}

\begin{align}
\mathbf{\Delta}\mathbf{\hat{V}}=\hat{L}\Delta U =
  \begin{bmatrix}
    \frac{1}{2} \left( \Delta h - \frac{\hat{h}\Delta u_\perp}{\hat{a}} \right) \\[.5em]
    \hat{h}u_\perp                                                          \\[.5em]
    \frac{1}{2} \left( \Delta h + \frac{\hat{h}\Delta u_\perp}{\hat{a}} \right) \\[.5em]
    (c_1 h)_{R} - (c_1 h)_{L} - \hat{c_1}(h_{R}-h_{L})                      \\[.5em]
    \vdots                                                                      \\[.5em]
    (c_P h)_{R} - (c_P h)_{L} - \hat{c_P}(h_{R}-h_{L})                      \\[.5em]
  \end{bmatrix}
\end{align}

Above,

* $a$ is the _celerity_ of a simple gravity wave, and
* $u_{\perp} = -u \sin \phi + v \cos \phi$ is the velocity perpendicular to
  the interface normal.

The quantities with a hat are Roe averages, which are calculated thus:

\begin{eqnarray}
\hat{h}   &=& \sqrt{h_L h_R}                                              \\
\hat{u}   &=& \frac{\sqrt{h_L}u_L + \sqrt{h_R}u_R}{\sqrt{h_L}+\sqrt{h_R}} \\
\hat{v}   &=& \frac{\sqrt{h_L}v_L + \sqrt{h_R}v_R}{\sqrt{h_L}+\sqrt{h_R}} \\
\hat{a}   &=& \sqrt{\frac{g}{2}(h_L + h_R)}                               \\
\hat{c_i} &=& \frac{\sqrt{h_L}c_{i,L}+\sqrt{h_R} c_{i,R}}{\sqrt{h_L}+\sqrt{h_R}}
\end{eqnarray}

The asterisks indicate that the eigenvalues $\hat{\lambda}_1=\hat{u}_{\parallel}-\hat{a}$
and $\hat{\lambda}_{3}=\hat{u}_{\parallel}+\hat{a}$ are adjusted because Roe's
method does not provide correct fluxes for critical flow:

\begin{equation}
|\hat{\lambda}_1|^{*} = \frac{\hat{\lambda}_1^{2}}{\Delta \lambda} + \frac{\Delta \lambda}{4}
\quad if \quad -\Delta \lambda /2 < \hat{\lambda}_1 < \Delta \lambda /2
\end{equation}

\begin{equation}
|\hat{\lambda}_{3}|^{*} = \frac{\hat{\lambda}_{3}^{2}}{\Delta \lambda} + \frac{\Delta \lambda}{4}
\quad if \quad -\Delta \lambda /2 < \hat{\lambda}_{3} < \Delta \lambda /2
\end{equation}

with $\Delta \lambda = 4(\lambda_{R}-\lambda_{L})$.

## References

* Hairsine, P. B., and C. W. Rose (1991). Rainfall detachment and deposition:
Sediment transport in the absence of flow-driven processes, Soil Sci. Soc. Am. J., 55(2), 320–324.
* Hairsine, P. B., and C. W. Rose (1992). Modeling water erosion due to overland flow using physical principles: 1. Sheet flow, Water Resour. Res., 28(1), 237–243.
* Kim, J., V. Y. Ivanov, and N. D. Katopodes (2013). Modeling erosion and sedimentation coupled with hydrological and overland flow processes at the watershed scale, Water Resour. Res., 49, 5134–5154, doi:10.1002/wrcr.20373.
