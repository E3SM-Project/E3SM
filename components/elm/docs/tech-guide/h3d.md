# Hybrid-3D hillslope hydrological model


### Area-Weighted Parameter Averaging

$$
r_{sub,\,top}^{default}
= \sum_{c=c_0}^{c_0+N-1}
  tmp_{rsub,\,top}(c)
  \times
  hs_{dA}(l,\,c-c_0+1)
$$

$$
zwt_{h3d}^{avg}
= \frac{1}{hs_{area}(l)}
  \sum_{c=c_0}^{c_0+N-1}
  zwt_{h3d}(c)
  \times
  hs_{dA}(l,\,c-c_0+1)
$$

$$
r_{sub,\,top,\,h3d}^{max,\,avg}
= \frac{1}{hs_{area}(l)}
  \sum_{c=c_0}^{c_0+N-1}
  r_{sub,\,top,\,h3d}^{max}(c)
  \times
  hs_{dA}(l,\,c-c_0+1)
$$

$$
fff^{avg}
= \frac{1}{hs_{area}(l)}
  \sum_{c=c_0}^{c_0+N-1}
  fff(c)
  \times
  hs_{dA}(l,\,c-c_0+1)
$$

### Default Subsurface Runoff Calculation

$$
r_{sub,\,top}^{default}
= r_{sub,\,top,\,h3d}^{max,\,avg}
  \times
  \exp\!\left(-fff^{avg} \times zwt_{h3d}^{avg}\right)
  \times
  hs_{area}(l)
$$

### Area-Weighted Storage Change

$$
\Delta S_{sat}^{tot}
= \sum_{c=c_0}^{c_0+N-1}
  \Delta S_{sat}(c)
  \times
  hs_{dA}(l,\,c-c_0+1)
$$

### Area-Weighted Drainage Sum

$$
Q_{drain}^{tot}
= \sum_{c=c_0}^{c_0+N-1}
  q_{drain}(c)
  \times
  hs_{dA}(l,\,c-c_0+1)
$$

Variable Definitions

| Variable              | Description                                          | Units |
| :-------------------- | :--------------------------------------------------- | :---- |
| `hs_dA(l,i)`          | Surface area of h3D soil column *i* in land unit *l* | m²    |
| `hs_area(l)`          | Total surface area of hillslope in land unit *l*     | m²    |
| `l`                   | Land unit index                                      | –     |
| `c₀`                  | First column index in land unit *l*                  | –     |
| `N = nh3dc_per_lunit` | Number of columns per land unit                      | –     |


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

## Finite-Difference Discretization

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

## Subsurface Runoff and Storage Change

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

Adaptive time-stepping:

$$
\Delta t_{\text{new}} = 0.5\,\Delta t_{\text{old}}
$$


## Tridiagonal Matrix Solution

Forward elimination:

$$
\begin{aligned}
\beta_1 &= b_1, \quad \tilde{r}_1 = \frac{r_1}{\beta_1}, \\
\gamma_i &= \frac{c_{i-1}}{\beta_{i-1}}, \\
\beta_i &= b_i - a_i\gamma_i, \\
\tilde{r}_i &= \frac{r_i - a_i\tilde{r}_{i-1}}{\beta_i}, \quad (i=2,\dots,N)
\end{aligned}
$$

Back substitution:

$$
h_N = \tilde{r}_N, \qquad
h_i = \tilde{r}_i - \gamma_{i+1} h_{i+1}, \quad (i=N-1,\dots,1)
$$

## Boundary Conditions

| Boundary       | Condition                             | Implementation           |
| :------------- | :------------------------------------ | :----------------------- |
| Lower (stream) | Head-dependent outflow                | Robin-type term in $r_1$ |
| Upper (divide) | No-flow ($\partial h/\partial x = 0$) | Set $c_N = 0$            |



