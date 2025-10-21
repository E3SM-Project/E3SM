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

Each land unit represents a single hillslope, which has consistent topographic and geometric properties:
- **Overall slope angle:** $\theta$ — mean hillslope angle (rad)
- **Width function:** $w(x)$ — lateral width distribution along the hillslope (m)
- **Distance function:** $x(i)$ — distance from the stream outlet to node $i$ (m)
- **Total hillslope area:** $A_{hs} = \displaystyle \int_0^L w(x)\,dx$ (m²)


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


# DrainageH3D → H3D_DRI → LateralResponse

The h3D pathway computes (a) a **default topographic drainage propensity** from SIMTOP-style scalings, (b) a **within-landunit lateral response** by solving a **1-D hillslope Dupuit–Boussinesq** problem for saturated thickness, and (c) converts the updated saturated thickness to **subsurface runoff** and **diagnostics**.

---

## 1. Area-Weighted Hillslope Inputs (per landunit \(l\))

$$
r_{sub,\,top}^{default}
=\sum_{c=c_0}^{c_0+N-1}
  tmp_{rsub,\,top}(c)\; hs_{dA}(l,\,c-c_0+1)
\qquad[\mathrm{m^3\,s^{-1}}]
$$

$$
\overline{zwt}_{h3d}
=rac{1}{hs_{area}(l)}
 \sum_{c=c_0}^{c_0+N-1}
 zwt_{h3d}(c)\; hs_{dA}(l,\,c-c_0+1)
\qquad[\mathrm{m}]
$$

$$
\overline{r_{sub,\,top,\,h3d}^{max}}
=rac{1}{hs_{area}(l)}
 \sum_{c=c_0}^{c_0+N-1}
 r_{sub,\,top,\,h3d}^{max}(c)\; hs_{dA}(l,\,c-c_0+1)
\qquad[\mathrm{m\,s^{-1}}]
$$

$$
\overline{fff}
=rac{1}{hs_{area}(l)}
 \sum_{c=c_0}^{c_0+N-1}
 fff(c)\; hs_{dA}(l,\,c-c_0+1)
\qquad[\mathrm{m^{-1}}]
$$

$$
r_{sub,\,top}^{default}
=
\overline{r_{sub,\,top,\,h3d}^{max}}\;
\exp\!\left(-\,\overline{fff}\;\overline{zwt}_{h3d}
ight)\;
hs_{area}(l)
\qquad[\mathrm{m^3\,s^{-1}}]
$$

---

## 2. Saturated Thickness and Bounds

$$
h_{sat}(c) \equiv z_{wt,\,bed}(c) - z_{wt}(c),
\qquad h_{sat}\ge 0
\qquad[\mathrm{m}]
$$

$$
z_{wt}(c) \in [\,0,\;z_{wt,\,bed}(c)\,],\quad
z_{wt,\,bed}(c)=egin{cases}
z_{ibed}(c), & 	ext{variable soil thickness}\
z_{ibed}(c)+25\,\mathrm{m}, & 	ext{otherwise}
\end{cases}
$$

$$
h_{sat}^{begin}(c)=z_{wt,\,bed}(c)-z_{wt}^{begin}(c)
$$

---

## 3. Specific Yield \(f_{	ext{drain}}\) (Brooks–Corey)

Let \(idx=\min\{jwt(c)+1,\;n_{lev\,bed}(c)\}\).

$$
f_{drain}(c)
=
\max\!\left[
0.02,\;
	heta_{sat}(c,idx)\left(
1-\left(
1+rac{1000\;\max\!\{0,\;z_{wt,\,bed}(c)-h_{sat}(c)\}}
{\psi_{sat}(c,idx)}

ight)^{-1/b(c,idx)}

ight)

ight]
\qquad[-]
$$

---

## 4. Transmissivity Along the Hillslope (Interface Form)

For hillslope node \(k=1,\dots,N\) (lower boundary \(k=1\) at stream; upper \(k=N\) at divide), with anisotropy \(f_{aniso}=100\):

$$
wK\!H(k)
=
f_{aniso}\;
hs_{w,\,nod}(l,k)\;
rac{hk_{sat}(c_k,idx_k)}{1000}\;
h_{sat}^{prev}(k)
\qquad[\mathrm{m^3\,s^{-1}}]
$$

where \(c_k=c_0+k-1\) and \(idx_k=\min\{jwt(c_k)+1,\,n_{lev\,bed}(c_k)\}\).

---

## 5. Discrete Hillslope PDE (Implicit, Tridiagonal in \(h_{sat}\))

Let \(	heta_k = h3d\_slope(c_k)\cdot\pi/180\), and geometry \(hs_{dx}(l,k)\), \(hs_{dx,\,nod}(l,k)\), \(hs_{w,\,nod}(l,k)\).

### 5.1 Interior Nodes \(k=2,\dots,N-1\)

$$
a_k
= -rac{wK\!H(k)\cos	heta_k\;\Delta t_{h3d}}
        {hs_{dx,\,nod}(l,k)\;hs_{dx}(l,k)\;hs_{w,\,nod}(l,k)},
\qquad
c_k
= -rac{wK\!H(k+1)\cos	heta_k\;\Delta t_{h3d}}
        {hs_{dx,\,nod}(l,k+1)\;hs_{dx}(l,k)\;hs_{w,\,nod}(l,k)},
$$

$$
b_k = f_{drain}(c_k) - (a_k + c_k),
$$

$$
r_k
= f_{drain}(c_k)\;h_{sat}^{old}(k)
+ rac{\Delta t_{h3d}\,\sin	heta_k}{hs_{w,\,nod}(l,k)\;hs_{dx}(l,k)}
  \Big[wK\!H(k+1) - wK\!H(k)\Big],
$$

$$
a_k\,h_{sat}(k-1) + b_k\,h_{sat}(k) + c_k\,h_{sat}(k+1) = r_k.
$$

### 5.2 Lower Boundary \(k=1\) (Stream; Head-Dependent Outflow)

$$
a_1=0,\qquad
c_1
= -rac{wK\!H(2)\cos	heta_1\;\Delta t_{h3d}}
        {hs_{dx,\,nod}(l,2)\;hs_{dx}(l,1)\;hs_{w,\,nod}(l,1)},
$$

$$
b_1=f_{drain}(c_1)-(a_1+c_1),
$$

$$
r_1=
f_{drain}(c_1)\;h_{sat}^{old}(1)
+rac{\Delta t_{h3d}}{hs_{w,\,nod}(l,1)\;hs_{dx}(l,1)}
\left[
\sin	heta_1\;wK\!H(2)
-rac{\cos	heta_1}{hs_{dx}(l,1)}
\,hs_{w,\,nod}(l,1)\,
rac{f_{aniso}\,hk_{sat}(c_1,idx_1)}{1000}\,ig(h_{sat}^{prev}(1ig)^2

ight].
$$

### 5.3 Upper Boundary \(k=N\) (Divide; No-Flow)

$$
a_N
= -rac{wK\!H(N)\cos	heta_N\;\Delta t_{h3d}}
        {hs_{dx,\,nod}(l,N)\;hs_{dx}(l,N)\;hs_{w,\,nod}(l,N)},
\qquad
c_N=0,
$$

$$
b_N=f_{drain}(c_N)-(a_N+c_N),
$$

$$
r_N=
f_{drain}(c_N)\;h_{sat}^{old}(N)
-rac{\Delta t_{h3d}\,\sin	heta_N}{hs_{w,\,nod}(l,N)\;hs_{dx}(l,N)}
\,wK\!H(N).
$$

---

## 6. Linear Solve and Iteration

### Tridiagonal System (Thomas Algorithm)

$$
a_k\,h_{sat}(k-1) + b_k\,h_{sat}(k) + c_k\,h_{sat}(k+1) = r_k,
\qquad k=1,\dots,N.
$$

### Picard Loop & Convergence

$$
\max_{kig|h_{sat}^{(m+1)}(k)-h_{sat}^{(m)}(kig|
< h_{sat}^{thres}=10^{-4}\ \mathrm{m}.
$$

If not converged in \(n_{iter}^{max}\), halve \(\Delta t_{h3d}\) and retry:

$$
\Delta t_{h3d}\leftarrow 	frac12\,\Delta t_{h3d},
\qquad
\sum \Delta t_{h3d} = \Delta t_{\mathrm{ELM}}.
$$

Recover water table:

$$
z_{wt}(c_k) = z_{wt,\,bed}(c_k) - h_{sat}(k),
\qquad
z_{wt}\in[0,\,80\,\mathrm{m}]\ 	ext{(clamp)}.
$$

---

## 7. Storage Change and Subsurface Runoff (ELM Timestep \(\Delta t\))

$$
\Delta S_{sat}(c_k)= f_{drain}(c_k)\,\Big(h_{sat}(k)-h_{sat}^{begin}(k)\Big)
\qquad[\mathrm{m}],
$$

$$
R_{sub}(c_k)= -\,\Delta S_{sat}(c_k)\qquad[\mathrm{m}],
\qquad
Q_{sub}(c_k)= rac{R_{sub}(c_k)}{\Delta t}	imes 1000
\qquad[\mathrm{mm\,s^{-1}}].
$$

$$
qflx_{drain}(c_k) = Q_{sub}(c_k),
\qquad
qflx_{drain,\,h3d}(c_k)=Q_{sub}(c_k).
$$

---

## 8. Coupling Back to Column Hydrology

### Non-h3D (SIMTOP-like)

$$
r_{sub,\,top}(c)= imped(c)\; r_{sub,\,top}^{max}(c)\;
\exp\ig(-\,fff(c)\,z_{wt}(cig)
\qquad[\mathrm{mm\,s^{-1}}].
$$

### h3D Columns

$$
r_{sub,\,top}(c) = imped(c)\; qflx_{drain,\,h3d}(c)
\qquad[\mathrm{mm\,s^{-1}}].
$$

### Total Drainage

$$
qflx_{drain}(c)= qflx_{rsub,\,sat}(c) + r_{sub,\,top}(c).
$$

---

## 9. Units and Conversions

$$
hk_{sat}\,[\mathrm{mm\,s^{-1}}] \ 
ightarrow\ \mathrm{m\,s^{-1}} 	ext{ via } /1000,
\qquad
Q\,[\mathrm{m\,s^{-1}}] \ 
ightarrow\ \mathrm{mm\,s^{-1}} 	ext{ via } 	imes 1000.
$$

---

## 10. Indices and Mappings

- Node \(k\) ↔ column \(c_k=c_0+k-1\), \(k=1,\dots,N=nh3dc\_{per\_lunit}\).
- \(jwt(c)\) is the index above the water table; \(idx=\min(jwt+1,n_{lev\,bed})\).
- \(	heta_k\) uses degrees in code; equations use radians.

