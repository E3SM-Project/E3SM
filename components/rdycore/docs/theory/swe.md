# Shallow Water Equations

The two-dimensional shallow water equations can be written in the conservative
form

$$
\frac{\partial\mathbf{U}}{\partial t} + \frac{\partial \mathbf{E}}{\partial x} + \frac{\partial \mathbf{G}}{\partial y} = \mathbf{S}_r + \mathbf{S}_b + \mathbf{S}_f \tag{1}\label{1}
$$

Here,

\begin{align}
  \mathbf{U}
  =
  \begin{bmatrix}
  h \\[.5em]
  hu \\[.5em]
  hv
  \end{bmatrix},
\end{align}

is the _solution vector_, with components

* $h$, the flow depth
* $u$, the vertically-averaged velocity in the $x$ direction
* $v$, the vertically-averaged velocity in the $y$ direction

The terms $\mathbf{E}$ and $\mathbf{G}$ on the left hand side of $\eqref{1}$ are
the _flux vectors_ in the $x$ and $y$ spatial dimensions:

\begin{align}
  \mathbf{E}
  =
  \begin{bmatrix}
  hu \\[.5em]
  hu^2 + \frac{1}{2}gh^2 \\[.5em]
  huv
  \end{bmatrix},
\end{align}

\begin{align}
  \mathbf{G}
  =
  \begin{bmatrix}
  hv \\[.5em]
  huv \\[.5em]
  hv^2 + \frac{1}{2}gh^2
  \end{bmatrix},
\end{align}

where $g$ is the acceleration due to gravity.

The three source terms on the right hand side of $\eqref{1}$ represent
contributions to $\mathbf{U}$ from

* $\mathbf{S}_r$, _net runoff (water) production_
* $\mathbf{S}_b$, _bed elevation slope_
* $\mathbf{S}_f$, _bed friction roughness_.

These source terms are

\begin{align}
  \mathbf{S}_r
  =
  \begin{bmatrix}
  Q \\[.5em]
  0 \\[.5em]
  0
  \end{bmatrix},
\end{align}

\begin{align}
  \mathbf{S}_b
  =
  \begin{bmatrix}
  0 \\[.5em]
  -gh\frac{\partial z}{\partial x} \\[.5em]
  -gh\frac{\partial z}{\partial y}
  \end{bmatrix},
\end{align}

and

\begin{align}
  \mathbf{S}_f
  =
  \begin{bmatrix}
  0 \\[.5em]
  - C_D u \sqrt{u^2 + v^2} \\[.5em]
  - C_D v \sqrt{u^2 + v^2}
  \end{bmatrix}.
\end{align}

It is sometimes convenient to interpret the terms involving $\mathbf{E}$ and
$\mathbf{G}$ as the (two-dimensional) divergence of a multi-component spatial
vector called the _flux function_ $\mathbf{\vec{F}}$.

$$
\frac{\partial\mathbf{U}}{\partial t} + \vec{\nabla}\cdot\mathbf{\vec{F}}(\mathbf{U}) = \mathbf{S}(\mathbf{U})\tag{2}\label{2}
$$

where $\mathbf{\vec{F}} = (\mathbf{F}_x, \mathbf{F}_y) = (\mathbf{E}, \mathbf{G})$.

We have written $\mathbf{\vec{F}}$ and $\mathbf{S}$ in $\eqref{2}$ in a way that
emphasizes that these functions depend on the solution vector $\mathbf{U}$.

## Spatial Discretization

We can rewrite the shallow water equations in a form more convenient for
numerical treatment by defining a computational domain $\Omega$ bounded by a
piecewise linear closed curve $\Gamma = \partial\Omega$.

We create a discrete representation by partitioning $\Omega$ into disjoint
cells, with $\Omega_i$ representing cell $i$. The boundary of cell $i$, written
$\partial\Omega_i$, is the set of faces separating it from its neighboring
cells. Using this notation, we obtain a discrete set of equations for the
solution in cell $i$ by integrating $\eqref{1}$ over $\Omega_i$ and using
Green's theorem:

\begin{eqnarray}
\frac{\partial}{\partial t} \int_{\Omega_i} \mathbf{U} d\Omega_i +
\int_{\Omega_i} \left[ \frac{\partial\mathbf{E}}{\partial x} +
\frac{\partial\mathbf{G}}{\partial y} \right] d\Omega_i &=&
\int_{\Omega_i} \mathbf{S} d\Omega_i \nonumber\\
\frac{\partial}{\partial t} \int_{\Omega_i} \mathbf{U} d\Omega_i +
\oint_{\partial\Omega_i} \left( \mathbf{E}~dy - \mathbf{G}~dx \right) &=&
\int_{\Omega_i} \mathbf{S} d\Omega_i \tag{3}\label{3}
\end{eqnarray}

This equation can be used to approximate discontinuous flows, because all
quantities appear under integrals. By contrast, $\eqref{1}$ cannot be used
where derivatives of $\mathbf{U}$ don't exist.

We can interpret the line integral in $\eqref{3}$ in terms of the flux
$\mathbf{\vec{F}} = (\mathbf{F}_x, \mathbf{F}_y)$ between a cell $i$ and its
neighboring cells.

$$
 \frac{\partial}{\partial t} \int_{\Omega_i} \mathbf{U} d\Omega_i +
\oint_{\partial\Omega_i} \mathbf{\vec{F}} \cdot \vec{n}~dl =
\int_{\Omega_i} \mathbf{S} d\Omega_i \tag{4}\label{4}
$$

Here, we have defined a unit normal vector $\vec{n} = (n_x, n_y)$
pointing outward along the cell boundary $\partial\Omega_i$. $\eqref{4}$ is a
"surface integral" with a differential arc length $dl$ integrated over the
boundary of cell $i$. One obtains this surface integral by integrating $\eqref{2}$
over the domain $\Omega$ and applying the (two-dimensional) divergence theorem
to the flux term.

In the rest of this section, we use the flux form $\eqref{4}$ of the shallow
water equations.

We can obtain a finite volume method for these equations by defining
_horizontally-averaged_ quantities for flow depth and velocities:

\begin{eqnarray}
h_i &=& \frac{1}{A_i}\int_{\Omega_i} h d\Omega_i \\
u_i &=& \frac{1}{A_i}\int_{\Omega_i} u d\Omega_i \\
v_i &=& \frac{1}{A_i}\int_{\Omega_i} v d\Omega_i
\end{eqnarray}

where $A_i = \int_{\Omega_i}d\Omega_i$ is the area enclosed within cell $i$. We
also introduce the _horizontally-averaged solution vector_

$$
\mathbf{U}_i = [h_i, h_i u_i, h_i v_i]^T
$$

and the _horizontally-averaged source vector_

$$
\mathbf{S}_i = \frac{1}{A_i}\int_{\Omega_i} \left(\mathbf{S}_r + \mathbf{S}_b + \mathbf{S}_f\right) d\Omega_i.
$$

Finally, we define the _face-averaged normal flux vector_ between cell $i$ and
an adjoining cell $j$:

$$
\mathbf{F}_{ij} = \frac{1}{l_{ij}}\int_{\partial\Omega_i\bigcap\partial\Omega_j}\mathbf{\vec{F}}\cdot\vec{n}~dl \tag{5}\label{5}
$$

where $l_{ij}$ is the length of the face connecting cells $i$ and $j$.

With these definitions, the shallow water equations in cell $i$ are

$$
\frac{\partial\mathbf{U}_i}{\partial t} + \sum_j\mathbf{F}_{ij} l_{ij} = \mathbf{S}_i, \tag{6}\label{6}
$$

where the index $j$ in each term of the sum refers to a neighboring cell of cell $i$.

### Boundary conditions

To incorporate boundary conditions, we partition the domain boundary $\Gamma$
into disjoint line segments, each of which represents a _boundary face_
$\Gamma_i$. Every cell $j$ that touches the boundary $\Gamma$ has at least one
boundary face. Such a cell is a _boundary cell_. The boundary $\Gamma$ consists
entirely of faces of boundary cells.

In dealing with boundary conditions, we must distinguish between the faces a
boundary cell $i$ _does_ and _does not_ share with the boundary $\Gamma$:

\begin{eqnarray}
\frac{\partial\mathbf{U}_i}{\partial t} +
\sum_{j: \partial\Omega_j\subset\Gamma}\mathbf{F}_{ij}^{\Gamma} l_{ij} &+&
\sum_{j: \partial\Omega_j\not\subset\Gamma}\mathbf{F}_{ij} l_{ij}
&= \mathbf{S}_i. \tag{7}\label{7} \\
\text{(boundary)}& &\text{(interior)} &
\end{eqnarray}

To enforce boundary conditions, we must compute the effective boundary fluxes
$\mathbf{F}_{ij}^{\Gamma}$ that appear in the first sum in $\eqref{7}$. These
boundary fluxes have specific forms depending on their respective boundary
conditions.

### Evaluating normal fluxes

We have reduced the spatial discretization of our finite volume method to the
calculation of normal fluxes between neighboring cells and on boundary faces.
The normal flux function is

\begin{equation}
\mathbf{\vec{F}}\cdot\vec{n} = \mathbf{E} n_x + \mathbf{G} n_y.
\end{equation}

We can evaluate $\mathbf{F}_{ij} \approx \mathbf{\vec{F}}\cdot\vec{n}$, the
approximate normal flux at a face shared shared by _interior cells_ $i$ and $j$,
by solving the relevant Riemann problem using Roe's method. If we designate
$\phi$ as the angle between $\vec{n}$ and the $x$ axis and adopt $i$ and $j$
subscripts for quantities respectively in cells $i$ and $j$, we can approximate
the normal flux by the expression

\begin{equation}
\mathbf{F}_{ij} = \frac{1}{2} \left( \mathbf{\vec{F}}_i + \mathbf{\vec{F}}_j - \mathbf{\hat{R}} |\mathbf{\hat{\Lambda}| \mathbf{\Delta}\hat{V}} \right)
\end{equation}

where $\Delta f$ is the variation of the quantity of $f$ along a face. In
particular,

\begin{align}
  \mathbf{R}
  =
  \begin{bmatrix}
  1 & 0 & 1  \\[.5em]
  \hat{u} - \hat{a}\cos\phi & -\sin\phi & \hat{u} + \hat{a}\cos\phi  \\[.5em]
  \hat{v} - \hat{a}\sin\phi &  \cos\phi & \hat{v} + \hat{a}\sin\phi
  \end{bmatrix}
\end{align}

\begin{align}
  \mathbf{\Delta\hat{V}}
  =
  \begin{bmatrix}
  \frac{1}{2} \left( \Delta h - \frac{\hat{h}\Delta u_\parallel}{\hat{a}} \right) \\[.5em]
  \hat{h}u_\perp \\[.5em]
  \frac{1}{2} \left( \Delta h + \frac{\hat{h}\Delta u_\parallel}{\hat{a}} \right)
  \end{bmatrix}
\end{align}

\begin{align}
  |\mathbf{\hat{\Lambda}}|
  =
  \begin{bmatrix}
  | \hat{u}_\parallel - \hat{a} |^* & 0 & 0  \\[.5em]
  0                                 & |\hat{u}_\parallel| & 0 \\[.5em]  
  0                                 &                     & | \hat{u}_\parallel + \hat{a} |^* 
  \end{bmatrix}
\end{align}

where

\begin{eqnarray}
  \hat{h} & = & \sqrt{h_i h_j} \\
  \hat{u} & = & \frac{ \sqrt{h_i} \vec{u}_i + \sqrt{h_j} \vec{u}_j}{ \sqrt{h_i} + \sqrt{h_j}} \\
  \hat{v} & = & \frac{ \sqrt{h_i} \vec{v}_i + \sqrt{h_j} \vec{v}_j}{ \sqrt{h_i} + \sqrt{h_j}} \\
  \hat{a} & = & \sqrt{\frac{g}{2} \left( h_i + h_j \right)},
\end{eqnarray}

$\Delta f = f_j - f_i$ is the change in the quantity $f$ moving from cell
$i$ to $j$, and $w_\parallel = \vec{w}\cdot\vec{n}$ and
$w_\perp = |\vec{w} - w_\parallel \vec{n}|$ for any vector $\vec{w}$.

We have used asterisks in the expression for $|\mathbf{\hat{\Lambda}}|$ to
indicate that the eigenvalues 
$\hat{\lambda}_1 = \hat{u}_\perp - \hat{a}$ and
$\hat{\lambda}_3 = \hat{u}_\perp + \hat{a}$ 
must be adjusted, since Roe's method does not provide the correct flux for
critical flow.

\begin{eqnarray}
  |\hat{\lambda}|_1 &=& \frac{\hat{\lambda}^2_1}{\Delta \lambda} + \frac{\Delta \lambda}{4} \mbox{$~$ if $~ -\Delta \lambda/2 < \hat{\lambda}_1 < \Delta \lambda/2$} \\
  |\hat{\lambda}|_3 &=& \frac{\hat{\lambda}^2_2}{\Delta \lambda} + \frac{\Delta \lambda}{4} \mbox{$~$ if $~ -\Delta \lambda/2 < \hat{\lambda}_3 < \Delta \lambda/2$}
\end{eqnarray}
 
### Source terms

Recall that there are three source terms $\mathbf{S}_r, $\mathbf{S}_b$, and
$\mathbf{S}_f$. In this section, we write the source terms for the momentum
vector as $\mathbf{S}_b = [0, \mathbf{\vec{S}}_b]^T$ and $\mathbf{S}_f = [0, \mathbf{\vec{S}}_f]^T$
because their second and third components correspond to the spatial components
of the momentum vector.

#### Net runoff production $\mathbf{S}_r$

This source term contributes only to height of the water, and is expressed as
$\mathbf{S}_r = [Q_r, 0, 0]^T$, where $Q_r$, _the net runoff production_,
with units of water height per unit time, is

* a constant
* a spatially homogeneous time-dependent function $Q_r(t)$, or
* a spatially heterogeneous time-dependent function $Q_r(x, y, t)$.

In any case, we can approximate the integral of this term using the mean value
theorem of calculus:

\begin{equation}
\int_{\Omega_i} S_r~d\Omega_i = \int_{\Omega_i} [Q_r, 0, 0]^T~d\Omega \approx [Q_r A_i, 0, 0]^T,
\end{equation}

where $A_i$ is the area of cell $i$.

#### Bed elevation slope term $\mathbf{S}_b$

This term represents the force of gravity on the water and can be
approximated as

\begin{equation}
\int_{\Omega_i} \mathbf{\vec{S}}_b~d\Omega_i = \int_{\Omega_i} -gh\nabla z~d\Omega_i \approx -gh\left(\overline{\frac{\partial z}{\partial x}}, \overline{\frac{\partial z}{\partial y}}\right) A_i
\end{equation}

where $\nabla z = (\partial z/\partial x, \partial z/\partial y)$ is the two-dimensional
gradient of the bed elevation function $z$.

For a triangular grid cell,

\begin{eqnarray}
\overline{\frac{\partial z}{\partial x}} &=& \frac{(y_2 - y_0)(z_1 - z_0) - (y_1 - y_0)(z_2 - z_0)}{(y_2 - y_0)(x_1 - x_0) - (y_1 - y_0)(x_2 - x_0)} \\
\overline{\frac{\partial z}{\partial y}} &=& \frac{(x_2 - x_0)(z_1 - z_0) - (x_1 - x_0)(z_2 - z_0)}{(x_2 - y_0)(x_1 - x_0) - (x_1 - y_0)(y_2 - y_0)}.
\end{eqnarray}

#### Bed friction roughness term $\mathbf{S}_f$

Like the runoff term, this term involves only quantities within a single cell
and can be approximated by the mean value theorem:

\begin{equation}
\int_{\Omega_i} \mathbf{\vec{S}}_f~d\Omega_i = \int_{\Omega_i} C_D \vec{u} \sqrt{u^2 + v^2}~d\Omega_i \approx C_D \vec{u} \sqrt{u^2 + v^2} A_i
\end{equation}

where $\vec{u} = (u, v)$ is the flow velocity vector. The $x$ and $y$ spatial
components contribute to the second and third vector components of the integral
of $\mathbf{S}_r$, respectively.

## Temporal Discretization

The above spatial discretization produces a "semi-discrete" system of equations
that can be integrated using various methods for solving systems of ordinary
differential equations. The following methods of integration are provided by
PETSc and supported for RDycore:

* [Forward Euler](https://petsc.org/release/manualpages/TS/TSEULER/)

## References

* [Bradford, S. F., & Sanders, B. F. (2002). Finite-volume model for shallow-water flooding of arbitrary topography. Journal of hydraulic engineering, 128(3), 289-298.](https://ascelibrary.org/doi/10.1061/%28ASCE%290733-9429%282002%29128%3A3%28289%29)

* [Kim, J., Warnock, A., Ivanov, V. Y., & Katopodes, N. D. (2012).
Coupled modeling of hydrologic and hydrodynamic processes including
overland and channel flow. Advances in water resources, 37, 104-126.](https://www.sciencedirect.com/science/article/pii/S0309170811002211?via%3Dihub)


