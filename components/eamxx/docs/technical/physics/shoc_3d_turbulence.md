# SHOC-HOMME 3D Turbulence Coupling

This note documents the EAMxx 3D-turbulence path that couples SHOC and HOMME.
It has two main pieces:

1. SHOC diagnoses horizontal eddy diffusivities and passes them to HOMME so the
   dycore can apply horizontal turbulent diffusion.
2. HOMME provides horizontal velocity-gradient information back to SHOC so SHOC
   can use a three-dimensional shear-production term in the TKE equation.

## Horizontal Eddy Diffusivity Passed to HOMME

When `ctl_nl::do_3d_turbulence=true`, SHOC computes horizontal eddy
diffusivities on the physics grid after the main SHOC update. The formulation
is

$$
K_{H,h} = C_{kh,h}\,L_h\,\sqrt{\mathrm{TKE}},
\qquad
K_{M,h} = C_{km,h}\,L_h\,\sqrt{\mathrm{TKE}},
$$

where the horizontal length scale is

$$
L_h = \sqrt{\Delta x\,\Delta y}.
$$

In the current implementation, SHOC obtains `\Delta x` and `\Delta y` from the
grid geometry and uses the runtime parameters
`shoc::coeff_kh_horiz` and `shoc::coeff_km_horiz`, both of which default to
`0.1`.

After SHOC computes `eddy_diff_heat_horiz` and `eddy_diff_mom_horiz`, EAMxx
forwards those fields to HOMME:

1. On the standard physics-grid path, the physics-to-dynamics remapper
   registers the SHOC fields with HOMME helper fields `Kh_dyn` and `Km_dyn`.
2. On pgN finite-volume physics grids, the same fields are supplied to
   `GllFvRemap::run_fv_phys_to_dyn`, which remaps them into HOMME's dynamics
   representation.
3. HOMME then binds the remapped fields into its derived-state turbulence slots
   `m_turb_diff_heat` and `m_turb_diff_mom`, which are consumed by the dycore's
   horizontal turbulent-diffusion machinery.

This split keeps the turbulence closure in SHOC while leaving the actual
horizontal diffusion operator in HOMME, where the horizontal discretization and
metric terms already live.

## Three-Dimensional Shear Production of TKE

In the legacy one-dimensional SHOC path, shear production uses only vertical
shear of the horizontal wind. In the 3D-turbulence path, SHOC instead forms the
full local strain invariant and uses it in the TKE tendency.

The shear-production contribution is

$$
\begin{align}
  -P_s
  &= 2K_M \Bigg[
      \left(\frac{\partial u}{\partial x}\right)^{2}
    + \left(\frac{\partial v}{\partial y}\right)^{2}
    + \left(\frac{\partial w}{\partial z}\right)^{2} \notag\\
  &\quad
    + \tfrac{1}{2}\Big(
        \left(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\right)^{2}
      + \left(\frac{\partial u}{\partial z}+\frac{\partial w}{\partial x}\right)^{2}
      + \left(\frac{\partial v}{\partial z}+\frac{\partial w}{\partial y}\right)^{2}
      \Big)
  \Bigg].
\end{align}
$$

In code, SHOC stores the bracketed strain invariant as
`tke_shear_strain3d = 2 S_{ij}S_{ij}` and then advances TKE with

$$
\text{shear production} = C_k K_M \left(2 S_{ij}S_{ij}\right),
$$

with `C_k = 0.1`, consistent with the existing SHOC TKE closure.

The calculation is split across HOMME and SHOC:

1. After each HOMME dynamics step, HOMME computes horizontal derivatives of the
   local Cartesian velocity components on the dynamics grid.
2. Those derivatives are projected back to the local basis to form the six
   horizontal pieces of the velocity-gradient tensor:
   `du/dx`, `du/dy`, `dv/dx`, `dv/dy`, `dw/dx`, and `dw/dy`.
3. EAMxx remaps those six components onto the physics grid as
   `tke_shear_strain3d_components`.
4. Inside SHOC, vertical derivatives `du/dz`, `dv/dz`, and `dw/dz` are computed
   on interface levels and then interpolated to midpoint levels, matching
   SHOC's native staggered-grid treatment.
5. SHOC assembles the full symmetric strain tensor and forms the invariant
   `2 S_{ij}S_{ij}`, which is then used in the TKE equation.

The same assembled 3D strain is also converted back to a midpoint shear
magnitude for SHOC's cold-surface fallback in the eddy-diffusivity calculation,
so the 3D path remains consistent with the rest of the closure.
