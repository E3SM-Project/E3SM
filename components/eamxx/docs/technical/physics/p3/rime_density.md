# P3 Rime Density

This document describes the P3 routine `calc_rime_density`. The routine
computes the mass-weighted mean ice fall speed used in cloud-liquid riming and
the density assigned to newly collected cloud liquid. It also describes the
physical-property tests that verify the implementation.

## Purpose

The routine computes:

- `vtrmi1`: mass-weighted mean ice fall speed used in cloud-liquid riming
- `rho_qm_cloud`: density assigned to newly collected cloud liquid

These outputs are used during prognostic ice updates to determine how collected
cloud water changes both the rimed ice mass and the rimed ice volume.

## P3 Closure Formulation

P3 predicts rimed ice mass and rimed ice volume separately. The bulk rime
density is diagnosed from their ratio,

$$
\rho_{\mathrm{rime}} = \frac{q_{\mathrm{rime}}}{b_{\mathrm{rime}}},
$$

where $q_{\mathrm{rime}}$ is rime mass and $b_{\mathrm{rime}}$ is rime volume.
In the implementation this corresponds to the prognostic variables `qm` and
`bm`, with the bulk density later diagnosed in `calc_bulk_rho_rime` and then
bounded for consistency before being used by lookup-table ice properties.

The `calc_rime_density` routine provides the incremental density closure for
newly collected cloud liquid. For cloud-liquid riming, the added rime volume is

$$
\Delta b_{\mathrm{rime,cloud}} = \frac{\Delta q_{\mathrm{rime,cloud}}}{\rho_{qm,\mathrm{cloud}}},
$$

or, in tendency form,

$$
\frac{db_{\mathrm{rime,cloud}}}{dt} =
\frac{dq_{\mathrm{rime,cloud}}/dt}{\rho_{qm,\mathrm{cloud}}}.
$$

Using the implementation names,

- $dq_{\mathrm{rime,cloud}}/dt \rightarrow$ `qc2qi_collect_tend`
- $\rho_{qm,\mathrm{cloud}} \rightarrow$ diagnosed density of collected cloud liquid
- `bm` $\rightarrow$ prognostic rime volume
- `qm` $\rightarrow$ prognostic rime mass

So the cloud-liquid collection contribution to the prognostic rime-volume
tendency is

$$
\left.\frac{d\,bm}{dt}\right|_{\mathrm{cloud\;riming}} =
\frac{qc2qi\_collect\_tend}{\rho\_qm\_cloud}.
$$

In `update_prognostic_ice`, this closure appears as an update of `bm` using the
collected cloud-liquid rime tendency divided by `rho_qm_cloud`. A larger
`rho_qm_cloud` adds less rime volume for the same collected mass and therefore
produces denser rimed ice. A smaller `rho_qm_cloud` adds more volume and
therefore produces lower-density rime.

## Physical Context

Cloud-liquid riming occurs when ice collects supercooled cloud droplets. The
density of the accreted rime depends on droplet size, relative impact speed,
and temperature. Denser rime is expected for:

- larger cloud droplets
- larger ice-droplet impact speed
- temperatures closer to freezing

Lower-density rime is expected for weaker impacts and colder conditions. The
implementation maps these controls through the impact parameter $Ri$, clamps
that parameter to the interval $[1, 12]$, and then maps it to
`rho_qm_cloud`.

## Governing Local Formulas

The active collection mask is

$$
\mathrm{active} =
(qc2qi\_collect\_tend \ge Q_{\mathrm{small}})
\land (T_{\mathrm{atm}} < T_{\mathrm{zerodegc}})
\land \mathrm{context}.
$$

For active lanes,

$$
vtrmi1 = \mathrm{table\_val\_qi\_fallspd} \cdot \mathrm{rhofaci}.
$$

When cloud liquid is also active, $qc_{\mathrm{incld}} \ge Q_{\mathrm{small}}$,
the routine computes

$$
V_{t,qc} =
acn \frac{\Gamma(4+bcn+\mu_c)}{\lambda_c^{bcn} \Gamma(4+\mu_c)}.
$$

Because $bcn = 2$, this simplifies to

$$
V_{t,qc} = acn \frac{(\mu_c+4)(\mu_c+5)}{\lambda_c^2}.
$$

The remaining local quantities are

$$
D_c = \frac{4+\mu_c}{\lambda_c},
\qquad
V_{\mathrm{impact}} = |vtrmi1 - V_{t,qc}|,
$$

$$
\mathrm{inv\_Tc} = \frac{1}{\min(-0.001, T_{\mathrm{atm}} - T_{\mathrm{zerodegc}})},
$$

$$
Ri = \max\left(1,\min\left(-0.5 \times 10^6 D_c V_{\mathrm{impact}}\,
\mathrm{inv\_Tc}, 12\right)\right).
$$

The routine then applies the clamp $1 \le Ri \le 12$ and maps it to rime
density using a piecewise fit.

For $Ri \le 8$,

$$
\rho_{qm,\mathrm{cloud}} = 1000\left(0.051 + 0.114 Ri - 0.0055 Ri^2\right).
$$

For $Ri > 8$,

$$
\rho_{qm,\mathrm{cloud}} = 611 + 72.25(Ri - 8).
$$

When the lane is active but $qc_{\mathrm{incld}} < Q_{\mathrm{small}}$,

$$
\rho_{qm,\mathrm{cloud}} = 400.
$$

When the active mask is false but `context` is true,

$$
vtrmi1 = 0,
\qquad
\rho_{qm,\mathrm{cloud}} = 400.
$$

When `context` is false, the routine does not modify either output lane. The
incoming values of `vtrmi1` and `rho_qm_cloud` are preserved.

The default value is therefore `rho_qm_cloud = 400 kg m^-3`, and the wet-growth
upper endpoint of the fit is `rho_qm_cloud = 900 kg m^-3` at `Ri = 12`.

## Mathematical Properties

The density function is monotone increasing over $1 \le Ri \le 12$. For the
quadratic branch,

$$
\rho_{qm,\mathrm{cloud}}(Ri) = -5.5 Ri^2 + 114 Ri + 51,
$$

so

$$
\frac{d\rho_{qm,\mathrm{cloud}}}{dRi} = 114 - 11 Ri.
$$

On the implemented interval $1 \le Ri \le 8$, this derivative is positive,
ranging from `103` at `Ri = 1` to `26` at `Ri = 8`. The linear branch has
slope

$$
\frac{d\rho_{qm,\mathrm{cloud}}}{dRi} = 72.25.
$$

Thus the density is monotone increasing on both branches. The value is
continuous at `Ri = 8`, but the derivative jumps from `26` to `72.25`, so the
tests should require continuity and monotonicity, not smoothness. The
magnitude of the temperature difference in the denominator is bounded below by
`0.001 K`. For temperatures within `0.001 K` of freezing from below, the
routine uses an effective denominator of `-0.001 K`. The routine is inactive
at exactly `T_zerodegc` because the active mask uses `T_atm < T_zerodegc`.

Useful exact values of the piecewise fit are:

- `Ri = 1` gives `rho_qm_cloud = 159.5 kg m^-3`
- `Ri = 4` gives `rho_qm_cloud = 419.0 kg m^-3`
- `Ri = 8` gives `rho_qm_cloud = 611.0 kg m^-3`
- `Ri = 10` gives `rho_qm_cloud = 755.5 kg m^-3`
- `Ri = 12` gives `rho_qm_cloud = 900.0 kg m^-3`

## Property-Test Strategy

The unit tests live in
`components/eamxx/src/physics/p3/tests/p3_calc_rime_density_unit_tests.cpp`.

The `run_phys()` implementation uses separate Catch2 sections for:

- `activation_and_defaults`
- `velocity_identity`
- `targeted_ri_density_regimes`
- `rime_density_increases_with_impact_speed`
- `rime_density_increases_toward_freezing`
- `impact_speed_symmetry`
- `rime_density_bounds`
- `temperature_floor_prevents_singularity`
- `context_mask_preserves_inactive_lanes`

These tests construct deterministic lanes, run the production
`Functions::calc_rime_density` kernel, and compare against independent helper
formulas or one-sided bounds. The BFB path remains separate in `run_bfb()`.

## Tolerance Philosophy

The tests use identity tolerances for exact algebraic formulas such as the
`vtrmi1` identity and the default-value branches. One-sided inequality
tolerances are used for monotonicity and bound checks. Near-zero floors are
kept separate from identity tolerances so that exact-zero branches are not
confused with scaled relative checks.
