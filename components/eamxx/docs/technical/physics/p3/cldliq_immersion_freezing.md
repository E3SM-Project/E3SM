# P3 Cloud-Liquid Immersion Freezing

This document describes the P3 routine `cldliq_immersion_freezing`. The
routine computes heterogeneous immersion freezing tendencies for cloud-liquid
mass and cloud-droplet number, and it summarizes the physical-property tests
that verify the implementation.

## Purpose

The routine computes:

- `qc2qi_hetero_freeze_tend`: cloud-liquid mass converted to ice by
  heterogeneous immersion freezing
- `nc2ni_immers_freeze_tend`: cloud-droplet number converted to ice number by
  immersion freezing

These tendencies are used when supercooled cloud droplets freeze into ice in
the P3 microphysics scheme.

## Physical Context

Cloud-liquid immersion freezing represents heterogeneous freezing of
supercooled cloud droplets. For the current positive runtime exponent, the
rate increases with supercooling and depends on the assumed cloud-droplet size
distribution. The mass tendency and number tendency are different moments of
that same size distribution, so their ratio encodes the mean mass of newly
frozen droplets.

## Active Mask

The routine is active when all of the following conditions hold:

- `qc_incld >= QSMALL`
- `T_atm <= T_rainfrz`
- `context == true`

Both thresholds are inclusive. If any of these conditions fail, the routine
does not modify the output lane.

## Governing Formulas

For active lanes, define

$$
\theta = T_{\mathrm{zerodegc}} - T_{\mathrm{atm}},
\qquad
a = \texttt{runtime\_options.immersion\_freezing\_exponent},
$$

$$
\mathrm{expAimmDt} = e^{a\theta},
\qquad
\lambda = \lambda_c,
\qquad
\mu = \mu_c,
\qquad
c = \texttt{cdist1}.
$$

The implementation computes

$$
qc2qi_{\mathrm{hetero\_freeze\_tend}} =
CONS6 \; c \; \Gamma(\mu + 7) \; e^{a\theta} \; \lambda^{-6},
$$

$$
nc2ni_{\mathrm{immers\_freeze\_tend}} =
CONS5 \; c \; \Gamma(\mu + 4) \; e^{a\theta} \; \lambda^{-3}.
$$

The code forms these powers through `inv_lamc3 = (1/lamc)^3`, so the mass
tendency uses `square(inv_lamc3)` and the number tendency uses `inv_lamc3`.

## Moment Identity

Let

$$
M = qc2qi_{\mathrm{hetero\_freeze\_tend}},
\qquad
N = nc2ni_{\mathrm{immers\_freeze\_tend}}.
$$

Then the mass-to-number ratio simplifies exactly to

$$
\frac{M}{N} =
\frac{\mathrm{CONS6}}{\mathrm{CONS5}}
\frac{(\mu_c + 4)(\mu_c + 5)(\mu_c + 6)}{\lambda_c^3}.
$$

This follows from

$$
\frac{\Gamma(\mu + 7)}{\Gamma(\mu + 4)}
= (\mu + 4)(\mu + 5)(\mu + 6).
$$

Using

$$
\mathrm{CONS5} = \frac{\pi}{6}\mathrm{BIMM},
\qquad
\mathrm{CONS6} = \left(\frac{\pi}{6}\right)^2 \rho_w\mathrm{BIMM},
$$

the constant ratio is

$$
\frac{\mathrm{CONS6}}{\mathrm{CONS5}} = \frac{\pi \rho_w}{6}.
$$

This ratio is the mean mass per newly frozen droplet implied by the assumed
cloud-droplet size distribution. Since
\(\mathrm{CONS6} / \mathrm{CONS5} = \pi \rho_w / 6\), the ratio has the same
structure as the mass of a spherical water droplet, with the effective
droplet-volume factor supplied by the Gamma-distribution moment ratio.
Temperature, the runtime exponent, and the distribution prefactor `cdist1`
cancel out because they multiply the mass and number moments identically.

## Volume-Weighted Freezing Interpretation

The implemented mass-to-number tendency ratio can be interpreted as the mean
mass per newly frozen droplet. Because

$$
\frac{CONS6}{CONS5} = \frac{\pi \rho_w}{6},
$$

the ratio has the structure of the mass of a spherical water droplet,
\(m = (\pi \rho_w / 6)D^3\).

For an ordinary Gamma cloud droplet number distribution, the ordinary mean
cloud droplet mass is

$$
m_{\mathrm{cloud}} =
\frac{\pi \rho_w}{6}
\frac{(\mu_c + 1)(\mu_c + 2)(\mu_c + 3)}{\lambda_c^3}.
$$

For the newly frozen droplets implied by the P3 immersion-freezing closure,

$$
m_{\mathrm{frozen}} =
\frac{\pi \rho_w}{6}
\frac{(\mu_c + 4)(\mu_c + 5)(\mu_c + 6)}{\lambda_c^3}.
$$

Therefore,

$$
\frac{m_{\mathrm{frozen}}}{m_{\mathrm{cloud}}}
= \frac{(\mu_c + 4)(\mu_c + 5)(\mu_c + 6)}
  {(\mu_c + 1)(\mu_c + 2)(\mu_c + 3)}
> 1
\qquad \text{for } \mu_c \ge 0.
$$

This expresses a general microphysical property of volume-proportional
immersion freezing: newly frozen droplets are biased toward larger droplets
than the ordinary cloud droplet population, because larger droplets have more
volume in which to contain an immersed ice-nucleating site.

## Implementation Notes

`qc_incld` is currently a gate only. Above threshold it does not enter the
active-lane formulas directly. Cloud-water dependence enters through the
diagnosed distribution parameters `lamc`, `mu_c`, and `cdist1`.
This does not mean cloud water is physically irrelevant; it means this small
kernel receives cloud-water dependence through the precomputed droplet-spectrum
parameters.

`inv_qc_relvar` is currently inactive in this routine because the code uses

```cpp
// sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, 2);
sgs_var_coef = 1;
```

If subgrid-variance scaling is restored, the physical-property tests should be
updated from invariance to the corresponding exact scaling law.

When `context` is false, or when the lane fails the `qc_incld` or temperature
activation gate, the routine leaves the incoming output values unchanged.
The property tests therefore seed outputs with sentinel values when checking
inactive lanes.

## Mathematical Properties

For active lanes, the separable structure implies exact logarithmic scaling
identities:

$$
\frac{d\ln M}{d\theta} = a,
\qquad
\frac{d\ln N}{d\theta} = a,
$$

$$
\frac{d\ln M}{d\ln \lambda_c} = -6,
\qquad
\frac{d\ln N}{d\ln \lambda_c} = -3,
$$

$$
\frac{d\ln M}{d\ln cdist1} = 1,
\qquad
\frac{d\ln N}{d\ln cdist1} = 1.
$$

These identities encode the key physics:

- for positive `a`, stronger supercooling increases both tendencies
  exponentially
- smaller `lamc` corresponds to larger characteristic droplets, with a stronger
  effect on frozen mass than on frozen number because mass samples a higher
  moment of the droplet distribution
- `cdist1` scales both tendencies linearly
- `cdist1 = 0` gives zero mass and number freezing tendencies for otherwise
  active lanes

## Physical-Property Tests

The unit tests live in
`components/eamxx/src/physics/p3/tests/p3_cldliq_imm_freezing_unit_tests.cpp`.

The `run_phys()` implementation uses separate Catch2 sections for:

- `activation_and_thresholds`
- `temperature_supercooling_scaling`
- `runtime_exponent_sensitivity`
- `lambda_power_law_scaling`
- `mass_number_moment_identity`
- `volume_nucleation_prefers_larger_than_mean_droplets`
- `distribution_prefactor_scaling`
- `zero_distribution_prefactor_gives_zero`
- `scalar_multiplier_cancellation_in_mass_number_ratio`
- `qc_gate_only_behavior`
- `inv_qc_relvar_currently_inactive`
- `context_mask_preserves_inactive_lanes`
- `finite_nonnegative_outputs`

These tests use deterministic inputs, launch the production
`Functions::cldliq_immersion_freezing` kernel through Kokkos, copy the results
back to host, and compare the outputs against independent identities, exact
ratios, and one-sided physical inequalities. The BFB regression path remains
separate in `run_bfb()`.

## Tolerance Philosophy

The tests use identity tolerances for exact algebraic relations such as the
moment identity and the separable scaling ratios. Exact equality is used for
inactive-lane preservation checks because the kernel should not touch those
lanes. Near-zero floors remain separate from scaled relative checks so that
inactive preservation and positive-tendency assertions are not conflated. The
volume-weighted size-bias test also uses this identity tolerance, because the
comparison still goes through production-kernel outputs involving `tgamma`,
`exp`, and ratio formation.
