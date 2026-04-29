# P3 Liquid Relaxation Timescale

This document describes the P3 routine
`calc_liq_relaxation_timescale`. Despite the routine name, the outputs `epsr`
and `epsc` are relaxation-rate coefficients used in the liquid adjustment
calculation. This document also describes the physical-property tests that
verify the implementation.

## Purpose

The routine computes two liquid-phase relaxation-rate coefficients used by P3:

- `epsr`, the rain relaxation-rate coefficient
- `epsc`, the cloud-liquid relaxation-rate coefficient

The rain expression combines an analytic gamma-distribution term with a
ventilation-enhanced lookup-table term. The cloud-liquid expression is a direct
closed-form product.

## Physical Formulation

For active rain lanes, defined by

$$
q_{r,\mathrm{incld}} \ge Q_{\mathrm{small}}
\quad\text{and}\quad
\mathrm{context}=\mathrm{true},
$$

define

$$
T(\mu_r,\lambda_r)
=
\mathrm{apply\_table}(\mathrm{revap\_table\_vals},
\mathrm{lookup}(\mu_r,\lambda_r)),
$$

$$
A = \frac{\Gamma(\mu_r+2)}{\lambda_r},
\qquad
B = \left(\frac{\rho}{\mu}\right)^{1/2} Sc^{1/3} T(\mu_r,\lambda_r).
$$

Then the rain relaxation-rate coefficient is

$$
\epsilon_r =
2\pi C_r \rho D_v
\left[
f_{1r} A + f_{2r} B
\right].
$$

For active cloud-liquid lanes, defined by

$$
q_{c,\mathrm{incld}} \ge Q_{\mathrm{small}}
\quad\text{and}\quad
\mathrm{context}=\mathrm{true},
$$

the cloud-liquid relaxation-rate coefficient is

$$
\epsilon_c = 2\pi \rho D_v C_c.
$$

For cloud-small lanes inside context,

$$
q_{c,\mathrm{incld}} < Q_{\mathrm{small}},
$$

the implementation sets

$$
\epsilon_c = 0.
$$

For rain-small lanes, and for lanes outside context, the implementation returns
`epsr = 0` because it initializes the entire `epsr` pack to zero before
overwriting active rain lanes.

## Variable Definitions

- $\rho$: air density
- $D_v$: vapor diffusivity, `dv`
- $\mu$: dynamic viscosity used in the ventilation term
- $Sc$: Schmidt number, `sc`
- $\mu_r$: rain size-distribution shape parameter
- $\lambda_r$: rain size-distribution slope parameter, `lamr`
- $C_r$: rain prefactor, `cdistr`
- $C_c$: cloud-liquid prefactor, `cdist`
- $f_{1r}$: analytic rain coefficient
- $f_{2r}$: ventilation rain coefficient
- $Q_{\mathrm{small}}$: hydrometeor activity threshold used separately for
  `qr_incld` and `qc_incld`
- $T(\mu_r,\lambda_r)$: rain evaporation lookup-table value

## Implementation Details

- Source:
  `components/eamxx/src/physics/p3/impl/p3_calc_liq_relaxation_timescale_impl.hpp`
- The threshold comparison is inclusive:
  - `qr_incld >= QSMALL`
  - `qc_incld >= QSMALL`
- `epsr` is initialized globally before the context-filtered rain overwrite:

  ```cpp
  epsr = 0;
  ```

  This means out-of-context lanes are guaranteed to return `epsr = 0`.

- `epsc` is only modified through `.set(... && context, ...)` operations.
  Therefore, when `context == false`, the routine does not guarantee that
  `epsc` is overwritten.

## Lookup-Table Role

The lookup-table contribution enters only the rain expression through
$T(\mu_r,\lambda_r)$. This term depends on both `mu_r` and `lamr`, so tests
that vary those inputs in the full rain expression can mix algebraic scaling
with interpolation effects.

To avoid false positives, the property tests isolate exact algebraic identities
by holding `mu_r` and `lamr` fixed whenever the lookup-table value must remain
unchanged.

## Mathematical Properties

For active rain lanes, the implementation can be decomposed into an analytic
gamma-distribution term and a ventilation-enhanced table term:

$$
\epsilon_r =
2\pi C_r \rho D_v
\left[
f_{1r}\frac{\Gamma(\mu_r+2)}{\lambda_r}
+
f_{2r}\left(\frac{\rho}{\mu}\right)^{1/2}
Sc^{1/3} T(\mu_r,\lambda_r)
\right].
$$

This decomposition implies exact additivity between the analytic and
ventilation terms. It also implies exact linear scaling with $D_v$ and $C_r$.
When the ventilation term is disabled, the expression is linear in $\rho$,
$f_{1r}$, and $\lambda_r^{-1}$. When the analytic term is disabled and the
lookup-table term is held fixed, the expression scales as $\rho^{3/2}$,
$\mu^{-1/2}$, and $Sc^{1/3}$.

For active cloud-liquid lanes,

$$
\epsilon_c = 2\pi \rho D_v C_c,
$$

so the cloud-liquid coefficient is exactly linear in $\rho$, $D_v$, and
$C_c$.

The mixed rain expression also yields a robust density bound when both terms
are active and non-negative:

$$
2\epsilon_r(\rho)
\le
\epsilon_r(2\rho)
\le
2\sqrt{2}\,\epsilon_r(\rho).
$$

### Connection to Bulk Microphysics

The cloud-liquid expression has the form of a diffusion-limited vapor exchange
coefficient. It is linear in vapor diffusivity, air density, and the cloud
geometric distribution factor. The rain expression adds a ventilation-enhanced
term. The factor $(\rho/\mu)^{1/2}Sc^{1/3}$ represents the expected
Reynolds/Schmidt-number dependence of ventilation corrections in bulk
microphysics.

This structure implies that increasing diffusivity, air density, or the
distribution prefactors increases the relaxation-rate coefficients. Increasing
viscosity damps the ventilation contribution, while increasing Schmidt number
enhances it. The tests encode these expectations as exact scaling identities
when the lookup-table term is held fixed.

## Property Tests

The unit tests live in
`components/eamxx/src/physics/p3/tests/p3_calc_liq_relaxation_timescale_unit_tests.cpp`.

### Test Organization

The Catch2 test case `p3_calc_liq_relaxation_timescale` has two top-level
sections:

- `properties`
- `bfb`

Inside `properties`, the nested sections cover:

- context and threshold masking
- cloud-liquid identity and exact scaling
- rain additivity between the analytic and ventilation terms
- full rain prefactor scaling with `dv` and `cdistr`
- analytic-only rain identities
- ventilation-only rain power-law identities
- mixed microphysical scaling with `rho`, `mu`, and `sc`
- mixed-term density bound
- separation of cloud-liquid and rain-only factors
- zero-coefficient and coefficient-locality behavior
- non-negativity and finite-output checks

### Tolerance Philosophy

| Test Category | Tolerance Type | Value | Rationale |
| --- | --- | --- | --- |
| Exact identities and scaling | Identity | `100 * epsilon` | Tight algebraic checks with floating-point slack |
| Bounds and non-negativity | Inequality | `100 * epsilon`-scaled | Stable comparisons for one-sided checks |
| Expected-zero branches | Absolute | `1e-30` | Numerical zero floor for thresholded and zero-coefficient cases |

Here, `epsilon` is machine epsilon for EAMxx `Real`
(`std::numeric_limits<Real>::epsilon()`).

### Context-Mask Coverage

The suite explicitly documents the asymmetric output behavior:

- `epsr` is zero for out-of-context lanes because the routine zeroes the full
  pack before applying the rain mask.
- `epsc` is not guaranteed to be written outside context, so the tests
  initialize it with a sentinel value and verify that the sentinel is
  preserved.

### Table-Sensitive Identities

The suite avoids asserting lookup-sensitive identities that would vary with the
interpolated table value. In particular:

- all non-coefficient inputs are held fixed for the rain-additivity check
- `mu_r` and `lamr` are held fixed for ventilation-only scaling checks
- the exact `lamr^{-1}` scaling is tested only in the analytic-only regime

### BFB Coverage

The BFB regression path is preserved to ensure the property tests do not change
production behavior. It uses representative randomized inputs spanning active
and inactive threshold states and compares the serialized outputs against the
established baseline format.
