# P3 Autoconversion

This document details the implementation and verification of the cloud water to
rain autoconversion process in the P3 microphysics scheme.

## Physical Formulation

The autoconversion parameterization follows the formulation of **Khairoutdinov
and Kogan (2000)**. It represents the collision and coalescence of cloud
droplets to form rain drops.

### Autoconversion Rate

The time tendency of rain mass mixing ratio due to autoconversion,
$\left(\frac{\partial q_r}{\partial t}\right)_{\text{auto}}$, is calculated as a
power law function of the cloud water mixing ratio ($q_c$) and the cloud droplet
number concentration ($N_c$).

$$
\left(\frac{\partial q_r}{\partial t}\right)_{\text{auto}} = A \, q_c^{\alpha} \, (N_{c, \text{vol}})^{-\beta}
$$

Where:

* $A$ is the autoconversion prefactor ( `autoconversion_prefactor` ).
* $\alpha$ is the cloud water exponent ( `autoconversion_qc_exponent` ).
* $\beta$ is the cloud number exponent ( `autoconversion_nc_exponent` ).
  * **Note**: $\beta$ is stored as a positive value (1.79 in default
     configuration) following KK2000 convention. In the implementation this
     appears as
     `pow(nc_incld * sp(1.e-6) * rho, -autoconversion_nc_exponent)`, which is
     the code equivalent of the analytic form
     $\mathrm{pow}(N_{c, \text{vol}}, -\beta)$.
* $N_{c, \text{vol}}$ is the cloud droplet number concentration in units of
  $\text{cm}^{-3}$.

In the implementation, $N_{c, \text{vol}}$ is derived from the in-cloud number
mixing ratio $N_c$ ($\#/\text{kg}$) and air density $\rho$ ($\text{kg}/\text{m}^3$),
where $N_c$ corresponds to the variable `nc_incld`:

$$
N_{c,\text{vol}} = N_c \cdot \rho \cdot 10^{-6}
\;\equiv\; \texttt{nc\_incld} \cdot \texttt{sp(1.e-6)} \cdot \rho
$$

\equiv \texttt{nc_incld} \cdot \texttt{sp(1.e-6)} \cdot \rho
### Number Tendencies

Consistent with the 2-moment approach, the parameterization calculates the
tendencies for both mass and number.

#### 1. Cloud-to-Rain Number Transfer

The rate at which cloud droplet number converts to rain droplet number is
determined by preserving the specific number concentration:

$$
\left|\frac{\partial N_c}{\partial t}\right|_{\text{auto}} = \left(\frac{\partial q_r}{\partial t}\right)_{\text{auto}} \frac{N_c}{q_c}
$$

**Implementation**: The variable `nc2nr_autoconv_tend` stores the positive
magnitude of this transfer rate. The actual tendency of cloud droplet number is
the negative of this value:

$$
\frac{\partial N_c}{\partial t}
=
-\texttt{nc2nr\_autoconv\_tend}
$$

#### 2. Rain Droplet Number Source

The source of rain droplet number is determined by the characteristic mass of
newly formed rain embryos:

$$
\left(\frac{\partial N_r}{\partial t}\right)_{\text{auto}} = \frac{\left(\frac{\partial q_r}{\partial t}\right)_{\text{auto}}}{m_{\text{drop}}(r_{\text{auto}})}
$$

where $m_{\text{drop}}(r) = \frac{4}{3}\pi\rho_w r^3$ and $r_{\text{auto}}$ is
the characteristic autoconversion radius (default: 25 $\mu\text{m}$).

**Implementation**: The variable `ncautr` stores this rain number source. It is
computed as:

```cpp
ncautr = qc2qr_autoconv_tend * CONS3
```

where `CONS3` is computed locally as:

$$
\text{CONS3} = \frac{1}{\frac{4\pi}{3}\rho_w r_{\text{auto}}^3} = \frac{1}{m_{\text{drop}}(r_{\text{auto}})}
$$

**Important**: Note that `ncautr` (rain number source) is **not equal** to
`nc2nr_autoconv_tend` (cloud number transfer magnitude). The relationship is:
$$
\texttt{ncautr}
=
\texttt{nc2nr\_autoconv\_tend}
\times \frac{q_c}{N_c}
\times \texttt{CONS3}
$$

## Implementation Details

* **Reference**: This logic is implemented in `cloud_water_autoconversion()`
  within [p3_autoconversion_impl.hpp](https://github.com/E3SM-Project/E3SM/blob/master/components/eamxx/src/physics/p3/impl/p3_autoconversion_impl.hpp).
* **Thresholds**: The process is active only when $q_c > 10^{-8} \, \text{kg}/\text{kg}$.
* **Consistency Enforcement**: The implementation ensures that if either the
  mass or number tendency is zero, both are set to zero, preventing numerical
  artifacts in extreme regimes:

  ```cpp
  nc2nr_autoconv_tend.set(qc2qr_autoconv_tend == 0 && context, 0);
  qc2qr_autoconv_tend.set(nc2nr_autoconv_tend == 0 && context, 0);
  ```

* **Subgrid Variance**: The code includes a placeholder for subgrid variance
  scaling (Jensen's inequality effect), but it is currently disabled:

  ```cpp
  // sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, sp(2.47) );
  sgs_var_coef = 1;
  ```

  **Note**: When enabled, the scaling factor uses the same exponent
  ($\alpha = 2.47$) as the $q_c$ dependency, consistent with Jensen's
  inequality for power-law functions:

  $$
  \langle q_c^{\alpha} \rangle > \langle q_c \rangle^{\alpha}
  $$

  when variance exists.

## Property Tests

The unit testing suite (`p3_autoconversion_unit_tests.cpp`) employs a "Physics
Property Test" strategy. We validate that the implementation adheres to
fundamental physical principles across a wide parameter space.

### Tolerance Philosophy

The test suite uses two categories of tolerances based on what is being
validated:

#### 1. Identity Tolerances (Precision-Dependent)

These validate that the implementation follows its documented mathematical
formulas exactly. Since floating-point precision differs between single and
double precision builds, these tolerances are conditionally compiled:

* **Double precision**: `1e-14` (~100× machine epsilon)
* **Single precision**: `1e-7` (~1000× machine epsilon)

**Examples**:

* Verifying `nc2nr_autoconv_tend = qc2qr_autoconv_tend × (Nc/qc)`
* Conservation laws with exact arithmetic

**Rationale**: These are code correctness checks. Failure indicates
implementation bugs, not physics issues. The tolerance accounts for ~10
floating-point operations while remaining strict enough to catch formula errors.

#### 2. Physics Tolerances (Precision-Independent)

These validate physical approximations and configuration parameters. They do not
depend on floating-point precision:

* **Standard**: `1%` (0.01)

**Examples**:

* Embryo radius consistency (configuration parameter)
* KK2000 parameterization accuracy
* Monotonicity detection thresholds

**Rationale**: Physics uncertainties exceed numerical precision. A 1%
parameterization error is 1% regardless of whether calculations use float or
double.

#### 3. Regime Thresholds (Physical Relevance)

These define when values are physically negligible rather than mathematically
zero:

* **Small droplet regime floor**: `1e-15` $\text{kg}/\text{kg}/\text{s}$ (physically irrelevant
  rate for $r < 1 \, \mu\text{m}$)
* **Absolute floor**: `1e-30` (proxy for numerical zero)

**Rationale**: Power-law formulas produce non-zero values even in physically
irrelevant regimes. The small droplet threshold distinguishes when mean droplet
sizes are too small for efficient collision-coalescence ($r < 1 \, \mu\text{m}$,
representing incipient cloud or sub-threshold conditions) from regimes where
autoconversion is physically active.

### Test Strategy

We perform a dense sampling of the phase space:

* **Cloud Water ($q_c$)**: Logarithmic sweep from $5 \times 10^{-9}$
  (sub-threshold regime) to $10^{-2} \, \text{kg}/\text{kg}$ (heavy precipitation).
  * Below $10^{-8}$: autoconversion inactive (incipient cloud)
  * Above $10^{-8}$: active cloud-to-rain transition
* **Cloud Number ($N_c$)**: Logarithmic sweep from $10^{6}$ to $10^{9} \, \#/\text{m}^3$.

### Parameter Validation

Before running property checks, the test suite verifies that the runtime
configuration parameters match the expected Khairoutdinov and Kogan (2000)
constants (e.g., $A=1350$, $\alpha=2.47$, $r_{\text{auto}}=25 \, \mu\text{m}$).

### 1. Threshold Behavior

For $q_c < 10^{-8} \, \text{kg}/\text{kg}$, the test verifies that the autoconversion
rate is strictly zero.

### 2. Monotonicity & Sensitivity

We verify the sign of the partial derivatives with respect to the input state
variables:

* **Colloidal Stability (Inverse $N_c$ Dependency)**:
  Increasing droplet number concentration while holding water content fixed must
  strictly decrease the rate. The test enforces this with **no relative
  tolerance**, as the physics dictates a strict inverse relationship.

  $$
  \frac{\partial}{\partial N_c} \left(\frac{\partial q_r}{\partial t}\right)_{\text{auto}} < 0
  $$

* **Water Content Sensitivity (Positive $q_c$ Dependency)**:
  Increasing water content while holding number concentration fixed should
  increase the rate.

  $$
  \frac{\partial}{\partial q_c} \left(\frac{\partial q_r}{\partial t}\right)_{\text{auto}} > 0
  $$

### 3. Consistency Constraints

We check that the derived number tendencies match their physical definitions.

* **Specific Loss Conservation** (Mathematical Identity):
  The code checks that the number loss rate satisfies the conservation
  relationship: $(dNc/dt)_{\text{auto}} = (dqr/dt)_{\text{auto}} \cdot (N_c / q_c)$.

  **Tolerance**: Precision-dependent (`1e-14` for double, `1e-7` for float).
  This is an exact mathematical identity in the implementation, so failures
  indicate code bugs rather than physics issues.

* **Rain Embryo Mass** (Configuration Consistency):
  The implicit mass of the newly formed rain drops is recovered from the ratio
  of mass tendency to number tendency:

  $$
  m_{\text{effective}} = \frac{\left(\partial q_r / \partial t\right)}{\left(\partial N_r / \partial t\right)}
  $$

  The test verifies that this mass corresponds to the characteristic radius
  $r_{\text{auto}}$ (25 $\mu\text{m}$) within **1% tolerance** (physics-based,
  precision-independent). This confirms the `autoconversion_radius`
  configuration parameter is correctly preserved through the calculation.

### 4. Physical Limits

* **Small Droplet Regime Limit** (Physical Relevance Check):
  When the mean droplet radius is very small ($r < 1 \, \mu\text{m}$, characteristic of
  incipient cloud or sub-threshold conditions), autoconversion must be
  physically negligible due to low collision efficiency. **The test uses a
  threshold of $10^{-15} \, \text{kg}/\text{kg}/\text{s}$** to distinguish "physically irrelevant"
  from "mathematically zero." This accounts for the power-law formula producing
  non-zero values in regimes where collision-coalescence has no physical
  significance. This is a regime threshold, not a precision-dependent tolerance.

### 5. Variance Scaling (Disabled)

The test suite checks for subgrid variance scaling.

* **Principle**: Due to the nonlinearity of the autoconversion rate, subgrid
  variability in $q_c$ should enhance the grid-mean rate.
* **Status**: Since the feature is currently disabled in the implementation
  (`sgs_var_coef = 1`), the test detects this state and reports it as
  "DISABLED" without failing. If the rates differ significantly but do not show
  enhancement, it flags a failure.

### Tolerance Summary

| Test Category | Tolerance Type | Value (Double/Single) | Rationale |
| ------------- | -------------- | --------------------- | --------- |
| **Specific Loss Conservation** | Identity | `1e-14` / `1e-7` | Mathematical identity (2 FLOPs) |
| **Rain Embryo Mass** | Physics | `1%` (both) | Configuration parameter accuracy |
| **Monotonicity (Nc)** | Strict | No tolerance | Physics requires strict inequality |
| **Monotonicity (qc)** | Physics | `0.1%` | Physical sensitivity detection |
| **Small Droplet Regime** | Physical | `1e-15` $\text{kg}/\text{kg}/\text{s}$ | Relevance threshold ($r < 1 \, \mu\text{m}$) |
| **Variance Detection** | Feature | `0.1%` | Distinguish on/off state |
