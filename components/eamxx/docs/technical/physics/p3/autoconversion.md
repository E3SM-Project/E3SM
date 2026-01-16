# P3 Autoconversion

This document details the implementation and verification of the cloud water to rain autoconversion process in the P3 microphysics scheme.

## Physical Formulation

The autoconversion parameterization follows the formulation of **Khairoutdinov and Kogan (2000)**. It represents the collision and coalescence of cloud droplets to form rain drops.

### Autoconversion Rate

The time tendency of rain mass mixing ratio due to autoconversion, $ \left(\frac{\partial q_r}{\partial t}\right)_{auto} $, is calculated as a power law function of the cloud water mixing ratio ($q_c$) and the cloud droplet number concentration ($N_c$).

$$ \left(\frac{\partial q_r}{\partial t}\right)_{auto} = A \, q_c^{\alpha} \, (N_{c, vol})^{-\beta} $$

Where:

*   $ A $ is the autoconversion prefactor ( `autoconversion_prefactor` ).
*   $ \alpha $ is the cloud water exponent ( `autoconversion_qc_exponent` ).
*   $ \beta $ is the cloud number exponent ( `autoconversion_nc_exponent` ).
    *   **Note**: $\beta$ is stored as a positive value (1.79 in default configuration) following KK2000 convention. The implementation uses `pow(Nc_vol, -1.79)`.
*   $ N_{c, vol} $ is the cloud droplet number concentration in units of $ cm^{-3} $.

In the implementation, $ N_{c, vol} $ is derived from the in-cloud number mixing ratio $ N_c $ ($ \#/kg $) and air density $ \rho $ ($ kg/m^3 $):

$$ N_{c, vol} = N_c \cdot \rho \cdot 10^{-6} $$

### Number Tendencies

Consistent with the 2-moment approach, the parameterization calculates the tendencies for both mass and number.

**1. Cloud Droplet Number Loss**
The loss of cloud droplet number ($ \partial N_c / \partial t $) is assumed to be proportional to the loss of cloud mass. The specific loss rate is conserved:

$$ \frac{1}{q_c} \left(\frac{\partial q_c}{\partial t}\right)_{auto} = \frac{1}{N_c} \left(\frac{\partial N_c}{\partial t}\right)_{auto} $$

Thus:

$$ \left(\frac{\partial N_c}{\partial t}\right)_{auto} = \left(\frac{\partial q_r}{\partial t}\right)_{auto} \frac{N_c}{q_c} $$

Note: The implementation variable `nc2nr_autoconv_tend` typically stores the magnitude of this loss term.

**2. Rain Number Source**
The source of rain number ($ \partial N_r / \partial t $) is determined by assuming that all new rain drops formed by autoconversion have a characteristic radius, $ r_{auto} $ (typically $ 25 \mu m $, configurable via `autoconversion_radius`).

$$ \left(\frac{\partial N_r}{\partial t}\right)_{auto} = \frac{\left(\frac{\partial q_r}{\partial t}\right)_{auto}}{m_{drop}(r_{auto})} $$

Where $ m_{drop}(r) = \frac{4}{3} \pi \rho_w r^3 $. In the code, this mass is related to `CONS2`: $ m_{drop}(r) = \text{CONS2} \times r^3 $, where $ \text{CONS2} = \frac{4\pi}{3}\rho_w $.

## implementation Details

*   **Reference**: This logic is implemented in `cloud_water_autoconversion()` within [p3_autoconversion_impl.hpp](../../../../src/physics/p3/impl/p3_autoconversion_impl.hpp).
*   **Thresholds**: The process is active only when $ q_c > 10^{-8} \, kg/kg $.
*   **Subgrid Variance**: The code includes a placeholder for subgrid variance scaling (Jensen's inequality effect), but it is currently disabled:
    ```cpp
    // sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, sp(2.47) );
    sgs_var_coef = 1;
    ```
    This means the rate is calculated using grid-mean values without enhancement factors.

## Property Tests

The unit testing suite (`p3_autoconversion_unit_tests.cpp`) employs a "Physics Property Test" strategy. We validate that the implementation adheres to fundamental physical principles across a wide parameter space.

### Test Strategy
We perform a dense sampling of the phase space:
*   **Cloud Water ($q_c$)**: Logarithmic sweep from $ 5 \times 10^{-9} $ (below threshold) to $ 10^{-2} \, kg/kg $.
*   **Cloud Number ($N_c$)**: Logarithmic sweep from $ 10^{6} $ to $ 10^{9} \, \#/m^3 $.

### Parameter Validation
Before running property checks, the test suite verifies that the runtime configuration parameters match the expected Khairoutdinov and Kogan (2000) constants (e.g., $A=1350$, $\alpha=2.47$, $r_{auto}=25\mu m$).

### 1. Threshold Behavior
For $q_c < 10^{-8} \, kg/kg$, the test verifies that the autoconversion rate is strictly zero.

### 2. Monotonicity & Sensitivity

We verify the sign of the partial derivatives with respect to the input state variables (using a relative tolerance check):

*   **Colloidal Stability (Inverse $N_c$ Dependency)**:
    Increasing droplet number concentration while holding water content fixed should decrease the rate.
    $$ \frac{\partial}{\partial N_c} \left(\frac{\partial q_r}{\partial t}\right)_{auto} < 0 $$

*   **Water Content Sensitivity (Positive $q_c$ Dependency)**:
    Increasing water content while holding number concentration fixed should increase the rate.
    $$ \frac{\partial}{\partial q_c} \left(\frac{\partial q_r}{\partial t}\right)_{auto} > 0 $$

### 3. Consistency Constraints

We check that the derived number tendencies match their physical definitions.

*   **Specific Loss Conservation**:
    The code checks that the number loss rate satisfies the conservation relationship: $ (dNc/dt)_{auto} = (dqr/dt)_{auto} \cdot (N_c / q_c) $.
    
*   **Rain Embryo Mass**:
    The implicit mass of the newly formed rain drops is recovered from the ratio of mass tendency to number tendency:
    $$ m_{effective} = \frac{\left(\partial q_r / \partial t\right)}{\left(\partial N_r / \partial t\right)} $$
    The test verifies that this mass corresponds to a drop size of roughly $20-30 \mu m$ (consistent with $r_{auto} = 25 \mu m$).

### 4. Physical Limits

*   **Haze Limit**:
    If the mean droplet radius corresponds to haze ($ r < 1 \mu m $), the autoconversion rate must be negligible.

### 5. Variance Scaling (Disabled)

The test suite checks for subgrid variance scaling.
*   **Principle**: Due to the nonlinearity of the autoconversion rate, subgrid variability in $ q_c $ should enhance the grid-mean rate.
*   **Status**: Since the feature is currently disabled in the implementation (`sgs_var_coef = 1`), the test detects this state and reports it as "DISABLED" without failing. If the rates differ significantly but do not show enhancement, it flags a failure.
