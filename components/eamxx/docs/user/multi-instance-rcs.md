# Multi-instance and NBFB testing capabilities

For NBFB testing, we utilize a "System Test" called `RCS`
which stands for Reproducible Climate Statistics. The premise
of RCS is comparing the statistics of two (potentially differing)
populations to determine if they are statistically identical.
To enable RCS, we utilize CIME's multi-instance capability.

## Multi-instance

A given case can be made to house a number of instances
by modifying the `NINST` variable. For a given component `CMP`,
active or not, `NINST_CMP` can be set to a value greater than 1.
After CIME's setup operation (`case.setup`), `CMP` will have
`NINST_CMP` instances of itself available, with their own runtime
options that can be changed.

For the `scream` component, the support is incomplete because
of the departure from the convention of `user_nl_scream` to
utilize a more readable YAML-based configuration. As a result,
the configuration of `NINST_ATM` when `ATM` is `scream`, the user
must replicate `NINST_ATM` copies of `data/scream_input.yaml` with
`data/scream_input.yaml_0001` where the last digits represent the
number of instances (with added leading zeros). With the new input
files, the users can change each individual file's content so that
ensemble member `_0001` reflects a user-specified configuration.

Beyond `NINST`, the user can also choose to utilize the multi-driver
capability. If `MULTI_DRIVER` is `FALSE` (default), then all instances
will be launched from the same coupler instance, which can be problematic
and can result in out-of-memory errors. Instead, the user can set
`MULTI_DRIVER` to `TRUE` which will result in a coupler instance
for each `NINST`.

## RCS testing

The Reproducible Climate Statistics (RCS) test performs a rigorous
statistical comparison between two ensembles of one-year runs to
verify climate reproducibility in not-bit-for-bit (NBFB) cases
using advanced statistical methods.

### Available Statistical Tests

RCS provides a comprehensive suite of 11 two-sample statistical tests organized
into three categories; all tests are implemented in SciPy's stats module
([SciPy stats](https://docs.scipy.org/doc/scipy/reference/stats.html)).

#### Distribution Tests (Compare Entire Distributions)

These tests assess whether two samples come from the same probability distribution:

##### `ks` - Kolmogorov-Smirnov (default)

- Sensitivity: Moderate
- Best for: General-purpose distribution comparison
- Assumption: Continuous distributions
- **Use this for routine testing**

##### `ad` - Anderson-Darling

- Sensitivity: High (especially in distribution tails)
- Best for: Detecting subtle distributional differences
- Assumption: Continuous distributions
- Note: Uses stricter alpha (0.001) due to high sensitivity
- **Use this for sensitivity analysis**

##### `cvm` - Cramér-von Mises

- Sensitivity: Moderate to High
- Best for: Overall distribution comparison with tail sensitivity
- Assumption: Continuous distributions

##### `epps` - Epps-Singleton

- Sensitivity: Moderate
- Best for: Detecting differences in both mean AND variance
- Assumption: Works well for non-normal distributions

##### `energy` - Energy Distance

- Sensitivity: High
- Best for: Detecting any type of distributional difference
- Assumption: None (distribution-free)
- Note: Computationally intensive, uses permutation testing

#### Location Tests (Compare Means/Medians)

These tests focus on differences in central tendency:

##### `mw` - Mann-Whitney U Test

- Sensitivity: Moderate
- Best for: Comparing medians of non-normal distributions
- Assumption: None (non-parametric)
- **Use when distributions may be non-normal**

##### `ttest` - Welch's t-test

- Sensitivity: Moderate to High
- Best for: Comparing means of approximately normal distributions
- Assumption: Approximately normal distributions (robust to violations)

##### `brunner` - Brunner-Munzel

- Sensitivity: Moderate to High
- Best for: Robust alternative to t-test for ordinal data
- Assumption: None (non-parametric)

#### Scale Tests (Compare Variances/Spread)

These tests focus on differences in variability:

##### `levene` - Levene's Test

- Sensitivity: Moderate
- Best for: Testing equality of variances
- Assumption: Robust to non-normality

##### `ansari` - Ansari-Bradley

- Sensitivity: Moderate
- Best for: Non-parametric scale comparison
- Assumption: Samples differ primarily in scale, not location

##### `mood` - Mood's Test

- Sensitivity: Moderate
- Best for: Non-parametric dispersion comparison
- Assumption: None (distribution-free)

### Statistical Test Methodology

RCS supports two complementary analysis modes for comprehensive climate validation:

#### Spatiotemporal Analysis (Default, Recommended)

This mode computes area-weighted global spatial means at each timestep, then
performs a single statistical test per variable comparing the distributions of
these global means across all instances and timesteps.

**Procedure:**

1. For each instance in each ensemble:
   - Variables with vertical dimensions (lev/ilev) are averaged vertically
   - Compute spatial mean at each timestep using grid cell areas as weights
2. Concatenate all timesteps from all instances
3. Remove NaN values from both distributions
4. Perform selected statistical test comparing the two distributions
5. Variable fails if p-value < α (default: 0.01, or 0.001 for Anderson-Darling)

**Characteristics:**

- Sensitive to global systematic biases
- Single test per variable (more conservative, reduces multiple testing issues)
- Requires area weights for proper spatial averaging
- Handles NaN values from land/ocean masks gracefully
- **Recommended for most use cases**

#### Temporal Analysis (Alternative)

This mode computes temporal means at each spatial column, then performs
per-column statistical tests comparing the distributions across instances.

**Procedure:**

1. For each instance in each ensemble:
   - Variables with vertical dimensions are averaged vertically
   - Compute temporal mean at each spatial column
2. For each column, perform selected statistical test comparing distributions
3. Apply multiple testing correction (Bonferroni or FDR) across all columns
4. Variable fails if more than CRITICAL_FRACTION (default: 0.1%) of columns reject

**Characteristics:**

- Detects spatially-localized differences
- More tests per variable (requires multiple testing correction)
- Can identify regional changes that global averages might miss
- **Use for detecting localized numerical changes**

### Configuration Parameters

All configuration parameters can be adjusted via command-line arguments when
running `rcs_stats.py` standalone, or programmatically when calling
`run_stats_comparison()`. Parameters control test sensitivity, multiple testing
corrections, and failure thresholds.

#### Core Statistical Parameters

##### `--alpha` (default: 0.01 for most tests, 0.001 for Anderson-Darling)

- Significance level for hypothesis tests
- Lower values make tests stricter (fewer false positives)
- Higher values increase power (fewer false negatives)
- Example: `--alpha 0.001` for very strict testing
  
##### `--test_type` (default: `ks`)

- Statistical test identifier
- Distribution tests: `ks`, `ad`, `cvm`, `epps`, `energy`
- Location tests: `mw`, `ttest`, `brunner`
- Scale tests: `levene`, `ansari`, `mood`
- Example: `--test_type ad` for high sensitivity

##### `--analysis_type` (default: `spatiotemporal`)

- Analysis mode for data aggregation
- `spatiotemporal`: Area-weighted global means (recommended)
- `temporal`: Per-column temporal means
- Example: `--analysis_type temporal` for localized detection

#### Multiple Testing Correction

When testing hundreds of variables simultaneously, the chance of false positives
increases. Multiple testing corrections adjust significance thresholds to control
error rates across all tests.

##### `--correction_method` (default: `bonferroni`)

Method for correcting multiple comparisons.
Example: `--correction_method bonferroni`

###### `bonferroni`: Conservative, controls family-wise error rate (FWER)

- Divides alpha by number of tests
- Use for: Strict validation, production BFB testing
- Guarantees overall false positive rate ≤ alpha

###### `fdr`: False Discovery Rate (Benjamini-Hochberg procedure)

- Controls expected proportion of false discoveries
- Less conservative than Bonferroni, better power
- Use for: Exploratory analysis, detecting subtle changes
- Allows more true positives while controlling false positives

###### `none`: No correction applied

- Use for: Single variable analysis, pre-screened tests
- Caution: High false positive rate with many tests

#### Failure Thresholds

##### `--critical_fraction` (default: 0.001)

- Maximum fraction of failed sub-tests per variable
- Only used in temporal analysis (per-column tests)
- Variable fails if more than this fraction of columns reject null hypothesis
- Range: 0.0 to 1.0
- Example: `--critical_fraction 0.01` allows 1% of columns to fail
  
##### `--max_failed_vars` (default: 0)

- Maximum number of variables allowed to fail before overall test fails
- Overall test status = FAIL if failed_vars > max_failed_vars
- Use 0 for strict BFB testing (no failures allowed)
- Use higher values for regression testing or exploratory work
- Example: `--max_failed_vars 5` allows up to 5 variable failures

#### Effect Size Filtering

##### `--magnitude_threshold` (default: None)

- Minimum relative difference to consider significant
- Requires BOTH statistical significance (p < alpha) AND practical
  significance (relative difference > threshold)
- Computed as: |mean1 - mean2| / ((|mean1| + |mean2|) / 2)
- Range: 0.0 to 1.0 (e.g., 0.01 = 1% difference)
- Use to filter out statistically significant but tiny differences
- Example: `--magnitude_threshold 0.01` requires >1% mean difference

### Variable Selection

The test automatically identifies suitable variables based on:

- Must have `time` dimension
- Must have spatial dimensions (`ncol`)
- Must not be coordinate variables (time, lat, lon, lev, etc.)
- Must not be entirely NaN or constant-valued

This ensures all physically meaningful prognostic and diagnostic variables
are tested without manual specification.

### NaN Handling

The implementation robustly handles missing values (NaN) throughout:

- All mean calculations use `nanmean`/`nansum`
- NaN values are filtered before K-S tests
- Columns/variables with insufficient valid data are skipped with warnings
- Spatial means properly account for partial land/ocean masks

### Output

The test produces comprehensive diagnostic information:

#### Console output: Summary with configuration details and test results

- Test type, analysis mode, and all configuration parameters
- Alpha level, correction method, and failure thresholds
- Number of instances and variables tested
- Summary counts of passed/failed variables

It also offers detailed statistics for failed variables including:

- Sample sizes and descriptive statistics (mean, std, median, quartiles)
- Mean differences (absolute and percentage)
- Standard deviation ratios
- Human-readable reasons for pass/fail decisions
- Correction method effects (if applicable)

#### Test log: Detailed comments appended to CIME test status

#### JSON file (`{test_type}_test_results.json`): Complete structured results

##### configuration: All parameters used for the test

- `alpha`: Significance level
- `correction_method`: Multiple testing correction method
- `critical_fraction`: Failure threshold for sub-tests
- `max_failed_vars`: Maximum allowed variable failures
- `magnitude_threshold`: Effect size threshold (if set)

##### summary: Overall test results

- `passed`: Number of variables that passed
- `failed`: Number of variables that failed
- `total`: Total number of variables tested
- `test_status`: Overall PASS/FAIL status

##### details: Per-variable comprehensive statistics

- Test statistic and p-value
- Sample1/Sample2 descriptive statistics (n, mean, median, std, min, max, quartiles)
- Differences (mean_diff, mean_diff_pct, median_diff, std_ratio)
- Hypothesis result (PASS/FAIL) with explanation
- Correction metadata (corrected_alpha, fdr_critical_value, correction_method)

##### failed_variables: List of variable names that failed

##### passed_variables: List of variable names that passed

### Running RCS Tests

#### Within CIME System Test Framework

In order to run RCS as a CIME system test, you must request multiple instances.
This can be achieved by adjusting runtime settings as discussed above.
Or if using the `RCS` "System Test", you can simply append the name
of the test with `_N#` (same driver) or `_C#` (multiple drivers).

The user must also enable a perturbation across instances.
A simple addition to `scream` configuration to perturb the initial
condition file (e.g., initial_conditions::perturbed_fields="T_mid")
should suffice. The RCS test will then ensure each instance has a
different seed, and thus follow a different trajectory.
RCS is designed such that it returns identical seeds and thus identical results,
unless code or configuration changes introduce numerical or climate differences.

For convenience, there exists a "testmod" that can enable the perturbation for
the user and can set a monthly average output stream that the RCS test will
copy across instances. With CIME's create_test, the following is recommended:

```shell
./cime/scripts/create_test RCS_P4_C4.$RES.$COMPSET.$MACH.eamxx-perturb
```

where `RCS_P4_C4` will result in 4 multi-driver instances all using a pelayout
of 4, and will use the eamxx-perturb testmod as a helper in the setup phase.
The rest of the options ($RES, $COMPSET, $MACH) should be familiar to users.

#### Standalone Command-Line Usage

The statistical comparison can also be run independently from the command line
for custom analysis or debugging. This is useful for:

- Testing different statistical methods on existing data
- Adjusting significance thresholds without re-running simulations
- Analyzing archived test results
- Developing and validating new test configurations

**Basic Usage:**

```text
# Default: Kolmogorov-Smirnov test with spatiotemporal analysis
rcs_stats.py /run/dir /base/dir
```

```text
# Specify different statistical test
rcs_stats.py /run/dir /base/dir \
    --test_type ad
```

```text
# Use temporal analysis instead of spatiotemporal
rcs_stats.py /run/dir /base/dir \
    --analysis_type temporal
```

```text
# Custom significance level
rcs_stats.py /run/dir /base/dir \
    --test_type ks --alpha 0.001
```

```text
# Combine options
rcs_stats.py /run/dir /base/dir \
    --test_type mw --analysis_type temporal --alpha 0.005
```

**Available Options:**

- `run_dir`: Directory containing current run ensemble output files
- `base_dir`: Directory containing baseline ensemble output files

**Statistical Test Selection:**

`--test_type`: Statistical test identifier (default: `ks`)

- Distribution tests: `ks`, `ad`, `cvm`, `epps`, `energy`
- Location tests: `mw`, `ttest`, `brunner`
- Scale tests: `levene`, `ansari`, `mood`

`--analysis_type`: Analysis mode (default: `spatiotemporal`)

- `spatiotemporal`: Area-weighted global means
- `temporal`: Per-column temporal means

`--alpha`: Significance level (default: 0.01, or 0.001 for AD)

**Multiple Testing Correction:**

`--correction_method`: Correction method (default: `bonferroni`)

- `bonferroni`: Conservative, controls family-wise error rate
- `fdr`: False Discovery Rate (Benjamini-Hochberg), less conservative
- `none`: No correction applied

`--critical_fraction`: Fraction of sub-tests allowed to fail (default: 0.001)

**Failure Thresholds:**

`--max_failed_vars`: Maximum variables allowed to fail (default: 0)
`--magnitude_threshold`: Minimum relative difference threshold (default: None)

**File Pattern Customization:**

`--run_file_pattern`: File pattern for run ensemble files
(default: `*.scream_????.h.AVERAGE.*.nc`)

- Use `????` as placeholder for 4-digit instance number
- Supports wildcards (`*`) for flexible matching

`--base_file_pattern`: File pattern for baseline ensemble files
(default: `*.scream_????.h.AVERAGE.*.nc`)

- Use `????` as placeholder for 4-digit instance number
- Supports wildcards (`*`) for flexible matching

**Complete Examples:**

```text
# Basic usage with defaults
rcs_stats.py /run/dir /base/dir
```

```text
# High-sensitivity Anderson-Darling test
rcs_stats.py /run/dir /base/dir \
    --test_type ad
```

```text
# Use FDR correction for better power
rcs_stats.py /run/dir /base/dir \
    --test_type ks \
    --correction_method fdr \
    --alpha 0.05
```

```text
# Combine multiple options for custom analysis
rcs_stats.py /run/dir /base/dir \
    --test_type mw \
    --analysis_type temporal \
    --correction_method fdr \
    --alpha 0.01 \
    --max_failed_vars 3
```

```text
# Require practical significance (>1% difference)
rcs_stats.py /run/dir /base/dir \
    --test_type ks \
    --magnitude_threshold 0.01 \
    --correction_method none
```

```text
# Exploratory analysis with relaxed thresholds
rcs_stats.py /run/dir /base/dir \
    --test_type energy \
    --alpha 0.05 \
    --correction_method fdr \
    --max_failed_vars 10
```

```text
# Compare ensembles with different output formats
rcs_stats.py /run/dir /base/dir \
    --run_file_pattern "*.eam_????.h0.*.nc" \
    --base_file_pattern "*.scream_????.h.AVERAGE.*.nc"
```

```text
# Custom instance number patterns
rcs_stats.py /run/dir /base/dir \
    --run_file_pattern "output.????.nc" \
    --base_file_pattern "baseline.????.nc"
```

**Getting Help:**

```text
python rcs_stats.py --help
```

This displays complete documentation including all available tests with
descriptions, usage examples, and parameter explanations.
