#!/usr/bin/env python3

"""
Modular statistical comparison framework for RCS system test.

This module provides a comprehensive suite of two-sample statistical tests
for comparing ensemble climate model simulations to determine if they are
statistically equivalent.

Test Categories:
----------------
1. DISTRIBUTION TESTS: Compare entire probability distributions
   - Kolmogorov-Smirnov (ks): General-purpose, moderate sensitivity
   - Anderson-Darling (ad): High sensitivity, emphasizes tails
   - Cramer-von Mises (cvm): Moderate-high sensitivity
   - Epps-Singleton (epps): Detects location and scale differences
   - Energy Distance (energy): Powerful for any type of difference

2. LOCATION TESTS: Compare central tendencies (mean/median)
   - Mann-Whitney U (mw): Non-parametric, compares medians
   - Welch's t-test (ttest): Parametric, compares means
   - Brunner-Munzel (brunner): Robust alternative to t-test

3. SCALE TESTS: Compare variability/spread
   - Levene (levene): Tests equality of variances
   - Ansari-Bradley (ansari): Non-parametric scale test
   - Mood (mood): Non-parametric dispersion test

Usage:
------
From command line:
    python new_test.py /run/dir /base/dir --test_type ks

From Python code:
    from new_test import run_stats_comparison
    comments, status = run_stats_comparison(
        run_dir, base_dir,
        analysis_type='spatiotemporal',
        test_type='ks'
    )

Choosing a Test:
----------------
- For general comparison: Use 'ks' (Kolmogorov-Smirnov) - balanced sensitivity
- For detecting subtle differences: Use 'ad' (Anderson-Darling) - very sensitive
- For comparing means only: Use 'ttest' or 'mw'
- For comparing variability: Use 'levene' or 'ansari'
- When unsure about distribution: Use 'mw' or 'energy' (non-parametric)

References:
-----------
All tests implemented using scipy.stats:
https://docs.scipy.org/doc/scipy/reference/stats.html
"""

import os
import glob
import json
import sys
import logging
import warnings
from abc import ABC, abstractmethod

sys.path.append(os.path.join(os.path.dirname(__file__), "../../scripts"))

try:
    from utils import _ensure_pylib_impl

    _ensure_pylib_impl("xarray")
    _ensure_pylib_impl("dask")
    _ensure_pylib_impl("scipy")
    _ensure_pylib_impl("statsmodels", min_version="0.14.0")

    import numpy as np
    import xarray as xr
    from scipy import stats
except ImportError as e:
    raise ImportError(f"Could not ensure Python packages: {e}") from e

logger = logging.getLogger(__name__)


# ==========================================================
# Abstract Base Test
# ==========================================================

# pylint: disable=too-few-public-methods
class StatisticalTest(ABC):
    """Base class for statistical comparison tests."""

    def __init__(self, alpha=0.01, magnitude_threshold=None):
        self.alpha = alpha
        self.magnitude_threshold = magnitude_threshold

    @abstractmethod
    def _compute_test_statistic(self, a, b):
        """Return (statistic, pvalue)."""

    def compare(self, a, b):
        """Run the test with preprocessing and standardized output."""
        a_orig, b_orig = a, b
        a, b = self._clean_data(a, b)

        # Compute descriptive statistics
        desc_stats = self._compute_descriptive_stats(a, b, a_orig, b_orig)

        if len(a) < 2 or len(b) < 2:
            return {
                "statistic": np.nan,
                "pvalue": 1.0,
                "hypothesis": "PASS",
                **desc_stats,
                "reason": "Insufficient samples",
            }

        if np.allclose(a, b, equal_nan=True):
            return {
                "statistic": 0.0,
                "pvalue": 1.0,
                "hypothesis": "PASS",
                **desc_stats,
                "reason": "Samples identical",
            }

        stat, pval = self._compute_test_statistic(a, b)
        reject = pval < self.alpha

        # Check magnitude threshold
        magnitude_check = None
        if self.magnitude_threshold is not None:
            rel_diff = self._relative_diff(a, b)
            magnitude_check = {
                "relative_difference": float(rel_diff),
                "magnitude_threshold": self.magnitude_threshold,
                "exceeds_threshold": rel_diff > self.magnitude_threshold,
            }
            reject = reject and (rel_diff > self.magnitude_threshold)

        result = {
            "statistic": float(stat),
            "pvalue": float(pval),
            "hypothesis": "FAIL" if reject else "PASS",
            **desc_stats,
        }

        if magnitude_check:
            result["magnitude_check"] = magnitude_check

        # Add interpretive reason
        if reject:
            result["reason"] = self._explain_failure(pval, desc_stats)
        else:
            result["reason"] = self._explain_pass(pval, desc_stats)

        return result

    @staticmethod
    def _clean_data(a, b):
        """Remove NaNs and flatten arrays."""
        a = np.asarray(a).ravel()
        b = np.asarray(b).ravel()
        return a[~np.isnan(a)], b[~np.isnan(b)]

    @staticmethod
    def _relative_diff(a, b):
        """Compute relative mean difference."""
        mean1, mean2 = np.nanmean(a), np.nanmean(b)
        denom = (abs(mean1) + abs(mean2)) / 2.0 + 1e-20
        return abs(mean1 - mean2) / denom

    @staticmethod
    def _compute_descriptive_stats(a, b, a_orig, b_orig):
        """Compute comprehensive descriptive statistics for both samples."""
        return {
            "sample1": {
                "n": len(a),
                "n_original": len(np.asarray(a_orig).ravel()),
                "n_nan": len(np.asarray(a_orig).ravel()) - len(a),
                "mean": float(np.mean(a)) if len(a) > 0 else np.nan,
                "median": float(np.median(a)) if len(a) > 0 else np.nan,
                "std": float(np.std(a, ddof=1)) if len(a) > 1 else np.nan,
                "min": float(np.min(a)) if len(a) > 0 else np.nan,
                "max": float(np.max(a)) if len(a) > 0 else np.nan,
                "q25": float(np.percentile(a, 25)) if len(a) > 0 else np.nan,
                "q75": float(np.percentile(a, 75)) if len(a) > 0 else np.nan,
            },
            "sample2": {
                "n": len(b),
                "n_original": len(np.asarray(b_orig).ravel()),
                "n_nan": len(np.asarray(b_orig).ravel()) - len(b),
                "mean": float(np.mean(b)) if len(b) > 0 else np.nan,
                "median": float(np.median(b)) if len(b) > 0 else np.nan,
                "std": float(np.std(b, ddof=1)) if len(b) > 1 else np.nan,
                "min": float(np.min(b)) if len(b) > 0 else np.nan,
                "max": float(np.max(b)) if len(b) > 0 else np.nan,
                "q25": float(np.percentile(b, 25)) if len(b) > 0 else np.nan,
                "q75": float(np.percentile(b, 75)) if len(b) > 0 else np.nan,
            },
            "difference": {
                "mean_diff": (
                    float(np.mean(a) - np.mean(b))
                    if len(a) > 0 and len(b) > 0
                    else np.nan
                ),
                "mean_diff_pct": (
                    float(100 * (np.mean(a) - np.mean(b)) /
                          (abs(np.mean(b)) + 1e-20))
                    if len(a) > 0 and len(b) > 0
                    else np.nan
                ),
                "median_diff": (
                    float(np.median(a) - np.median(b))
                    if len(a) > 0 and len(b) > 0
                    else np.nan
                ),
                "std_ratio": (
                    float(np.std(a, ddof=1) / (np.std(b, ddof=1) + 1e-20))
                    if len(a) > 1 and len(b) > 1
                    else np.nan
                ),
            },
        }

    def _explain_failure(self, pval, desc_stats):
        """Generate human-readable explanation for test failure."""
        diff = desc_stats["difference"]

        parts = [f"p={pval:.4e} < alpha={self.alpha}"]

        # Add context about differences
        if not np.isnan(diff["mean_diff"]):
            parts.append(
                f"mean difference: {diff['mean_diff']:.4e} "
                f"({diff['mean_diff_pct']:.2f}%)"
            )

        if not np.isnan(diff["std_ratio"]) and abs(diff["std_ratio"] - 1) > 0.1:
            parts.append(f"std ratio: {diff['std_ratio']:.3f}")

        return "; ".join(parts)

    def _explain_pass(self, pval, desc_stats):
        """Generate human-readable explanation for test pass."""
        diff = desc_stats["difference"]

        parts = [f"p={pval:.4e} >= alpha={self.alpha}"]

        if not np.isnan(diff["mean_diff_pct"]):
            parts.append(f"mean diff: {diff['mean_diff_pct']:.2f}%")

        return "; ".join(parts)


# ==========================================================
# Concrete Implementations
# ==========================================================

# Distribution-Based Tests (compare entire distributions)
# --------------------------------------------------------


class KSTest(StatisticalTest):
    """
    Kolmogorov-Smirnov two-sample test.
    Tests if two samples come from the same distribution.
    - Sensitivity: Moderate
    - Best for: Detecting overall distribution differences
    - Assumption: Continuous distributions
    """

    def _compute_test_statistic(self, a, b):
        # Suppress expected warning about switching to asymptotic method
        # This is normal for large samples and doesn't affect validity
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                'ks_2samp: Exact calculation unsuccessful'
            )
            return stats.ks_2samp(a, b)


class AndersonDarlingTest(StatisticalTest):
    """
    Anderson-Darling k-sample test.
    More sensitive than K-S, especially to differences in tails.
    - Sensitivity: High (especially in distribution tails)
    - Best for: Detecting subtle distributional differences
    - Assumption: Continuous distributions
    - Note: Returns approximate p-values in range [0.001, 0.25]
    """

    def _compute_test_statistic(self, a, b):
        if np.unique(a).size < 2 or np.unique(b).size < 2:
            return np.nan, 1.0

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            res = stats.anderson_ksamp([a, b])

        # anderson_ksamp returns significance_level as a percentage (0.1 to 25)
        # Convert to p-value (0.001 to 0.25). Values outside this range are
        # capped by scipy, so we respect those bounds.
        pvalue = res.significance_level / 100.0

        # Note: AD test is inherently more sensitive than KS.
        # Consider using a stricter alpha (e.g., 0.0001) for AD tests.
        return float(res.statistic), float(pvalue)


class CramerVonMisesTest(StatisticalTest):
    """
    Cramér-von Mises two-sample test.
    Similar to K-S but gives more weight to tail differences.
    - Sensitivity: Moderate to High
    - Best for: Overall distribution comparison with tail sensitivity
    - Assumption: Continuous distributions
    """

    def _compute_test_statistic(self, a, b):
        res = stats.cramervonmises_2samp(a, b)
        return res.statistic, res.pvalue


class EppsSingletonTest(StatisticalTest):
    """
    Epps-Singleton two-sample test.
    Tests for differences in both location and scale.
    - Sensitivity: Moderate
    - Best for: Detecting differences in mean AND variance
    - Assumption: Works well for non-normal distributions
    """

    def _compute_test_statistic(self, a, b):
        res = stats.epps_singleton_2samp(a, b)
        return res.statistic, res.pvalue


# Location Tests (compare medians/means)
# ----------------------------------------


class MannWhitneyUTest(StatisticalTest):
    """
    Mann-Whitney U test (Wilcoxon rank-sum test).
    Non-parametric test for difference in central tendency.
    - Sensitivity: Moderate
    - Best for: Comparing medians of non-normal distributions
    - Assumption: None (distribution-free)
    """

    def _compute_test_statistic(self, a, b):
        res = stats.mannwhitneyu(a, b, alternative="two-sided")
        return res.statistic, res.pvalue


class TTest(StatisticalTest):
    """
    Welch's t-test (unequal variance t-test).
    Parametric test for difference in means.
    - Sensitivity: Moderate to High (for mean differences)
    - Best for: Comparing means of approximately normal distributions
    - Assumption: Approximately normal distributions (robust to violations)
    """

    def _compute_test_statistic(self, a, b):
        res = stats.ttest_ind(a, b, equal_var=False)
        return res.statistic, res.pvalue


# Scale/Variance Tests
# ----------------------


class LeveneTest(StatisticalTest):
    """
    Levene's test for equality of variances.
    Tests if two samples have equal variances.
    - Sensitivity: Moderate
    - Best for: Detecting differences in variability/spread
    - Assumption: Robust to non-normality
    """

    def _compute_test_statistic(self, a, b):
        stat, pval = stats.levene(a, b)
        return stat, pval


class AnsariBradleyTest(StatisticalTest):
    """
    Ansari-Bradley test for equal scale parameters.
    Non-parametric test for differences in spread.
    - Sensitivity: Moderate
    - Best for: Comparing variability when distributions have similar medians
    - Assumption: Samples differ primarily in scale, not location
    """

    def _compute_test_statistic(self, a, b):
        res = stats.ansari(a, b)
        return res.statistic, res.pvalue


class MoodTest(StatisticalTest):
    """
    Mood's test for equal scale parameters.
    Non-parametric alternative to variance comparison.
    - Sensitivity: Moderate
    - Best for: Comparing spread/dispersion
    - Assumption: None (distribution-free)
    """

    def _compute_test_statistic(self, a, b):
        res = stats.mood(a, b)
        return res.statistic, res.pvalue


# Energy-Based Test
# ------------------


class EnergyDistanceTest(StatisticalTest):
    """
    Energy distance test (statistical energy).
    Powerful test based on energy statistics.
    - Sensitivity: High
    - Best for: Detecting any type of distributional difference
    - Assumption: None (works for any distribution)
    - Note: Computationally intensive for large samples
    """

    def _compute_test_statistic(self, a, b):
        res = stats.energy_distance(a, b)
        # energy_distance returns only the statistic, need permutation
        # for p-value. Use bootstrap approximation.
        n_permutations = 1000
        observed = res

        # Combine samples for permutation test
        combined = np.concatenate([a, b])
        n_a = len(a)

        # Permutation test
        count = 0
        for _ in range(n_permutations):
            perm_combined = np.random.permutation(combined)
            perm_a = perm_combined[:n_a]
            perm_b = perm_combined[n_a:]
            perm_stat = stats.energy_distance(perm_a, perm_b)
            if perm_stat >= observed:
                count += 1

        pvalue = (count + 1) / (n_permutations + 1)
        return observed, pvalue


# Quantile-Based Test
# ---------------------


class BrunnerMunzelTest(StatisticalTest):
    """
    Brunner-Munzel test (generalized Wilcoxon test).
    Non-parametric test for stochastic equality.
    - Sensitivity: Moderate to High
    - Best for: Robust alternative to t-test for ordinal data
    - Assumption: None (distribution-free)
    """

    def _compute_test_statistic(self, a, b):
        res = stats.brunnermunzel(a, b)
        return res.statistic, res.pvalue


# ==========================================================
# Test Factory
# ==========================================================


# pylint: disable=too-many-return-statements
def get_test(test_type: str, alpha=0.01, magnitude_threshold=None):
    """
    Factory function for selecting statistical test type.

    Available tests organized by category:

    DISTRIBUTION TESTS (compare entire distributions):
    - 'ks': Kolmogorov-Smirnov (recommended default, moderate sensitivity)
    - 'ad': Anderson-Darling (high sensitivity, especially in tails)
    - 'cvm': Cramér-von Mises (moderate-high sensitivity)
    - 'epps': Epps-Singleton (detects location and scale differences)
    - 'energy': Energy distance (high sensitivity, any difference type)

    LOCATION TESTS (compare means/medians):
    - 'mw': Mann-Whitney U (non-parametric, compares medians)
    - 'ttest': Welch's t-test (parametric, compares means)
    - 'brunner': Brunner-Munzel (robust alternative to t-test)

    SCALE TESTS (compare variances/spread):
    - 'levene': Levene's test (variance equality)
    - 'ansari': Ansari-Bradley (non-parametric scale test)
    - 'mood': Mood's test (non-parametric dispersion test)

    Recommended alpha values:
    - Most tests: 0.01 (default)
    - Anderson-Darling: 0.001 (auto-adjusted due to high sensitivity)
    - Energy distance: 0.01 (already uses permutation testing)
    """
    test_type = test_type.lower()

    # Distribution tests
    if test_type in ("ks", "kolmogorov-smirnov"):
        return KSTest(alpha, magnitude_threshold)
    if test_type in ("ad", "anderson-darling"):
        ad_alpha = alpha if alpha != 0.01 else 0.001
        return AndersonDarlingTest(ad_alpha, magnitude_threshold)
    if test_type in ("cvm", "cm", "cramer-von-mises", "cramer"):
        return CramerVonMisesTest(alpha, magnitude_threshold)
    if test_type in ("epps", "epps-singleton"):
        return EppsSingletonTest(alpha, magnitude_threshold)
    if test_type in ("energy", "energy-distance"):
        return EnergyDistanceTest(alpha, magnitude_threshold)

    # Location tests
    if test_type in ("mw", "mannwhitney", "mann-whitney"):
        return MannWhitneyUTest(alpha, magnitude_threshold)
    if test_type in ("ttest", "t-test", "welch"):
        return TTest(alpha, magnitude_threshold)
    if test_type in ("brunner", "brunnermunzel", "brunner-munzel"):
        return BrunnerMunzelTest(alpha, magnitude_threshold)

    # Scale tests
    if test_type in ("levene",):
        return LeveneTest(alpha, magnitude_threshold)
    if test_type in ("ansari", "ansari-bradley"):
        return AnsariBradleyTest(alpha, magnitude_threshold)
    if test_type in ("mood",):
        return MoodTest(alpha, magnitude_threshold)

    # Error message with helpful suggestions
    raise ValueError(
        f"Unknown test type: '{test_type}'\n"
        f"Available tests:\n"
        f"  Distribution: ks, ad, cvm, epps, energy\n"
        f"  Location: mw, ttest, brunner\n"
        f"  Scale: levene, ansari, mood\n"
        f"Run with --help for detailed descriptions."
    )


# ==========================================================
# Multiple Testing Correction
# ==========================================================


# pylint: disable=too-many-branches
def apply_multiple_testing_correction(results, alpha, method="bonferroni"):
    """
    Apply multiple testing correction to p-values.

    Args:
        results: Dictionary of {variable: test_result_dict}
        alpha: Significance level
        method: Correction method - 'bonferroni', 'fdr', or 'none'

    Returns:
        Updated results dictionary with corrected hypothesis decisions
    """
    if method == "none":
        return results

    # Extract p-values and variable names
    var_names = []
    pvalues = []
    for var, res in results.items():
        if "pvalue" in res and not np.isnan(res["pvalue"]):
            var_names.append(var)
            pvalues.append(res["pvalue"])

    if len(pvalues) == 0:
        return results

    n_tests = len(pvalues)
    pvalues = np.array(pvalues)

    if method == "bonferroni":
        # Bonferroni correction: divide alpha by number of tests
        corrected_alpha = alpha / n_tests
        reject = pvalues < corrected_alpha

        # Update results with corrected decisions
        for i, var in enumerate(var_names):
            original_result = results[var]["hypothesis"]
            results[var]["corrected_alpha"] = corrected_alpha
            results[var]["correction_method"] = "bonferroni"

            # Re-evaluate hypothesis with corrected alpha
            if reject[i]:
                results[var]["hypothesis"] = "FAIL"
                if original_result == "PASS":
                    results[var]["reason"] += " (failed after Bonferroni correction)"
            else:
                results[var]["hypothesis"] = "PASS"
                if original_result == "FAIL":
                    results[var]["reason"] += " (passed after Bonferroni correction)"

    elif method == "fdr":
        # Benjamini-Hochberg FDR correction
        # Sort p-values in ascending order
        sorted_indices = np.argsort(pvalues)
        sorted_pvalues = pvalues[sorted_indices]

        # Find largest i where p(i) <= (i/m) * alpha
        reject = np.zeros(n_tests, dtype=bool)
        for i in range(n_tests - 1, -1, -1):
            if sorted_pvalues[i] <= ((i + 1) / n_tests) * alpha:
                # Reject this and all smaller p-values
                reject[sorted_indices[: i + 1]] = True
                break

        # Update results with FDR-corrected decisions
        for i, var in enumerate(var_names):
            original_result = results[var]["hypothesis"]
            critical_value = ((np.where(sorted_indices == i)
                              [0][0] + 1) / n_tests) * alpha
            results[var]["fdr_critical_value"] = float(critical_value)
            results[var]["correction_method"] = "fdr"

            # Re-evaluate hypothesis with FDR
            if reject[i]:
                results[var]["hypothesis"] = "FAIL"
                if original_result == "PASS":
                    results[var]["reason"] += " (failed after FDR correction)"
            else:
                results[var]["hypothesis"] = "PASS"
                if original_result == "FAIL":
                    results[var]["reason"] += " (passed after FDR correction)"

    return results


# ==========================================================
# Main entry: run_stats_comparison
# ==========================================================


# pylint: disable=too-many-locals, too-many-statements
# pylint: disable=too-many-arguments, too-many-positional-arguments
def run_stats_comparison(
    run_dir,
    base_dir,
    analysis_type="spatiotemporal",
    test_type="ks",
    alpha=None,
    critical_fraction=0.001,
    correction_method="bonferroni",
    max_failed_vars=0,
    magnitude_threshold=None,
    run_file_pattern=None,
    base_file_pattern=None,
):
    """
    Compare ensembles using configurable statistical tests.

    Args:
        run_dir: Directory containing run ensemble output
        base_dir: Directory containing baseline ensemble output
        analysis_type: 'spatiotemporal' (area-weighted means) or
                      'temporal' (per-column means)
        test_type: Statistical test identifier:
            - 'ks': Kolmogorov-Smirnov (recommended, moderate sensitivity)
            - 'mw': Mann-Whitney U (non-parametric, moderate sensitivity)
            - 'ad': Anderson-Darling (very sensitive, use with caution)
            - 'cvm': Cramér-von Mises (moderate to high sensitivity)
        alpha: Significance level (default: 0.01 for most tests, 0.001 for AD)
        critical_fraction: Fraction of variables that can fail before overall
                          test fails (default: 0.001)
        correction_method: Multiple testing correction method:
            - 'bonferroni': Conservative, controls family-wise error rate
            - 'fdr': False Discovery Rate (Benjamini-Hochberg), less conservative
            - 'none': No correction applied
        max_failed_vars: Maximum number of variables allowed to fail
                        (default: 0)
        magnitude_threshold: Minimum relative difference to consider
                            significant (default: None, all differences count)
        run_file_pattern: File pattern for run ensemble files. Use '????' as
                         placeholder for instance number
                         (default: '*.scream_????.h.AVERAGE.*.nc')
        base_file_pattern: File pattern for baseline ensemble files.
                          Use '????' as placeholder for instance number
                          (default: '*.scream_????.h.AVERAGE.*.nc')

    Note: Anderson-Darling is significantly more sensitive than K-S and will
          detect smaller distributional differences. It uses a stricter default
          alpha (0.001 vs 0.01) to compensate.
    """
    # Configuration parameters
    ALPHA = alpha if alpha is not None else 0.01
    CRITICAL_FRACTION = critical_fraction
    CORRECTION_METHOD = correction_method.lower()
    MAX_FAILED_VARS = max_failed_vars
    MAGNITUDE_THRESHOLD = magnitude_threshold
    RUN_FILE_PATTERN = run_file_pattern
    BASE_FILE_PATTERN = base_file_pattern

    # Map test type to full name for reporting
    test_names = {
        "ks": "Kolmogorov-Smirnov",
        "ad": "Anderson-Darling",
        "cvm": "Cramer-von Mises",
        "epps": "Epps-Singleton",
        "energy": "Energy Distance",
        "mw": "Mann-Whitney U",
        "ttest": "Welch's t-test",
        "brunner": "Brunner-Munzel",
        "levene": "Levene",
        "ansari": "Ansari-Bradley",
        "mood": "Mood",
    }
    test_full_name = test_names.get(test_type.lower(), test_type.upper())

    # Map analysis type to full description
    analysis_names = {
        "spatiotemporal": "Spatiotemporal (area-weighted global means)",
        "temporal": "Temporal (per-column means)",
    }
    analysis_full_name = analysis_names.get(
        analysis_type, analysis_type
    )

    # Map correction method to full description
    correction_names = {
        "bonferroni": "Bonferroni (family-wise error rate control)",
        "fdr": "FDR - Benjamini-Hochberg (false discovery rate control)",
        "none": "None (no multiple testing correction)",
    }
    correction_full_name = correction_names.get(
        CORRECTION_METHOD, CORRECTION_METHOD.upper()
    )

    comments = [
        "",
        "=" * 70,
        "TEST CONFIGURATION",
        "=" * 70,
        f"Statistical Test: {test_full_name}",
        "  (Two-sample test comparing ensemble distributions)",
        "",
        f"Analysis Type: {analysis_full_name}",
        "  (Method for aggregating spatial and temporal data)",
        "",
        f"Significance Level (Alpha): {ALPHA}",
        "  (Probability threshold for rejecting null hypothesis)",
        "",
        f"Multiple Testing Correction: {correction_full_name}",
        "  (Adjusts p-value thresholds when testing many variables)",
        "",
        f"Critical Fraction: {CRITICAL_FRACTION}",
        "  (Maximum fraction of sub-tests allowed to fail per variable)",
        "",
        f"Max Failed Variables: {MAX_FAILED_VARS}",
        "  (Maximum number of variables that can fail before test fails)",
    ]
    if MAGNITUDE_THRESHOLD is not None:
        comments.extend([
            "",
            f"Magnitude Threshold: {MAGNITUDE_THRESHOLD}",
            "  (Minimum relative difference required for practical ",
            "significance)",
        ])

    comments.extend([
        "",
        f"Run File Pattern: {RUN_FILE_PATTERN}",
        "  (Glob pattern for run ensemble files)",
        "",
        f"Baseline File Pattern: {BASE_FILE_PATTERN}",
        "  (Glob pattern for baseline ensemble files)",
    ])

    # Load ensemble datasets using configurable patterns
    run_files = glob.glob(os.path.join(run_dir, RUN_FILE_PATTERN))
    base_files = glob.glob(os.path.join(base_dir, BASE_FILE_PATTERN))

    if not run_files or not base_files:
        raise FileNotFoundError(
            f"Missing scream output in {run_dir} or {base_dir}")

    run_instances = _extract_instance_nums(run_files, RUN_FILE_PATTERN)
    base_instances = _extract_instance_nums(base_files, BASE_FILE_PATTERN)

    comments.extend([
        "",
        "=" * 70,
        "ENSEMBLE INFORMATION",
        "=" * 70,
        f"Run Ensemble Instances: {len(run_instances)}",
        "  (Number of perturbed ensemble members in current run)",
        "",
        f"Baseline Ensemble Instances: {len(base_instances)}",
        "  (Number of perturbed ensemble members in baseline)",
    ])

    ensemble_1 = {
        i: xr.open_mfdataset(
            os.path.join(run_dir, RUN_FILE_PATTERN.replace("????", i)),
            decode_times=False,
            data_vars="all",
        )
        for i in run_instances
    }
    ensemble_2 = {
        i: xr.open_mfdataset(
            os.path.join(base_dir, BASE_FILE_PATTERN.replace("????", i)),
            decode_times=False,
            data_vars="all",
        )
        for i in base_instances
    }

    sample_ds = list(ensemble_1.values())[0]
    test_vars = _get_testable_variables(sample_ds)

    comments.extend([
        "",
        f"Variables to Test: {len(test_vars)}",
        "  (Number of suitable variables with time and spatial dimensions)",
        "",
        "=" * 70,
    ])

    area_weights = (
        _get_area_weights(
            sample_ds) if analysis_type == "spatiotemporal" else None
    )
    test_obj = get_test(test_type, alpha=ALPHA,
                        magnitude_threshold=MAGNITUDE_THRESHOLD)
    results = {}

    # Simplified test loop
    for var in test_vars:
        try:
            a, b = _extract_data_for_var(
                var, ensemble_1, ensemble_2, analysis_type, area_weights
            )
            results[var] = test_obj.compare(a, b)
        except (ValueError, KeyError, OSError, IndexError, TypeError) as e:
            # Log expected data- or IO-related errors and continue testing
            # other variables; allow truly unexpected exceptions to propagate
            # for visibility.
            logger.warning("Error testing %s: %s", var, e)

    # Apply multiple testing correction if requested
    if CORRECTION_METHOD != "none":
        results = apply_multiple_testing_correction(
            results, ALPHA, method=CORRECTION_METHOD
        )

    failed = [v for v, r in results.items() if r["hypothesis"] == "FAIL"]
    passed = [v for v, r in results.items() if r["hypothesis"] == "PASS"]

    # Build detailed comments with statistics
    comments.append("\n" + "=" * 70)
    comments.append("TEST RESULTS SUMMARY")
    comments.append("=" * 70)
    comments.append(
        f"Variables Passed: {len(passed)} "
        f"({100.0 * len(passed) / len(results):.1f}%)"
    )
    comments.append(
        "  (Variables where null hypothesis was NOT rejected)"
    )
    comments.append("")
    comments.append(
        f"Variables Failed: {len(failed)} "
        f"({100.0 * len(failed) / len(results):.1f}%)"
    )
    comments.append(
        "  (Variables where distributions are statistically different)"
    )

    # Show details for failed variables
    if failed:
        comments.append("\n" + "=" * 70)
        comments.append("FAILED VARIABLES (detailed)")
        comments.append("=" * 70)
        for var in failed[:10]:  # Limit to first 10
            r = results[var]
            comments.append(f"\n{var}:")
            comments.append(f"  Result: {r.get('reason', 'Failed')}")

            # Show sample statistics
            s1 = r.get("sample1", {})
            s2 = r.get("sample2", {})
            diff = r.get("difference", {})

            if s1 and s2:
                comments.append(
                    f"  Sample 1 (run):      "
                    f"n={s1.get('n', 0):6d}, "
                    f"mean={s1.get('mean', np.nan):12.6e}, "
                    f"std={s1.get('std', np.nan):12.6e}"
                )
                comments.append(
                    f"  Sample 2 (baseline): "
                    f"n={s2.get('n', 0):6d}, "
                    f"mean={s2.get('mean', np.nan):12.6e}, "
                    f"std={s2.get('std', np.nan):12.6e}"
                )

            if diff:
                comments.append(
                    f"  Difference: "
                    f"Δmean={diff.get('mean_diff', np.nan):12.6e} "
                    f"({diff.get('mean_diff_pct', np.nan):+.2f}%), "
                    f"std_ratio={diff.get('std_ratio', np.nan):.3f}"
                )

        if len(failed) > 10:
            comments.append(
                f"\n  ... and {len(failed) - 10} more failed variables")

    # Show sample of passed variables
    if passed and len(passed) <= 10:
        comments.append("\n" + "=" * 70)
        comments.append("PASSED VARIABLES (sample)")
        comments.append("=" * 70)
        for var in passed[:5]:
            r = results[var]
            diff = r.get("difference", {})
            comments.append(
                f"  {var}: {r.get('reason', 'Passed')} | "
                f"Δmean={diff.get('mean_diff_pct', np.nan):+.2f}%"
            )

    test_status = "PASS" if len(failed) <= MAX_FAILED_VARS else "FAIL"

    output_file = os.path.join(run_dir, f"{test_type}_test_results.json")
    with open(output_file, "w", encoding="utf-8") as f:
        output_data = {
            "test_type": test_type,
            "analysis_type": analysis_type,
            "configuration": {
                "alpha": ALPHA,
                "correction_method": CORRECTION_METHOD,
                "critical_fraction": CRITICAL_FRACTION,
                "max_failed_vars": MAX_FAILED_VARS,
            },
            "summary": {
                "passed": len(passed),
                "failed": len(failed),
                "total": len(results),
                "test_status": test_status,
            },
            "failed_variables": failed,
            "passed_variables": passed,
            "details": results,
        }
        if MAGNITUDE_THRESHOLD is not None:
            output_data["configuration"]["magnitude_threshold"] = MAGNITUDE_THRESHOLD
        json.dump(output_data, f, indent=2)
    comments.append(f"\nDetailed results saved to: {output_file}")

    return "\n".join(comments), test_status


# ==========================================================
# Shared helper functions (simplified versions)
# ==========================================================


def _extract_instance_nums(files, pattern):
    """
    Extract instance numbers from filenames based on a pattern.

    Args:
        files: List of file paths
        pattern: File pattern with '????' as placeholder for instance number

    Returns:
        Sorted list of unique instance numbers as 4-digit strings
    """
    instance_nums = []

    # Find the position of '????' in the pattern
    if "????" not in pattern:
        # If no placeholder, try legacy method
        for f in files:
            parts = os.path.basename(f).split(".scream_")
            if len(parts) > 1:
                inst = parts[1].split(".")[0]
                if inst.isdigit() and len(inst) == 4:
                    instance_nums.append(inst)
        return sorted(set(instance_nums))

    # Split pattern to get parts before and after the placeholder
    parts = pattern.split("????")
    if len(parts) != 2:
        raise ValueError(
            f"Pattern must contain exactly one '????' placeholder: {pattern}"
        )

    prefix, suffix = parts

    for f in files:
        basename = os.path.basename(f)

        # Remove the prefix
        if prefix.startswith("*"):
            # Handle wildcard prefix
            prefix_pattern = prefix.lstrip("*")
            if prefix_pattern and prefix_pattern not in basename:
                continue
            # Find where the fixed part of prefix starts
            if prefix_pattern:
                idx = basename.find(prefix_pattern)
                if idx == -1:
                    continue
                start_pos = idx + len(prefix_pattern)
            else:
                start_pos = 0
        else:
            if not basename.startswith(prefix):
                continue
            start_pos = len(prefix)

        # Extract potential instance number
        remaining = basename[start_pos:]

        # Remove the suffix
        if suffix.startswith("."):
            # Find the suffix part
            suffix_parts = suffix.split("*")
            # Get the fixed part after the instance number
            fixed_suffix = suffix_parts[0] if suffix_parts else suffix

            if fixed_suffix and fixed_suffix in remaining:
                end_pos = remaining.find(fixed_suffix)
                inst = remaining[:end_pos]
            else:
                # Try to extract first 4 digits
                inst = remaining[:4] if len(remaining) >= 4 else remaining
        else:
            inst = remaining[:4] if len(remaining) >= 4 else remaining

        # Validate instance number
        if inst.isdigit() and len(inst) == 4:
            instance_nums.append(inst)

    return sorted(set(instance_nums))


def _extract_data_for_var(var, ensemble_1, ensemble_2, analysis_type, area_weights):
    """Extract comparable arrays for given variable and analysis type."""
    data1, data2 = [], []

    for ens, coll in ((ensemble_1, data1), (ensemble_2, data2)):
        for ds in coll_from_dict(ens):
            if var not in ds:
                continue
            v = _prepare_variable_data(ds[var])

            if analysis_type == "spatiotemporal":
                if area_weights is not None and "ncol" in v.dims:
                    arr = (v * area_weights).sum(dim="ncol", skipna=True).values
                else:
                    spatial_dims = [d for d in v.dims if d != "time"]
                    arr = v.mean(dim=spatial_dims, skipna=True).values
            else:  # temporal
                # Compute temporal mean along time axis (axis=0)
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', 'Mean of empty slice')
                    arr = np.nanmean(v.values, axis=0)

                # Remove any remaining dimensions and flatten
                arr = arr.ravel()

                # Filter out NaN values for this instance
                arr = arr[~np.isnan(arr)]

            coll.append(arr)

    # For temporal analysis, we don't concatenate because arrays may have
    # different lengths due to masking. Instead, we need to ensure we're
    # comparing the same spatial locations.
    if analysis_type == "temporal":
        # Check if all arrays have the same length
        lengths1 = [len(a) for a in data1]
        lengths2 = [len(a) for a in data2]

        if len(set(lengths1)) > 1 or len(set(lengths2)) > 1:
            # Variable length arrays - likely due to different masking
            # Use the intersection of valid points
            min1 = min(lengths1)
            min2 = min(lengths2)
            min_len = min(min1, min2)
            logger.debug(
                "Variable %s has inconsistent spatial dimensions. "
                "Ensemble 1 lengths: %s, Ensemble 2 lengths: %s. "
                "Truncating to %d points.",
                var, lengths1, lengths2, min_len
            )
            data1 = [a[:min_len] for a in data1]
            data2 = [a[:min_len] for a in data2]

    # Check if we have any data after processing
    if not data1 or not data2:
        raise ValueError(f"No valid data for variable {var}")

    return np.concatenate(data1), np.concatenate(data2)


def coll_from_dict(d):
    """Helper for dataset iteration."""
    return [ds for _, ds in sorted(d.items())]


def _get_testable_variables(dataset):
    skip_vars = {"time", "lat", "lon", "ncol", "lev", "ilev", "area"}
    vars_out = []
    for v in dataset.data_vars:
        if v in skip_vars or "time" not in dataset[v].dims:
            continue
        try:
            arr = dataset[v].values
            if not np.all(np.isnan(arr)) and not np.allclose(arr, arr.flat[0]):
                vars_out.append(v)
        except (ValueError, TypeError):
            continue
    return vars_out


def _get_area_weights(dataset):
    if "area" in dataset:
        area = dataset["area"].values
        return area / np.nansum(area)
    return None


def _prepare_variable_data(var):
    if "lev" in var.dims:
        var = var.mean(dim="lev", skipna=True)
    if "ilev" in var.dims:
        var = var.mean(dim="ilev", skipna=True)
    return var


###############################################################################
def parse_command_line(args, description):
###############################################################################
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        usage="""\n{0} run_dir base_dir [options]
OR
{0} --help

\033[1mEXAMPLES:\033[0m
    \033[1;32m# Default: KS test with Bonferroni correction\033[0m
    > {0} /path/to/run /path/to/baseline

    \033[1;32m# Anderson-Darling with temporal analysis\033[0m
    > {0} /path/to/run /path/to/baseline --test_type ad \\
          --analysis_type temporal

    \033[1;32m# Custom significance level with FDR correction\033[0m
    > {0} /path/to/run /path/to/baseline --test_type ks \\
          --alpha 0.001 --correction_method fdr

    \033[1;32m# No multiple testing correction\033[0m
    > {0} /path/to/run /path/to/baseline \\
          --correction_method none --magnitude_threshold 0.01

    \033[1;32m# Custom file patterns\033[0m
    > {0} /path/to/run /path/to/baseline \\
          --run_file_pattern "*.eam_????.h0.*.nc" \\
          --base_file_pattern "*.scream_????.h.AVERAGE.*.nc"
""".format(Path(args[0]).name),
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
\033[1mAVAILABLE STATISTICAL TESTS:\033[0m

  \033[1mDISTRIBUTION TESTS\033[0m (compare entire distributions):
    ks          Kolmogorov-Smirnov (recommended default)
    ad          Anderson-Darling (high sensitivity, especially tails)
    cvm         Cramér-von Mises (moderate-high sensitivity)
    epps        Epps-Singleton (location + scale)
    energy      Energy distance (powerful, any difference)

  \033[1mLOCATION TESTS\033[0m (compare means/medians):
    mw          Mann-Whitney U (non-parametric median test)
    ttest       Welch's t-test (parametric mean test)
    brunner     Brunner-Munzel (robust alternative to t-test)

  \033[1mSCALE TESTS\033[0m (compare variances/spread):
    levene      Levene's test (variance equality)
    ansari      Ansari-Bradley (non-parametric scale)
    mood        Mood's test (non-parametric dispersion)
""",
    )

    parser.add_argument("run_dir", type=str,
                        help="Directory of the new run ensemble")
    parser.add_argument("base_dir", type=str,
                        help="Directory of the baseline ensemble")
    parser.add_argument(
        "--analysis_type",
        type=str,
        default="spatiotemporal",
        choices=["spatiotemporal", "temporal"],
        help="Analysis type: spatiotemporal (area-weighted "
        "global means) or temporal (per-column means)",
    )
    parser.add_argument(
        "--test_type",
        type=str,
        default="ks",
        choices=[
            "ks",
            "ad",
            "cvm",
            "epps",
            "energy",
            "mw",
            "ttest",
            "brunner",
            "levene",
            "ansari",
            "mood",
        ],
        help="Statistical test identifier (default: ks)",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=None,
        help="Significance level (default: 0.01 for most, " "0.001 for AD)",
    )
    parser.add_argument(
        "--critical_fraction",
        type=float,
        default=0.001,
        help="Fraction of variables allowed to fail (default: 0.001)",
    )
    parser.add_argument(
        "--correction_method",
        type=str,
        default="bonferroni",
        choices=["bonferroni", "fdr", "none"],
        help="Multiple testing correction: bonferroni (conservative), "
        "fdr (Benjamini-Hochberg), or none (default: bonferroni)",
    )
    parser.add_argument(
        "--max_failed_vars",
        type=int,
        default=0,
        help="Maximum number of variables allowed to fail (default: 0)",
    )
    parser.add_argument(
        "--magnitude_threshold",
        type=float,
        default=None,
        help="Minimum relative difference to consider significant "
        "(default: None, all differences count)",
    )
    parser.add_argument(
        "--run_file_pattern",
        type=str,
        default="*.scream_????.h.AVERAGE.*.nc",
        help="File pattern for run ensemble (use ???? for instance number, "
        "default: *.scream_????.h.AVERAGE.*.nc)",
    )
    parser.add_argument(
        "--base_file_pattern",
        type=str,
        default="*.scream_????.h.AVERAGE.*.nc",
        help="File pattern for baseline ensemble (use ???? for instance "
        "number, default: *.scream_????.h.AVERAGE.*.nc)",
    )

    return parser.parse_args(args[1:])


###############################################################################
def _main_func(description):
###############################################################################
    cli_comments, cli_status = run_stats_comparison(
        **vars(parse_command_line(sys.argv, description))
    )

    print("\n")
    print("=" * 70)
    print(f"OVERALL TEST STATUS: {cli_status}")
    print("=" * 70)

    print("\n")
    print("=" * 70)
    print("DETAILS BELOW ...")
    print("=" * 70)
    print("\n")
    print(cli_comments)


###############################################################################

if (__name__ == "__main__"):
    _main_func(__doc__)
