#include "catch2/catch.hpp"

#include "p3_functions.hpp"
#include "p3_test_data.hpp"
#include "p3_unit_tests_common.hpp"

#include "share/core/eamxx_types.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <limits>       // std::numeric_limits

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3CloudWaterAutoconversion : public UnitWrap::UnitTest<D>::Base
{

void cloud_water_autoconversion_unit_bfb_tests() {

  CloudWaterAutoconversionData cwadc[max_pack_size] = {
    // rho, qc_incld, nc_incld, inv_qc_relvar
    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},

    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},

    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},

    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},
  };

  // Sync to device
  view_1d<CloudWaterAutoconversionData> cwadc_device("cwadc", max_pack_size);
  auto cwadc_host = Kokkos::create_mirror_view(cwadc_device);

  // This copy only copies the input variables.
  std::copy(&cwadc[0], &cwadc[0] + max_pack_size, cwadc_host.data());
  Kokkos::deep_copy(cwadc_device, cwadc_host);

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < max_pack_size; ++i) {
      cwadc[i].read(Base::m_ifile);
    }
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack rho, inv_rho, qc_incld, nc_incld, qr_incld, mu_c, nu, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, inv_qc_relvar;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      rho[s]                 = cwadc_device(vs).rho;
      qc_incld[s]            = cwadc_device(vs).qc_incld;
      nc_incld[s]            = cwadc_device(vs).nc_incld;
      inv_qc_relvar[s]       = cwadc_device(vs).inv_qc_relvar;
      qc2qr_autoconv_tend[s] = cwadc_device(vs).qc2qr_autoconv_tend;
      nc2nr_autoconv_tend[s] = cwadc_device(vs).nc2nr_autoconv_tend;
      ncautr[s]              = cwadc_device(vs).ncautr;
    }

    Functions::cloud_water_autoconversion(
        rho, qc_incld, nc_incld, inv_qc_relvar, qc2qr_autoconv_tend,
        nc2nr_autoconv_tend, ncautr,
        p3::Functions<Real,DefaultDevice>::P3Runtime());

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      cwadc_device(vs).rho                 = rho[s];
      cwadc_device(vs).qc_incld            = qc_incld[s];
      cwadc_device(vs).nc_incld            = nc_incld[s];
      cwadc_device(vs).inv_qc_relvar       = inv_qc_relvar[s];
      cwadc_device(vs).qc2qr_autoconv_tend = qc2qr_autoconv_tend[s];
      cwadc_device(vs).nc2nr_autoconv_tend = nc2nr_autoconv_tend[s];
      cwadc_device(vs).ncautr              = ncautr[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(cwadc_host, cwadc_device);

  // Validate results
  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(cwadc[s].rho                  == cwadc_host(s).rho);
      REQUIRE(cwadc[s].qc_incld             == cwadc_host(s).qc_incld);
      REQUIRE(cwadc[s].nc_incld             == cwadc_host(s).nc_incld);
      REQUIRE(cwadc[s].inv_qc_relvar        == cwadc_host(s).inv_qc_relvar);
      REQUIRE(cwadc[s].qc2qr_autoconv_tend  == cwadc_host(s).qc2qr_autoconv_tend);
      REQUIRE(cwadc[s].nc2nr_autoconv_tend  == cwadc_host(s).nc2nr_autoconv_tend);
      REQUIRE(cwadc[s].ncautr               == cwadc_host(s).ncautr);
    }
  }
  else if (this->m_baseline_action == GENERATE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      cwadc_host(s).write(Base::m_ofile);
    }
  }
}

  void run_bfb() {
    cloud_water_autoconversion_unit_bfb_tests();
  }

  void validate_autoconversion_parameters() {
    std::cout << "\n=== Parameter Validation ===\n";
    
    auto runtime = p3::Functions<Real,DefaultDevice>::P3Runtime();
    
    const Real expected_prefactor = 1350.0;     // s^-1
    const Real expected_qc_exp = 2.47;          // dimensionless
    const Real expected_nc_exp = 1.79;          // dimensionless
    const Real expected_radius = 25.0e-6;       // m (25 um)
    
    const Real tolerance = 0.01;  // 1% tolerance
    
    // Validate prefactor
    Real prefactor_error = std::abs(runtime.autoconversion_prefactor - expected_prefactor) 
                          / expected_prefactor;
    if (prefactor_error > tolerance) {
        std::cout << "WARNING: Prefactor mismatch\n";
        std::cout << "  Expected: " << expected_prefactor << " s^-1\n";
        std::cout << "  Actual: " << runtime.autoconversion_prefactor << " s^-1\n";
    } else {
        std::cout << "  Prefactor: " << runtime.autoconversion_prefactor << " s^-1\n";
    }
    
    // Validate qc exponent
    Real qc_exp_error = std::abs(runtime.autoconversion_qc_exponent - expected_qc_exp) 
                       / expected_qc_exp;
    if (qc_exp_error > tolerance) {
        std::cout << "WARNING: qc exponent mismatch\n";
        std::cout << "  Expected: " << expected_qc_exp << "\n";
        std::cout << "  Actual: " << runtime.autoconversion_qc_exponent << "\n";
    } else {
        std::cout << "  qc exponent: " << runtime.autoconversion_qc_exponent << "\n";
    }
    
    // Validate Nc exponent
    Real nc_exp_error = std::abs(runtime.autoconversion_nc_exponent - expected_nc_exp) 
                       / expected_nc_exp;
    if (nc_exp_error > tolerance) {
        std::cout << "WARNING: Nc exponent mismatch\n";
        std::cout << "  Expected: " << expected_nc_exp << "\n";
        std::cout << "  Actual: " << runtime.autoconversion_nc_exponent << "\n";
    } else {
        std::cout << "  Nc exponent: " << runtime.autoconversion_nc_exponent << "\n";
    }
    
    // Validate characteristic radius
    Real radius_error = std::abs(runtime.autoconversion_radius - expected_radius) 
                       / expected_radius;
    if (radius_error > tolerance) {
        std::cout << "INFO: Characteristic radius differs from documented value\n";
        std::cout << "  Documentation: " << expected_radius*1e6 << " um\n";
        std::cout << "  Implementation: " << runtime.autoconversion_radius*1e6 << " um\n";
    } else {
        std::cout << "  Characteristic radius: " << runtime.autoconversion_radius*1e6 << " um\n";
    }
    
    std::cout << "===========================\n\n";
  }

  void run_physics() {
    validate_autoconversion_parameters();

    // =========================================================================
    // Test Tolerance Definitions
    // =========================================================================
    // Identity tolerances: precision-dependent (mathematical exactness)
    // These verify that the implementation follows its documented formulas.
    const Real identity_tol = 10 * std::numeric_limits<Real>::epsilon();

    // Physics tolerances: precision-independent (empirical/configuration)
    // These verify physical approximations and configuration parameters.
    const Real physics_tol = 0.01;      // 1% (matches KK2000 accuracy, config params)

    // Regime thresholds: physical relevance cutoffs (DO NOT MODIFY)
    const Real small_droplet_regime_floor = 1e-15;  // Physically negligible rate for r < 1μm
    const Real absolute_floor = 1e-30;         // Numerical zero proxy
    const Real detect_threshold = 1e-3;        // Feature detection (0.1%)

    // =========================================================================
    // P3 Autoconversion Physics Property Testing Plan
    // =========================================================================
    // This suite validates that the `cloud_water_autoconversion` implementation
    // adheres to the fundamental physical principles of the Khairoutdinov and 
    // Kogan (2000) parameterization.
    //
    // Strategy:
    // Uses a deterministic grid-sampling strategy that sweeps the parameter 
    // space logarithmically to cover regimes ranging from sub-threshold 
    // (qc < 1e-8) through incipient cloud to heavy precipitation, and 
    // pristine to polluted conditions.
    // =========================================================================

    // Grid Setup: qc from 5e-9 to 1e-2 kg/kg
    // Lower bound: sub-threshold regime (qc < 1e-8, autoconversion inactive)
    // Upper bound: heavy precipitation (qc ~ 1e-2, vigorous autoconversion)
    const int n_qc = 40;
    const int n_nc = 40;
    const int num_cases = n_qc * n_nc;

    const Real log_qc_min = std::log10(5.0e-9);   
    const Real log_qc_max = std::log10(1.0e-2);
    const Real log_nc_min = std::log10(1.0e6);
    const Real log_nc_max = std::log10(1.0e9);

    view_1d<Real> qc_dev("qc_dev", num_cases);
    view_1d<Real> nc_dev("nc_dev", num_cases);
    view_1d<Real> rho_dev("rho_dev", num_cases);
    view_1d<Real> inv_qc_relvar_dev("inv_qc_relvar_dev", num_cases);
    
    view_1d<Real> qc2qr_dev("qc2qr_dev", num_cases);
    view_1d<Real> nc2nr_dev("nc2nr_dev", num_cases);
    view_1d<Real> ncautr_dev("ncautr_dev", num_cases);

    auto qc_host = Kokkos::create_mirror_view(qc_dev);
    auto nc_host = Kokkos::create_mirror_view(nc_dev);
    auto rho_host = Kokkos::create_mirror_view(rho_dev);
    auto inv_qc_relvar_host = Kokkos::create_mirror_view(inv_qc_relvar_dev);

    for (int j = 0; j < n_nc; ++j) {
      for (int i = 0; i < n_qc; ++i) {
        int idx = j * n_qc + i;
        Real alpha_qc = (Real)i / (Real)(n_qc - 1);
        Real alpha_nc = (Real)j / (Real)(n_nc - 1);

        Real qc = std::pow(10.0, log_qc_min + alpha_qc * (log_qc_max - log_qc_min));
        Real nc = std::pow(10.0, log_nc_min + alpha_nc * (log_nc_max - log_nc_min));

        qc_host(idx) = qc;
        nc_host(idx) = nc;
        rho_host(idx) = 1.0; 
        inv_qc_relvar_host(idx) = 1.0;
      }
    }

    Kokkos::deep_copy(qc_dev, qc_host);
    Kokkos::deep_copy(nc_dev, nc_host);
    Kokkos::deep_copy(rho_dev, rho_host);
    Kokkos::deep_copy(inv_qc_relvar_dev, inv_qc_relvar_host);

    const int pack_size = Spack::n;
    const int num_packs = (num_cases + pack_size - 1) / pack_size;

    Kokkos::parallel_for("P3_Property_Test_Kernel", num_packs, KOKKOS_LAMBDA(const Int& ip) {
      const int offset = ip * pack_size;
      
      Spack rho, qc, nc, inv_qc_relvar;
      Spack qc2qr, nc2nr, ncautr;
      
      bool any_valid = false;
      for (int s = 0; s < pack_size; ++s) {
        int idx = offset + s;
        if (idx < num_cases) {
          rho[s] = rho_dev(idx);
          qc[s] = qc_dev(idx);
          nc[s] = nc_dev(idx);
          inv_qc_relvar[s] = inv_qc_relvar_dev(idx);
          any_valid = true;
        } else {
          rho[s] = 1.0;
          qc[s] = 0.0;
          nc[s] = 1.0e6;
          inv_qc_relvar[s] = 1.0;
        }
      }

      if (any_valid) {
          Functions::cloud_water_autoconversion(
              rho, qc, nc, inv_qc_relvar, 
              qc2qr, nc2nr, ncautr,
              p3::Functions<Real,DefaultDevice>::P3Runtime()
          );
      }
      
      for (int s = 0; s < pack_size; ++s) {
        int idx = offset + s;
        if (idx < num_cases) {
          qc2qr_dev(idx) = qc2qr[s];
          nc2nr_dev(idx) = nc2nr[s];
          ncautr_dev(idx) = ncautr[s];
        }
      }
    });

    Kokkos::fence();
    
    auto qc2qr_host = Kokkos::create_mirror_view(qc2qr_dev);
    auto nc2nr_host = Kokkos::create_mirror_view(nc2nr_dev);
    auto ncautr_host = Kokkos::create_mirror_view(ncautr_dev);

    Kokkos::deep_copy(qc2qr_host, qc2qr_dev);
    Kokkos::deep_copy(nc2nr_host, nc2nr_dev);
    Kokkos::deep_copy(ncautr_host, ncautr_dev);

    int failures = 0;

    for (int j = 0; j < n_nc; ++j) {
      for (int i = 0; i < n_qc; ++i) {
        int idx = j * n_qc + i;
        Real qc = qc_host(idx);
        Real nc = nc_host(idx);
        Real R = qc2qr_host(idx);
        Real N_loss = nc2nr_host(idx);
        Real R_ncautr = ncautr_host(idx); 

        // Check for non-zero below threshold
        if (qc < 1e-8) {
           if (R > 1e-30) {
              std::cout << "FAIL: Non-zero rate below threshold (1e-8). qc=" << qc << " R=" << R << "\n";
              failures++;
           }
           continue;
        }

        if (R < 1e-30) R = 0.0;
        if (N_loss < 1e-30) N_loss = 0.0;
        if (R_ncautr < 1e-30) R_ncautr = 0.0;
        
        // =====================================================================
        // 1. Physical Sensitivity Tests (Monotonicity)
        // =====================================================================
        // Uses strict inequality for Nc (physics requirement) and 0.1% tolerance for qc
        // (sensitivity detection). Both reference the global absolute_floor (1e-30).
        
        const Real relative_tolerance = 1e-3;  // 0.1% sensitivity threshold for qc monotonicity

        if (j < n_nc - 1) {
            int idx_next = (j + 1) * n_qc + i;
            Real R_next = qc2qr_host(idx_next);
            if (R_next < 1e-30) R_next = 0.0;
            
            // A. Colloidal Stability (Inverse Nc Dependency): dR/dNc < 0
            // R_next (high Nc) must be strictly less than R (low Nc).
            if (R > absolute_floor) {
                // FAIL if R_next is greater than or equal to R.
                // We do not use the loose 'relative_tolerance' here because the physics 
                // requires strict monotonicity.
                if (R_next >= R) {
                     std::cout << "Monotonicity Nc Fail: R increased or stayed same with Nc. qc=" << qc 
                               << " R=" << R << " R_next=" << R_next << "\n";
                     failures++;
                }
            }
        }
        
        if (i < n_qc - 1) {
            int idx_next = j * n_qc + (i + 1);
            Real R_next = qc2qr_host(idx_next);
             if (R_next < 1e-30) R_next = 0.0;
            
            // B. Water Content Sensitivity (Positive qc Dependency): dR/dqc > 0
            // R_next (high qc) should be greater than R (low qc).
            if (R_next > absolute_floor) {
                // If R_next < R, it's a failure. 
                Real min_allowed_R_next = R * (1.0 - relative_tolerance);
                if (R_next < min_allowed_R_next) {
                     std::cout << "Monotonicity qc Fail: R decreased with qc. qc=" << qc 
                               << " R=" << R << " R_next=" << R_next << "\n";
                     failures++;
                }
            }
        }
        
        // =====================================================================
        // 2. Consistency Tests (Constraint Satisfaction)
        // =====================================================================

        if (R > absolute_floor) {
             // A. Specific Loss Conservation (Mathematical Identity)
             // Verifies: nc2nr_autoconv_tend = qc2qr_autoconv_tend * (nc/qc)
             Real expected_N_loss = R * nc / qc;
             Real rel_error = std::abs(N_loss - expected_N_loss) / std::max(absolute_floor, expected_N_loss);
             
             if (rel_error > identity_tol) {
                 std::cout << "Specific Loss Conservation Fail (should be near machine precision):\n";
                 std::cout << "  Expected: " << expected_N_loss << "\n";
                 std::cout << "  Actual:   " << N_loss << "\n";
                 std::cout << "  Rel Error: " << rel_error << " (tolerance: " << identity_tol << ")\n";
                 failures++;
             }
        }

        // B. Rain Embryo Size Consistency
        if (R > 1e-20 && R_ncautr > 1e-20) {
             // B. Rain Embryo Size Consistency (Configuration Parameter)
             // Verifies: autoconversion_radius parameter is preserved through calculation
             Real mass_drop = R / R_ncautr;
             
             using C = scream::physics::Constants<Real>;
             const Real r_expected = 25.0e-6;  // Must match autoconversion_radius
             const Real expected_mass = (4.0/3.0) * C::Pi * 1000.0 * std::pow(r_expected, 3);
             
             if (std::abs(mass_drop - expected_mass) / expected_mass > physics_tol) {
                  std::cout << "Embryo Size Fail (configuration mismatch):\n";
                  std::cout << "  Computed mass: " << mass_drop << " kg\n";
                  std::cout << "  Expected mass: " << expected_mass << " kg (r=25um)\n";
                  std::cout << "  Tolerance: " << physics_tol*100 << "%\n";
                  failures++;
             }
        }
        
        // =====================================================================
        // 3. Limit & Singularity Tests
        // =====================================================================
        
        Real mean_mass = qc / nc;
        // Correct formula: r = (3*m / (4*pi*rho))^(1/3)
        // rho_w = 1000.
        using C = scream::physics::Constants<Real>;
        Real mean_rad = std::pow((3.0 * mean_mass) / (4.0 * C::Pi * 1000.0), 1.0/3.0);
        
        // A. Small Droplet Regime Limit
        if (mean_rad < 1e-6) {
            // A. Small Droplet Regime Limit (Physical Relevance Check)
            // For very small activated droplets (r < 1 μm), autoconversion rates 
            // should be physically negligible due to low collision efficiency.
            // This regime represents incipient cloud or sub-threshold conditions.
            // Note: Implementation threshold is qc < 1e-8, not radius-based.
            if (R > small_droplet_regime_floor) {
                std::cout << "Small Droplet Regime Fail: r=" << mean_rad << " m, R=" << R << " kg/kg/s\n";
                std::cout << "  (Rate exceeds physical relevance threshold: " << small_droplet_regime_floor << ")\n";
                failures++;
            }
        }
      }
    }

    // =========================================================================
    // CHECK 4: Subgrid Variance Scaling
    // =========================================================================
    
    view_1d<Real> var_rate_1("rate_1", n_qc);
    view_1d<Real> var_rate_2("rate_2", n_qc); 
    
    Real fixed_nc = 100.0e6; 

    Kokkos::parallel_for("Variance_Check", n_qc, KOKKOS_LAMBDA(const int& i) {
        Real alpha = (Real)i / (Real)(n_qc - 1);
        Real qc = std::pow(10.0, -6.0 + alpha * 4.0);
        
        Spack s_rho(1.0), s_qc(qc), s_nc(fixed_nc);
        Spack s_rate_1, s_rate_2, dummy;
        
        // Run 1: Homogeneous
        Functions::cloud_water_autoconversion(
            s_rho, s_qc, s_nc, Spack(1.0), 
            s_rate_1, dummy, dummy, 
            p3::Functions<Real,DefaultDevice>::P3Runtime()
        );
        
        // Run 2: High Variance
        Functions::cloud_water_autoconversion(
            s_rho, s_qc, s_nc, Spack(0.1), 
            s_rate_2, dummy, dummy, 
            p3::Functions<Real,DefaultDevice>::P3Runtime()
        );
        
        var_rate_1(i) = s_rate_1[0];
        var_rate_2(i) = s_rate_2[0];
    });
    
    Kokkos::fence();
    auto h_rate_1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), var_rate_1);
    auto h_rate_2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), var_rate_2);

    bool is_variance_scaling_active = false;
    int variance_suppression_failures = 0;
    const Real variance_detect_eps = detect_threshold;  // 0.1% change indicates active scaling
    const Real variance_validate_eps = 1e-12;           // Numerical noise floor for suppression check

    for(int i=0; i<n_qc; ++i) {
        if (h_rate_1(i) > 1e-30) {
            Real ratio = h_rate_2(i) / h_rate_1(i);

            // Detection: Do standard and variance rates differ?
            if (std::abs(ratio - 1.0) > variance_detect_eps) {
                is_variance_scaling_active = true;
                
                // Validation: If active, is the rate suppressed? (Physical violation)
                // Jensen's inequality for convex functions implies enhancement (ratio > 1).
                // We fail if ratio < 1.0 (minus tolerance).
                if (ratio < 1.0 - variance_validate_eps) {
                     variance_suppression_failures++;
                }
            }
        }
    }

    if (is_variance_scaling_active) {
        if (variance_suppression_failures > 0) {
             std::cout << "Variance Scaling FAIL: " << variance_suppression_failures << " cases suppressed rate (physical violation).\n";
             failures += variance_suppression_failures;
        } else {
             std::cout << "Variance Scaling: ACTIVE and working correctly (rates enhanced).\n";
        }
    } else {
         std::cout << "Variance Scaling: DISABLED (sgs_var_coef approx 1)\n";
    }

    REQUIRE(failures == 0);
  }

}; //  TestP3CloudWaterAutoconversion

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_cloud_water_autoconversion_test", "[p3_cloud_water_autoconversion_test]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion;

  T t;
  t.run_physics();
  t.run_bfb();
}

} // namespace
