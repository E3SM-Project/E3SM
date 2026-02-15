#include "catch2/catch.hpp"

#include "p3_functions.hpp"
#include "p3_test_data.hpp"
#include "p3_unit_tests_common.hpp"

#include "share/core/eamxx_types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace scream {
namespace p3 {
namespace unit_test {

namespace {

constexpr int kFracSweepSize = 15;
constexpr Real kGlaciatedFloor = 1e-4;
constexpr Real kAbsoluteFloor = 1e-30;

// Physical region tags used by back_to_cell_average scaling:
// L/R/I are single-phase cloud fractions; IL/IR/LR are pairwise intersections;
// GLACIATED_OR_I depends on runtime option use_separate_ice_liq_frac.
enum class ScaleKind {
  L,               // liquid cloud region
  R,               // rain cloud region
  I,               // ice cloud region
  IL,              // ice-liquid overlap
  IR,              // ice-rain overlap
  LR,              // liquid-rain overlap
  GLACIATED_OR_I,  // glaciated-only (or full ice when runtime option is off)
  UNCHANGED        // already cell-averaged; no additional scaling
};

// Single source of truth: each tendency variable and the cloud-fraction factor
// it should use. Property checks, kernel IO, and expected-value logic are all
// generated from this table to avoid mapping drift.
// To add a new tendency, add it here with its ScaleKind and, if serialized for
// BFB, add it to BACK_TO_CELL_AVERAGE_BFB_COMPARE_LIST as well.
#define BACK_TO_CELL_AVERAGE_TENDENCY_LIST(X)                                                   \
  X(qc2qr_accret_tend,        ScaleKind::LR)                                                    \
  X(qr2qv_evap_tend,          ScaleKind::R)                                                     \
  X(qc2qr_autoconv_tend,      ScaleKind::L)                                                     \
  X(nc_accret_tend,           ScaleKind::LR)                                                    \
  X(nc_selfcollect_tend,      ScaleKind::L)                                                     \
  X(nc2nr_autoconv_tend,      ScaleKind::L)                                                     \
  X(nr_selfcollect_tend,      ScaleKind::R)                                                     \
  X(nr_evap_tend,             ScaleKind::R)                                                     \
  X(ncautr,                   ScaleKind::LR)                                                    \
  X(qi2qv_sublim_tend,        ScaleKind::GLACIATED_OR_I)                                        \
  X(nr_ice_shed_tend,         ScaleKind::IL)                                                    \
  X(qc2qi_hetero_freeze_tend, ScaleKind::IL)                                                    \
  X(qr2qi_collect_tend,       ScaleKind::IR)                                                    \
  X(qc2qr_ice_shed_tend,      ScaleKind::IL)                                                    \
  X(qi2qr_melt_tend,          ScaleKind::I)                                                     \
  X(qc2qi_collect_tend,       ScaleKind::IL)                                                    \
  X(qr2qi_immers_freeze_tend, ScaleKind::R)                                                     \
  X(ni2nr_melt_tend,          ScaleKind::I)                                                     \
  X(nc_collect_tend,          ScaleKind::IL)                                                    \
  X(ncshdc,                   ScaleKind::IL)                                                    \
  X(nc2ni_immers_freeze_tend, ScaleKind::L)                                                     \
  X(nr_collect_tend,          ScaleKind::IR)                                                    \
  X(ni_selfcollect_tend,      ScaleKind::I)                                                     \
  X(qv2qi_vapdep_tend,        ScaleKind::GLACIATED_OR_I)                                        \
  X(nr2ni_immers_freeze_tend, ScaleKind::R)                                                     \
  X(ni_sublim_tend,           ScaleKind::GLACIATED_OR_I)                                        \
  X(qv2qi_nucleat_tend,       ScaleKind::UNCHANGED)                                             \
  X(ni_nucleat_tend,          ScaleKind::UNCHANGED)                                             \
  X(qc2qi_berg_tend,          ScaleKind::IL)                                                    \
  X(ncheti_cnt,               ScaleKind::L)                                                     \
  X(qcheti_cnt,               ScaleKind::L)                                                     \
  X(nicnt,                    ScaleKind::L)                                                     \
  X(qicnt,                    ScaleKind::L)                                                     \
  X(ninuc_cnt,                ScaleKind::L)                                                     \
  X(qinuc_cnt,                ScaleKind::L)

// BFB list matches the PTD_RW_SCALARS_ONLY serialized subset in
// BackToCellAverageData:
// - Includes legacy pass-through fields qcnuc and nc_nuceat_tend.
// - Excludes *_cnt fields, which are not serialized for BFB in test data.
#define BACK_TO_CELL_AVERAGE_BFB_COMPARE_LIST(X)                                                \
  X(qc2qr_accret_tend)                                                                           \
  X(qr2qv_evap_tend)                                                                             \
  X(qc2qr_autoconv_tend)                                                                         \
  X(nc_accret_tend)                                                                              \
  X(nc_selfcollect_tend)                                                                         \
  X(nc2nr_autoconv_tend)                                                                         \
  X(nr_selfcollect_tend)                                                                         \
  X(nr_evap_tend)                                                                                \
  X(ncautr)                                                                                      \
  X(qcnuc)                                                                                       \
  X(nc_nuceat_tend)                                                                              \
  X(qi2qv_sublim_tend)                                                                           \
  X(nr_ice_shed_tend)                                                                            \
  X(qc2qi_hetero_freeze_tend)                                                                    \
  X(qr2qi_collect_tend)                                                                          \
  X(qc2qr_ice_shed_tend)                                                                         \
  X(qi2qr_melt_tend)                                                                             \
  X(qc2qi_collect_tend)                                                                          \
  X(qr2qi_immers_freeze_tend)                                                                    \
  X(ni2nr_melt_tend)                                                                             \
  X(nc_collect_tend)                                                                             \
  X(ncshdc)                                                                                      \
  X(nc2ni_immers_freeze_tend)                                                                    \
  X(nr_collect_tend)                                                                             \
  X(ni_selfcollect_tend)                                                                         \
  X(qv2qi_vapdep_tend)                                                                           \
  X(nr2ni_immers_freeze_tend)                                                                    \
  X(ni_sublim_tend)                                                                              \
  X(qv2qi_nucleat_tend)                                                                          \
  X(ni_nucleat_tend)                                                                             \
  X(qc2qi_berg_tend)

} // anonymous namespace

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3BackToCellAverage : public UnitWrap::UnitTest<D>::Base {

  // One synthetic grid-point scenario for host-side validation.
  // Stores cloud fractions, derived overlaps, and pre/post tendencies.
  struct BackToCellAverageCase {
    Real cld_frac_l, cld_frac_r, cld_frac_i;
    Real ir_cldm, il_cldm, lr_cldm, cld_frac_glaciated;
    bool context;
#define DECL_CASE_FIELD(name, scale) Real name##_input, name##_output;
    BACK_TO_CELL_AVERAGE_TENDENCY_LIST(DECL_CASE_FIELD)
#undef DECL_CASE_FIELD
  };

  struct BackToCellAverageTestData {
    int n_frac;
    int num_cases;
    bool use_separate_ice_liq_frac;
    Real identity_tol;
    Real absolute_floor;
    typename view_1d<BackToCellAverageCase>::HostMirror cases;
  };

  using Runtime = typename Functions::P3Runtime;

  // Returns the physically expected scaling factor for each tendency.
  // GLACIATED_OR_I switches between glaciated-only and full-ice regions.
  KOKKOS_INLINE_FUNCTION
  static Real cloud_factor(const BackToCellAverageCase& c, const ScaleKind scale,
                           const bool use_separate_ice_liq_frac) {
    switch (scale) {
    case ScaleKind::L: return c.cld_frac_l;
    case ScaleKind::R: return c.cld_frac_r;
    case ScaleKind::I: return c.cld_frac_i;
    case ScaleKind::IL: return c.il_cldm;
    case ScaleKind::IR: return c.ir_cldm;
    case ScaleKind::LR: return c.lr_cldm;
    case ScaleKind::GLACIATED_OR_I:
      return use_separate_ice_liq_frac ? c.cld_frac_glaciated : c.cld_frac_i;
    case ScaleKind::UNCHANGED: return 1;
    }
    return 1;
  }

  KOKKOS_INLINE_FUNCTION
  static bool is_runtime_sensitive(const ScaleKind scale) {
    return scale == ScaleKind::GLACIATED_OR_I;
  }

  KOKKOS_INLINE_FUNCTION
  static bool is_scaled(const ScaleKind scale) {
    return scale != ScaleKind::UNCHANGED;
  }

  // Deterministic non-negative inputs with distinct magnitudes to ensure each
  // tendency mapping is independently testable.
  void set_case_inputs(BackToCellAverageCase& c, const Real base) const {
    int idx = 1;
#define INIT_INPUT(name, scale)                                                                  \
    c.name##_input = base * (1 + 0.05 * idx);                                                   \
    c.name##_output = -1;                                                                        \
    ++idx;
    BACK_TO_CELL_AVERAGE_TENDENCY_LIST(INIT_INPUT)
#undef INIT_INPUT
  }

  // Applies alternating signs so mapping checks exercise both source and sink
  // tendencies. Multiplicative scaling should be sign-agnostic.
  void apply_alternating_sign_pattern(BackToCellAverageCase& c, const int case_idx) const {
    int tidx = 0;
#define APPLY_SIGN(name, scale)                                                                  \
    if (((case_idx + tidx) % 2) == 1) {                                                         \
      c.name##_input = -c.name##_input;                                                         \
    }                                                                                            \
    ++tidx;
    BACK_TO_CELL_AVERAGE_TENDENCY_LIST(APPLY_SIGN)
#undef APPLY_SIGN
  }

  // Executes back_to_cell_average for a batch of synthetic cases and stores
  // both derived overlap fractions and output tendencies for host-side checks.
  void run_back_to_cell_average_kernel(const view_1d<BackToCellAverageCase>& cases_dev,
                                       const bool use_separate_ice_liq_frac) const {
    Runtime runtime_options;
    runtime_options.use_separate_ice_liq_frac = use_separate_ice_liq_frac;

    const Int num_cases = cases_dev.extent_int(0);
    const Int num_packs = (num_cases + Spack::n - 1) / Spack::n;

    Kokkos::parallel_for("p3_back_to_cell_average_sweep", num_packs,
                         KOKKOS_LAMBDA(const Int& ipack) {
      const Int offset = ipack * Spack::n;

      Spack cld_frac_l, cld_frac_r, cld_frac_i;
#define DECL_PACK(name, scale) Spack name;
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(DECL_PACK)
#undef DECL_PACK

      Smask context;

      for (Int s = 0; s < Spack::n; ++s) {
        const Int idx = offset + s;
        if (idx < num_cases) {
          const auto& d = cases_dev(idx);
          cld_frac_l[s] = d.cld_frac_l;
          cld_frac_r[s] = d.cld_frac_r;
          cld_frac_i[s] = d.cld_frac_i;
#define LOAD_INPUT(name, scale) name[s] = d.name##_input;
          BACK_TO_CELL_AVERAGE_TENDENCY_LIST(LOAD_INPUT)
#undef LOAD_INPUT
          context.set(s, d.context);
        } else {
          cld_frac_l[s] = 0;
          cld_frac_r[s] = 0;
          cld_frac_i[s] = 0;
#define LOAD_DEFAULT(name, scale) name[s] = 0;
          BACK_TO_CELL_AVERAGE_TENDENCY_LIST(LOAD_DEFAULT)
#undef LOAD_DEFAULT
          context.set(s, false);
        }
      }

      const Spack ir_cldm = min(cld_frac_i, cld_frac_r);
      const Spack il_cldm = min(cld_frac_i, cld_frac_l);
      const Spack lr_cldm = min(cld_frac_l, cld_frac_r);
      const Spack cld_frac_glaciated = max(kGlaciatedFloor, cld_frac_i - il_cldm);

      Functions::back_to_cell_average(
        cld_frac_l, cld_frac_r, cld_frac_i,
        qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend,
        nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend,
        nr_selfcollect_tend, nr_evap_tend, ncautr,
        qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend,
        qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend,
        qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend,
        nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend,
        nr_collect_tend, ni_selfcollect_tend, qv2qi_vapdep_tend,
        nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend,
        ni_nucleat_tend, qc2qi_berg_tend,
        ncheti_cnt, qcheti_cnt, nicnt, qicnt, ninuc_cnt, qinuc_cnt,
        context, runtime_options);

      for (Int s = 0; s < Spack::n; ++s) {
        const Int idx = offset + s;
        if (idx < num_cases) {
          auto& d = cases_dev(idx);
          d.ir_cldm = ir_cldm[s];
          d.il_cldm = il_cldm[s];
          d.lr_cldm = lr_cldm[s];
          d.cld_frac_glaciated = cld_frac_glaciated[s];
#define STORE_OUTPUT(name, scale) d.name##_output = name[s];
          BACK_TO_CELL_AVERAGE_TENDENCY_LIST(STORE_OUTPUT)
#undef STORE_OUTPUT
        }
      }
    });

    Kokkos::fence();
  }

  // Dense coverage of cloud-fraction phase space:
  // 15x15x15 combinations spanning the full [0,1] range for fi/fl/fr.
  BackToCellAverageTestData setup_parameter_sweep(const bool use_separate_ice_liq_frac = false,
                                                  const bool mixed_sign_inputs = false) {
    const int n_frac = kFracSweepSize;
    const int num_cases = n_frac * n_frac * n_frac;
    const Real identity_tol = 10 * std::numeric_limits<Real>::epsilon();
    const Real absolute_floor = kAbsoluteFloor;

    view_1d<BackToCellAverageCase> cases_dev("back_to_cell_average_cases", num_cases);
    auto cases_host = Kokkos::create_mirror_view(cases_dev);

    int idx = 0;
    for (int ii = 0; ii < n_frac; ++ii) {
      for (int il = 0; il < n_frac; ++il) {
        for (int ir = 0; ir < n_frac; ++ir) {
          auto& c = cases_host(idx);
          c.cld_frac_i = ii / (n_frac - 1.0);
          c.cld_frac_l = il / (n_frac - 1.0);
          c.cld_frac_r = ir / (n_frac - 1.0);
          c.context = true;

          const Real base = 1e-6 * (1 + 0.02 * (idx % 23));
          set_case_inputs(c, base);
          if (mixed_sign_inputs) {
            apply_alternating_sign_pattern(c, idx);
          }
          ++idx;
        }
      }
    }

    Kokkos::deep_copy(cases_dev, cases_host);
    run_back_to_cell_average_kernel(cases_dev, use_separate_ice_liq_frac);
    Kokkos::deep_copy(cases_host, cases_dev);

    std::cout << "\n=== P3 back_to_cell_average sweep ===\n"
              << "  n_frac: " << n_frac << "\n"
              << "  total cases: " << num_cases << "\n"
              << "  use_separate_ice_liq_frac: " << (use_separate_ice_liq_frac ? "true" : "false") << "\n"
              << "  identity_tol: " << identity_tol << "\n"
              << "  absolute_floor: " << absolute_floor << "\n"
              << "=====================================\n";

    return BackToCellAverageTestData{
      n_frac,
      num_cases,
      use_separate_ice_liq_frac,
      identity_tol,
      absolute_floor,
      cases_host
    };
  }

  // Deterministic BFB fixture with representative overlap regimes so baseline
  // behavior is stable and physically interpretable.
  void initialize_deterministic_bfb_inputs(std::array<BackToCellAverageData, max_pack_size>& data) const {
    const std::array<Real, max_pack_size> frac_l = {
      0.0, 1.0, 0.6, 0.3,
      0.9, 0.1, 0.4, 0.7,
      0.0, 1.0, 0.5, 0.2,
      0.8, 0.25, 0.65, 0.35
    };
    const std::array<Real, max_pack_size> frac_r = {
      0.0, 1.0, 0.2, 0.6,
      0.1, 0.9, 0.8, 0.4,
      1.0, 0.0, 0.5, 0.7,
      0.3, 0.65, 0.25, 0.95
    };
    const std::array<Real, max_pack_size> frac_i = {
      0.0, 1.0, 0.8, 0.1,
      0.5, 0.2, 0.9, 0.4,
      0.3, 0.0, 0.5, 1.0,
      0.75, 0.55, 0.15, 0.6
    };

    for (Int i = 0; i < max_pack_size; ++i) {
      auto& d = data[i];
      d.cld_frac_l = frac_l[i];
      d.cld_frac_r = frac_r[i];
      d.cld_frac_i = frac_i[i];
      d.use_hetfrz_classnuc = false;
      d.context = true;

      const Real base = 1e-6 * (i + 1);
      int idx = 1;
#define INIT_BFB_INPUT(name, scale) d.name = base * (1 + 0.03 * idx++);
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(INIT_BFB_INPUT)
#undef INIT_BFB_INPUT

      d.qcnuc = 0;
      d.nc_nuceat_tend = 0;
    }
  }

  // Bit-for-bit regression check for the baseline test infrastructure.
  // Kept separate from property tests that validate physics bookkeeping logic.
  void run_bfb() {
    static_assert(max_pack_size >= 16,
                  "Deterministic BFB fixture assumes max_pack_size >= 16.");

    std::array<BackToCellAverageData, max_pack_size> baseline_data;
    initialize_deterministic_bfb_inputs(baseline_data);

    view_1d<BackToCellAverageData> device_data("back_to_cell_average_bfb", max_pack_size);
    auto host_data = Kokkos::create_mirror_view(device_data);

    std::copy(baseline_data.begin(), baseline_data.end(), host_data.data());
    Kokkos::deep_copy(device_data, host_data);

    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        baseline_data[i].read(Base::m_ifile);
      }
    }

    Runtime runtime_options;
    runtime_options.use_separate_ice_liq_frac = false;

    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& ipack) {
      const Int offset = ipack * Spack::n;

      Spack cld_frac_l, cld_frac_r, cld_frac_i;
#define DECL_PACK(name, scale) Spack name;
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(DECL_PACK)
#undef DECL_PACK
      Smask context;

      for (Int s = 0; s < Spack::n; ++s) {
        const Int vs = offset + s;
        const auto& d = device_data(vs);
        cld_frac_l[s] = d.cld_frac_l;
        cld_frac_r[s] = d.cld_frac_r;
        cld_frac_i[s] = d.cld_frac_i;
#define LOAD_BFB_INPUT(name, scale) name[s] = d.name;
        BACK_TO_CELL_AVERAGE_TENDENCY_LIST(LOAD_BFB_INPUT)
#undef LOAD_BFB_INPUT
        context.set(s, d.context);
      }

      Functions::back_to_cell_average(
        cld_frac_l, cld_frac_r, cld_frac_i,
        qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend,
        nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend,
        nr_selfcollect_tend, nr_evap_tend, ncautr,
        qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend,
        qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend,
        qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend,
        nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend,
        nr_collect_tend, ni_selfcollect_tend, qv2qi_vapdep_tend,
        nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend,
        ni_nucleat_tend, qc2qi_berg_tend,
        ncheti_cnt, qcheti_cnt, nicnt, qicnt, ninuc_cnt, qinuc_cnt,
        context, runtime_options);

      for (Int s = 0; s < Spack::n; ++s) {
        const Int vs = offset + s;
        auto& d = device_data(vs);
#define STORE_BFB_OUTPUT(name, scale) d.name = name[s];
        BACK_TO_CELL_AVERAGE_TENDENCY_LIST(STORE_BFB_OUTPUT)
#undef STORE_BFB_OUTPUT
      }
    });

    Kokkos::deep_copy(host_data, device_data);

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
#define REQUIRE_BFB(name) REQUIRE(baseline_data[s].name == host_data[s].name);
        BACK_TO_CELL_AVERAGE_BFB_COMPARE_LIST(REQUIRE_BFB)
#undef REQUIRE_BFB
      }
    } else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        host_data(s).write(Base::m_ofile);
      }
    }
  }

  // Primary bookkeeping check: each tendency must use the cloud-fraction factor
  // defined by its process location in the implementation.
  void run_process_location_checks(const BackToCellAverageTestData& data) const {
    int failures = 0;

    for (int i = 0; i < data.num_cases; ++i) {
      const auto& c = data.cases(i);
#define CHECK_MAPPING(name, scale)                                                               \
      {                                                                                          \
        const Real expected = c.name##_input * cloud_factor(c, scale, data.use_separate_ice_liq_frac); \
        const Real actual = c.name##_output;                                                     \
        const Real rel_error = std::abs(actual - expected) / std::max(data.absolute_floor, std::abs(expected)); \
        if (rel_error > data.identity_tol) {                                                     \
          std::cout << "Process-location FAIL: " #name                                          \
                    << " case=" << i                                                            \
                    << " expected=" << expected                                                 \
                    << " actual=" << actual                                                     \
                    << " rel_error=" << rel_error                                               \
                    << " tol=" << data.identity_tol << "\n";                                 \
          ++failures;                                                                            \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_MAPPING)
#undef CHECK_MAPPING
    }

    REQUIRE(failures == 0);
  }

  // Physical and mathematical bounds:
  // non-negative inputs remain non-negative, and cell-avg tendencies do not
  // exceed in-cloud tendencies for scaled processes.
  void run_bounds_and_scaling_checks(const BackToCellAverageTestData& data) const {
    int failures = 0;

    for (int i = 0; i < data.num_cases; ++i) {
      const auto& c = data.cases(i);
#define CHECK_BOUNDS(name, scale)                                                                \
      {                                                                                          \
        const Real input = c.name##_input;                                                       \
        const Real output = c.name##_output;                                                     \
        if (is_scaled(scale)) {                                                                  \
          if (input >= 0 && output < -data.absolute_floor) {                                    \
            std::cout << "Non-negativity FAIL: " #name                                          \
                      << " case=" << i << " input=" << input                                   \
                      << " output=" << output << "\n";                                        \
            ++failures;                                                                          \
          }                                                                                      \
          const Real upper = input + data.identity_tol * std::max<Real>(1, std::abs(input));    \
          if (output > upper) {                                                                  \
            std::cout << "Scaling reduction FAIL: " #name                                        \
                      << " case=" << i << " input=" << input                                   \
                      << " output=" << output << " upper=" << upper << "\n";               \
            ++failures;                                                                          \
          }                                                                                      \
        } else {                                                                                 \
          const Real diff = std::abs(output - input);                                            \
          const Real tol = data.identity_tol * std::max<Real>(1, std::abs(input));              \
          if (diff > tol) {                                                                      \
            std::cout << "Unchanged-term identity FAIL: " #name                                 \
                      << " case=" << i << " input=" << input                                   \
                      << " output=" << output << " diff=" << diff                              \
                      << " tol=" << tol << "\n";                                              \
            ++failures;                                                                          \
          }                                                                                      \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_BOUNDS)
#undef CHECK_BOUNDS
    }

    REQUIRE(failures == 0);
  }

  // Validates overlap geometry (intersection <= constituents) and enforces the
  // glaciated minimum used for numerical robustness.
  void run_cloud_fraction_logic_checks(const BackToCellAverageTestData& data) const {
    int failures = 0;

    for (int i = 0; i < data.num_cases; ++i) {
      const auto& c = data.cases(i);
      const Real ir_max = std::min(c.cld_frac_i, c.cld_frac_r) + data.identity_tol;
      const Real il_max = std::min(c.cld_frac_i, c.cld_frac_l) + data.identity_tol;
      const Real lr_max = std::min(c.cld_frac_l, c.cld_frac_r) + data.identity_tol;

      if (c.ir_cldm > ir_max) {
        std::cout << "Intersection FAIL (ir): case=" << i
                  << " ir_cldm=" << c.ir_cldm << " max=" << ir_max << "\n";
        ++failures;
      }
      if (c.il_cldm > il_max) {
        std::cout << "Intersection FAIL (il): case=" << i
                  << " il_cldm=" << c.il_cldm << " max=" << il_max << "\n";
        ++failures;
      }
      if (c.lr_cldm > lr_max) {
        std::cout << "Intersection FAIL (lr): case=" << i
                  << " lr_cldm=" << c.lr_cldm << " max=" << lr_max << "\n";
        ++failures;
      }
      if (c.cld_frac_glaciated < kGlaciatedFloor - data.identity_tol) {
        std::cout << "Glaciated fraction floor FAIL: case=" << i
                  << " cld_frac_glaciated=" << c.cld_frac_glaciated << "\n";
        ++failures;
      }
    }

    REQUIRE(failures == 0);
  }

  // Isolates runtime-option behavior: only GLACIATED_OR_I tendencies should
  // respond to use_separate_ice_liq_frac; all others remain unchanged.
  void run_runtime_option_checks() {
    const auto data_default = setup_parameter_sweep(false);
    const auto data_separate = setup_parameter_sweep(true);

    int failures = 0;

    for (int i = 0; i < data_default.num_cases; ++i) {
      const auto& d0 = data_default.cases(i);
      const auto& d1 = data_separate.cases(i);

#define CHECK_RUNTIME_DIFF(name, scale)                                                          \
      {                                                                                          \
        const Real diff = std::abs(d0.name##_output - d1.name##_output);                        \
        const Real tol = data_default.identity_tol * std::max<Real>(1, std::abs(d0.name##_output)); \
        if (is_runtime_sensitive(scale)) {                                                       \
          const Real expected0 = d0.name##_input * d0.cld_frac_i;                               \
          const Real expected1 = d1.name##_input * d1.cld_frac_glaciated;                       \
          const Real err0 = std::abs(d0.name##_output - expected0) / std::max(data_default.absolute_floor, std::abs(expected0)); \
          const Real err1 = std::abs(d1.name##_output - expected1) / std::max(data_default.absolute_floor, std::abs(expected1)); \
          if (err0 > data_default.identity_tol || err1 > data_default.identity_tol) {           \
            std::cout << "Runtime branch mapping FAIL: " #name                                  \
                      << " case=" << i                                                          \
                      << " err_default=" << err0                                                \
                      << " err_separate=" << err1 << "\n";                                   \
            ++failures;                                                                          \
          }                                                                                      \
          /* If il_cldm==0, cld_frac_glaciated==cld_frac_i, so identical outputs */             \
          /* are expected. Require sensitivity only when scaling factors differ. */              \
          const Real scaling_diff = std::abs(d0.cld_frac_i - d1.cld_frac_glaciated);            \
          if (scaling_diff > data_default.identity_tol && diff <= tol) {                         \
            std::cout << "Runtime branch sensitivity FAIL: " #name                              \
                      << " case=" << i                                                          \
                      << " scaling_diff=" << scaling_diff                                       \
                      << " actual_diff=" << diff << "\n";                                     \
            ++failures;                                                                          \
          }                                                                                      \
        } else {                                                                                 \
          if (diff > tol) {                                                                      \
            std::cout << "Runtime branch isolation FAIL: " #name                                \
                      << " case=" << i                                                          \
                      << " diff=" << diff << " tol=" << tol << "\n";                        \
            ++failures;                                                                          \
          }                                                                                      \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_RUNTIME_DIFF)
#undef CHECK_RUNTIME_DIFF
    }

    REQUIRE(failures == 0);
  }

  // Verifies .set(context, ...) masking semantics: context=false lanes must
  // retain their input tendency values even when neighboring lanes are active.
  void run_context_mask_checks() {
    const int num_cases = 2 * Spack::n;
    const Real identity_tol = 10 * std::numeric_limits<Real>::epsilon();
    const Real absolute_floor = kAbsoluteFloor;

    view_1d<BackToCellAverageCase> cases_dev("back_to_cell_average_context_mask", num_cases);
    auto cases_host = Kokkos::create_mirror_view(cases_dev);

    for (int i = 0; i < num_cases; ++i) {
      auto& c = cases_host(i);
      c.cld_frac_i = (i % 3) * 0.3;
      c.cld_frac_l = ((i + 1) % 4) * 0.25;
      c.cld_frac_r = ((i + 2) % 5) * 0.2;
      c.context = (i < Spack::n) ? ((i % 2) == 0) : ((i % 2) == 1);
      set_case_inputs(c, 2e-6 * (i + 1));
      apply_alternating_sign_pattern(c, i);
    }
    // Explicitly pin one all-zero/cloud-fraction lane with context=false.
    cases_host(0).cld_frac_i = 0;
    cases_host(0).cld_frac_l = 0;
    cases_host(0).cld_frac_r = 0;
    cases_host(0).context = false;

    Kokkos::deep_copy(cases_dev, cases_host);
    run_back_to_cell_average_kernel(cases_dev, true);
    Kokkos::deep_copy(cases_host, cases_dev);

    int failures = 0;
    for (int i = 0; i < num_cases; ++i) {
      const auto& c = cases_host(i);
#define CHECK_CONTEXT_MASK(name, scale)                                                          \
      {                                                                                          \
        const Real expected = c.context                                                          \
                                ? c.name##_input * cloud_factor(c, scale, true)                 \
                                : c.name##_input;                                                \
        const Real diff = std::abs(c.name##_output - expected);                                 \
        const Real tol = identity_tol * std::max<Real>(1, std::abs(expected));                  \
        const Real rel = diff / std::max(absolute_floor, std::abs(expected));                   \
        if (diff > tol && rel > identity_tol) {                                                  \
          std::cout << "Context-mask FAIL: " #name                                               \
                    << " case=" << i << " context=" << c.context                                 \
                    << " input=" << c.name##_input                                               \
                    << " output=" << c.name##_output                                             \
                    << " expected=" << expected << "\n";                                         \
          ++failures;                                                                            \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_CONTEXT_MASK)
#undef CHECK_CONTEXT_MASK
    }

    REQUIRE(failures == 0);
  }

  // Signed tendencies are common in microphysics (sources and sinks). This
  // check confirms scaling uses output=input*factor for either sign.
  void run_signed_input_checks() {
    const auto data_default = setup_parameter_sweep(false, true);
    run_process_location_checks(data_default);

    const auto data_separate = setup_parameter_sweep(true, true);
    run_process_location_checks(data_separate);
  }

  // Spot checks for numerical robustness across extreme magnitude regimes.
  void run_extreme_value_checks() {
    const Real identity_tol = 10 * std::numeric_limits<Real>::epsilon();
    const Real absolute_floor = kAbsoluteFloor;
    const Real denorm = std::numeric_limits<Real>::denorm_min();
    const Real tiny = denorm > 0 ? 10 * denorm : 10 * std::numeric_limits<Real>::min();

    view_1d<BackToCellAverageCase> cases_dev("back_to_cell_average_extreme_values", 3);
    auto cases_host = Kokkos::create_mirror_view(cases_dev);

    for (int i = 0; i < 3; ++i) {
      auto& c = cases_host(i);
      if (i == 0) {
        c.cld_frac_i = 1e-10;
        c.cld_frac_l = 2e-10;
        c.cld_frac_r = 3e-10;
      } else if (i == 1) {
        c.cld_frac_i = 0.8;
        c.cld_frac_l = 0.4;
        c.cld_frac_r = 0.7;
      } else {
        c.cld_frac_i = 1.0;
        c.cld_frac_l = 1.0;
        c.cld_frac_r = 1.0;
      }
      c.context = true;
      set_case_inputs(c, 1.0);

      int tidx = 0;
#define SET_EXTREME(name, scale)                                                                 \
      if (i == 0) {                                                                              \
        c.name##_input = (tidx % 2 == 0) ? 1e15 : -1e15;                                        \
      } else if (i == 1) {                                                                       \
        c.name##_input = (tidx % 2 == 0) ? tiny : -tiny;                                        \
      } else {                                                                                   \
        c.name##_input = (tidx % 2 == 0) ? 1e5 : -1e5;                                          \
      }                                                                                          \
      c.name##_output = -1;                                                                      \
      ++tidx;
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(SET_EXTREME)
#undef SET_EXTREME
    }

    Kokkos::deep_copy(cases_dev, cases_host);
    run_back_to_cell_average_kernel(cases_dev, true);
    Kokkos::deep_copy(cases_host, cases_dev);

    int failures = 0;
    for (int i = 0; i < 3; ++i) {
      const auto& c = cases_host(i);
#define CHECK_EXTREME(name, scale)                                                               \
      {                                                                                          \
        if (!std::isfinite(c.name##_output)) {                                                   \
          std::cout << "Extreme-value finite FAIL: " #name                                       \
                    << " case=" << i << " output=" << c.name##_output << "\n";                  \
          ++failures;                                                                            \
        }                                                                                        \
        const Real expected = c.name##_input * cloud_factor(c, scale, true);                    \
        const Real diff = std::abs(c.name##_output - expected);                                 \
        const Real rel = diff / std::max(absolute_floor, std::abs(expected));                   \
        if (rel > identity_tol) {                                                                \
          std::cout << "Extreme-value mapping FAIL: " #name                                      \
                    << " case=" << i                                                             \
                    << " expected=" << expected                                                  \
                    << " actual=" << c.name##_output << "\n";                                   \
          ++failures;                                                                            \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_EXTREME)
#undef CHECK_EXTREME
    }

    REQUIRE(failures == 0);
  }

  // Explicit corner regimes (zero, full, and exclusive overlaps) complement
  // the dense sweep and make expected limiting behavior explicit.
  void run_edge_case_checks() {
    const Real identity_tol = 10 * std::numeric_limits<Real>::epsilon();
    const Real absolute_floor = kAbsoluteFloor;

    view_1d<BackToCellAverageCase> cases_dev("back_to_cell_average_edge_cases", 4);
    auto cases_host = Kokkos::create_mirror_view(cases_dev);

    // Case 0: all cloud fractions are zero.
    cases_host(0).cld_frac_l = 0;
    cases_host(0).cld_frac_r = 0;
    cases_host(0).cld_frac_i = 0;
    cases_host(0).context = true;
    set_case_inputs(cases_host(0), 2e-6);

    // Case 1: full coverage.
    cases_host(1).cld_frac_l = 1;
    cases_host(1).cld_frac_r = 1;
    cases_host(1).cld_frac_i = 1;
    cases_host(1).context = true;
    set_case_inputs(cases_host(1), 3e-6);

    // Case 2: ice without full liquid overlap.
    cases_host(2).cld_frac_l = 0.3;
    cases_host(2).cld_frac_r = 0.4;
    cases_host(2).cld_frac_i = 0.6;
    cases_host(2).context = true;
    set_case_inputs(cases_host(2), 4e-6);

    // Case 3: liquid without ice.
    cases_host(3).cld_frac_l = 0.5;
    cases_host(3).cld_frac_r = 0.2;
    cases_host(3).cld_frac_i = 0.0;
    cases_host(3).context = true;
    set_case_inputs(cases_host(3), 5e-6);

    Kokkos::deep_copy(cases_dev, cases_host);
    run_back_to_cell_average_kernel(cases_dev, false);
    Kokkos::deep_copy(cases_host, cases_dev);

    int failures = 0;

    // Zero-cloud case: scaled tendencies collapse to ~0, while UNCHANGED
    // tendencies pass through.
    {
      const auto& c = cases_host(0);
#define CHECK_ZERO_CASE(name, scale)                                                             \
      {                                                                                          \
        if (is_scaled(scale)) {                                                                  \
          if (std::abs(c.name##_output) > absolute_floor) {                                     \
            std::cout << "Edge zero-cloud FAIL: " #name                                         \
                      << " output=" << c.name##_output << "\n";                                 \
            ++failures;                                                                          \
          }                                                                                      \
        } else {                                                                                 \
          const Real diff = std::abs(c.name##_output - c.name##_input);                         \
          const Real tol = identity_tol * std::max<Real>(1, std::abs(c.name##_input));          \
          if (diff > tol) {                                                                      \
            std::cout << "Edge zero-cloud unchanged FAIL: " #name                               \
                      << " input=" << c.name##_input                                             \
                      << " output=" << c.name##_output << "\n";                                  \
            ++failures;                                                                          \
          }                                                                                      \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_ZERO_CASE)
#undef CHECK_ZERO_CASE
    }

    // Full-coverage case: outputs should equal inputs, intersections should be 1.
    {
      const auto& c = cases_host(1);
      if (std::abs(c.ir_cldm - 1) > identity_tol ||
          std::abs(c.il_cldm - 1) > identity_tol ||
          std::abs(c.lr_cldm - 1) > identity_tol) {
        std::cout << "Edge full-coverage intersection FAIL\n";
        ++failures;
      }
#define CHECK_FULL_CASE(name, scale)                                                             \
      {                                                                                          \
        const Real diff = std::abs(c.name##_output - c.name##_input);                           \
        const Real tol = identity_tol * std::max<Real>(1, std::abs(c.name##_input));            \
        if (diff > tol) {                                                                        \
          std::cout << "Edge full-coverage mapping FAIL: " #name                                \
                    << " diff=" << diff << " tol=" << tol << "\n";                          \
          ++failures;                                                                            \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_FULL_CASE)
#undef CHECK_FULL_CASE
    }

    // Ice-without-liquid case: il_cldm is min(0.6, 0.3) and glaciated >= 0.3.
    {
      const auto& c = cases_host(2);
      if (std::abs(c.il_cldm - 0.3) > identity_tol) {
        std::cout << "Edge ice-without-liquid FAIL: il_cldm=" << c.il_cldm << "\n";
        ++failures;
      }
      if (c.cld_frac_glaciated < 0.3) {
        std::cout << "Edge ice-without-liquid FAIL: cld_frac_glaciated="
                  << c.cld_frac_glaciated << "\n";
        ++failures;
      }
    }

    // Liquid-without-ice case: il_cldm=0 and all IL processes should be ~0.
    {
      const auto& c = cases_host(3);
      if (std::abs(c.il_cldm) > identity_tol) {
        std::cout << "Edge liquid-without-ice FAIL: il_cldm=" << c.il_cldm << "\n";
        ++failures;
      }
#define CHECK_IL_ZERO(name, scale)                                                               \
      if (scale == ScaleKind::IL) {                                                              \
        if (std::abs(c.name##_output) > absolute_floor) {                                       \
          std::cout << "Edge liquid-without-ice IL process FAIL: " #name                        \
                    << " output=" << c.name##_output << "\n";                                \
          ++failures;                                                                            \
        }                                                                                        \
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_IL_ZERO)
#undef CHECK_IL_ZERO
    }

    REQUIRE(failures == 0);
  }

  // Informational-only total-water closure diagnostic for this tendency set.
  // Non-gating: the system includes open-process terms, so this is reported
  // as a bookkeeping indicator rather than a strict pass/fail invariant.
  void run_bookkeeping_diagnostics(const BackToCellAverageTestData& data) const {
    Real max_abs_residual = 0;
    int max_idx = -1;

    for (int i = 0; i < data.num_cases; ++i) {
      const auto& c = data.cases(i);

      const Real dqc =
          -c.qc2qr_accret_tend_output
          -c.qc2qr_autoconv_tend_output
          -c.qc2qi_hetero_freeze_tend_output
          -c.qc2qi_collect_tend_output
          -c.qc2qr_ice_shed_tend_output
          -c.qc2qi_berg_tend_output;

      const Real dqr =
          +c.qc2qr_accret_tend_output
          +c.qc2qr_autoconv_tend_output
          -c.qr2qv_evap_tend_output
          -c.qr2qi_collect_tend_output
          +c.qi2qr_melt_tend_output
          -c.qr2qi_immers_freeze_tend_output
          +c.qc2qr_ice_shed_tend_output;

      const Real dqi =
          +c.qv2qi_vapdep_tend_output
          +c.qv2qi_nucleat_tend_output
          +c.qc2qi_berg_tend_output
          +c.qr2qi_collect_tend_output
          +c.qc2qi_collect_tend_output
          +c.qr2qi_immers_freeze_tend_output
          +c.qc2qi_hetero_freeze_tend_output
          -c.qi2qv_sublim_tend_output
          -c.qi2qr_melt_tend_output;

      const Real dqv =
          -c.qv2qi_vapdep_tend_output
          +c.qi2qv_sublim_tend_output
          -c.qv2qi_nucleat_tend_output
          +c.qr2qv_evap_tend_output;

      const Real residual = dqc + dqr + dqi + dqv;
      const Real abs_residual = std::abs(residual);
      if (abs_residual > max_abs_residual) {
        max_abs_residual = abs_residual;
        max_idx = i;
      }
    }

    std::cout << "\n=== Back-to-cell-average bookkeeping diagnostics (non-gating) ===\n"
              << "  Cases checked: " << data.num_cases << "\n"
              << "  Max |total-water residual|: " << max_abs_residual << "\n"
              << "  Case index: " << max_idx << "\n"
              << "===============================================================\n";

    REQUIRE(true);
  }
};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_back_to_cell_average", "[p3_back_to_cell_average]") {
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3BackToCellAverage;

  T t;

  // Each section runs independently so failures isolate quickly.
  SECTION("process_location") {
    const auto data = t.setup_parameter_sweep();
    t.run_process_location_checks(data);
  }

  SECTION("bounds_and_scaling") {
    const auto data = t.setup_parameter_sweep();
    t.run_bounds_and_scaling_checks(data);
  }

  SECTION("cloud_fraction_logic") {
    const auto data = t.setup_parameter_sweep();
    t.run_cloud_fraction_logic_checks(data);
  }

  SECTION("runtime_option") {
    t.run_runtime_option_checks();
  }

  SECTION("context_mask") {
    t.run_context_mask_checks();
  }

  SECTION("signed_inputs") {
    t.run_signed_input_checks();
  }

  SECTION("edge_cases") {
    t.run_edge_case_checks();
  }

  SECTION("extreme_values") {
    t.run_extreme_value_checks();
  }

  SECTION("bookkeeping_diagnostics") {
    const auto data = t.setup_parameter_sweep();
    t.run_bookkeeping_diagnostics(data);
  }

  SECTION("bfb") {
    t.run_bfb();
  }
}

} // anonymous namespace
