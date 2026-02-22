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

// For BFB scalar compares, treat +0 and -0 as equivalent while preserving
// exact equality for all nonzero values.
KOKKOS_INLINE_FUNCTION
bool bfb_equal_allow_signed_zero(const Real a, const Real b) {
  return a == b || (a == 0 && b == 0);
}

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
// - Keeps the historical compare scope used by prior baselines.
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
      const Spack cld_frac_glaciated = max(Real(1e-4), cld_frac_i - il_cldm);

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

  // Bit-for-bit regression check for the baseline test infrastructure.
  // Kept separate from property tests that validate physics bookkeeping logic.
  void run_bfb() {
    auto engine = Base::get_engine();

    std::array<BackToCellAverageData, max_pack_size> baseline_data;
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

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
#define REQUIRE_BFB(name) REQUIRE(bfb_equal_allow_signed_zero(baseline_data[s].name, host_data[s].name));
        BACK_TO_CELL_AVERAGE_BFB_COMPARE_LIST(REQUIRE_BFB)
#undef REQUIRE_BFB
      }
    } else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        host_data(s).write(Base::m_ofile);
      }
    }
  }

  // Primary mapping check over dense sweep cases. This uses the shared
  // ScaleKind/cloud_factor table and is paired with known-answer tests below
  // that hardcode expectations independently of this helper.
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

  // Hand-calculated known-answer checks independent of cloud_factor/ScaleKind.
  // These catch mirrored mapping bugs shared by implementation and helper table.
  void run_known_answer_tests() {
    const Real identity_tol = 10 * std::numeric_limits<Real>::epsilon();

    auto run_single_case = [&](const BackToCellAverageCase& c0,
                               const bool use_separate_ice_liq_frac) {
      view_1d<BackToCellAverageCase> case_dev("back_to_cell_average_known_answer", 1);
      auto case_host = Kokkos::create_mirror_view(case_dev);
      case_host(0) = c0;
      Kokkos::deep_copy(case_dev, case_host);
      run_back_to_cell_average_kernel(case_dev, use_separate_ice_liq_frac);
      Kokkos::deep_copy(case_host, case_dev);
      return case_host(0);
    };

    int failures = 0;

    // Case A: fl=0.4, fr=0.6, fi=0.5
    // lr=min(0.4,0.6)=0.4, ir=min(0.5,0.6)=0.5, il=min(0.5,0.4)=0.4,
    // glaciated=max(1e-4,0.5-0.4)=0.1.
    BackToCellAverageCase c_a;
    c_a.cld_frac_l = 0.4;
    c_a.cld_frac_r = 0.6;
    c_a.cld_frac_i = 0.5;
    c_a.context = true;
    set_case_inputs(c_a, 1.0);
    c_a.qc2qr_accret_tend_input = 1.0;   // LR -> 0.4
    c_a.qr2qv_evap_tend_input = 1.0;     // R  -> 0.6
    c_a.qc2qr_autoconv_tend_input = 1.0; // L  -> 0.4
    c_a.qi2qr_melt_tend_input = 1.0;     // I  -> 0.5
    c_a.qc2qi_collect_tend_input = 1.0;  // IL -> 0.4
    c_a.qr2qi_collect_tend_input = 1.0;  // IR -> 0.5
    c_a.qi2qv_sublim_tend_input = 1.0;   // GLACIATED_OR_I -> 0.5 or 0.1
    c_a.qv2qi_nucleat_tend_input = 1.0;  // UNCHANGED -> 1.0

    const auto c_a_default = run_single_case(c_a, false);
    const auto c_a_separate = run_single_case(c_a, true);

#define CHECK_KNOWN(name, actual, expected)                                                      \
    {                                                                                           \
      const Real diff = std::abs((actual) - (expected));                                        \
      if (diff > identity_tol) {                                                                 \
        std::cout << "Known-answer FAIL: " #name                                                 \
                  << " actual=" << (actual)                                                      \
                  << " expected=" << (expected)                                                  \
                  << " diff=" << diff << "\n";                                                   \
        ++failures;                                                                              \
      }                                                                                          \
    }
    CHECK_KNOWN(qc2qr_accret_tend, c_a_default.qc2qr_accret_tend_output, 0.4);
    CHECK_KNOWN(qr2qv_evap_tend, c_a_default.qr2qv_evap_tend_output, 0.6);
    CHECK_KNOWN(qc2qr_autoconv_tend, c_a_default.qc2qr_autoconv_tend_output, 0.4);
    CHECK_KNOWN(qi2qr_melt_tend, c_a_default.qi2qr_melt_tend_output, 0.5);
    CHECK_KNOWN(qc2qi_collect_tend, c_a_default.qc2qi_collect_tend_output, 0.4);
    CHECK_KNOWN(qr2qi_collect_tend, c_a_default.qr2qi_collect_tend_output, 0.5);
    CHECK_KNOWN(qi2qv_sublim_tend_default, c_a_default.qi2qv_sublim_tend_output, 0.5);
    CHECK_KNOWN(qi2qv_sublim_tend_separate, c_a_separate.qi2qv_sublim_tend_output, 0.1);
    CHECK_KNOWN(qv2qi_nucleat_tend_default, c_a_default.qv2qi_nucleat_tend_output, 1.0);
    CHECK_KNOWN(qv2qi_nucleat_tend_separate, c_a_separate.qv2qi_nucleat_tend_output, 1.0);

    // Case B: fl=0.6, fr=0.2, fi=0.6
    // il=min(0.6,0.6)=0.6 -> fi-il=0.0, so glaciated floor applies: 1e-4.
    // This case targets only GLACIATED_OR_I floor activation; other tendencies
    // are initialized by set_case_inputs for completeness but not checked here.
    BackToCellAverageCase c_b;
    c_b.cld_frac_l = 0.6;
    c_b.cld_frac_r = 0.2;
    c_b.cld_frac_i = 0.6;
    c_b.context = true;
    set_case_inputs(c_b, 1.0);
    c_b.qi2qv_sublim_tend_input = 2.0;

    const auto c_b_default = run_single_case(c_b, false);
    const auto c_b_separate = run_single_case(c_b, true);
    // Expected: input * f_i = 2.0 * 0.6 = 1.2.
    CHECK_KNOWN(qi2qv_sublim_floor_default, c_b_default.qi2qv_sublim_tend_output, 2.0 * 0.6);
    // Expected: input * glaciated_floor = 2.0 * 1e-4 = 2e-4.
    CHECK_KNOWN(qi2qv_sublim_floor_separate, c_b_separate.qi2qv_sublim_tend_output, 2.0 * kGlaciatedFloor);
#undef CHECK_KNOWN

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
          /* are expected. Require sensitivity only when expected change is */                   \
          /* measurably above floating-point noise for this precision/magnitude. */              \
          const Real scaling_diff = std::abs(d0.cld_frac_i - d1.cld_frac_glaciated);            \
          const Real expected_delta = std::abs(d0.name##_input) * scaling_diff;                 \
          const Real abs_meas_floor =                                                          \
              (std::numeric_limits<Real>::digits10 <= 7) ? Real(1e-7) : Real(1e-12);            \
          const Real meas_tol = std::max(abs_meas_floor, tol);                                   \
          if (expected_delta > meas_tol && diff <= meas_tol) {                                   \
            std::cout << "Runtime branch sensitivity FAIL: " #name                              \
                      << " case=" << i                                                          \
                      << " expected_delta=" << expected_delta                                   \
                      << " actual_diff=" << diff                                                \
                      << " tol=" << meas_tol << "\n";                                          \
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

  // Extreme-value spot checks focus on numerical robustness only.
  // Mapping algebra is already covered by process_location + known_answer_tests.
  void run_extreme_value_checks() {
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
      }
      BACK_TO_CELL_AVERAGE_TENDENCY_LIST(CHECK_EXTREME)
#undef CHECK_EXTREME
      if (!std::isfinite(c.ir_cldm) || !std::isfinite(c.il_cldm) ||
          !std::isfinite(c.lr_cldm) || !std::isfinite(c.cld_frac_glaciated)) {
        std::cout << "Extreme-value finite FAIL: derived cloud fractions"
                  << " case=" << i << "\n";
        ++failures;
      }
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
      if (std::abs(c.cld_frac_glaciated - kGlaciatedFloor) > identity_tol) {
        std::cout << "Edge liquid-without-ice glaciated-floor FAIL: cld_frac_glaciated="
                  << c.cld_frac_glaciated << " expected=" << kGlaciatedFloor << "\n";
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
};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_back_to_cell_average", "[p3_back_to_cell_average]") {
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3BackToCellAverage;

  T t;

  // Each section runs independently so failures isolate quickly.
  // Tautological algebra checks (e.g., generic bounds from multiplying by
  // factors in [0,1]) are intentionally omitted in favor of mapping-focused
  // and known-answer coverage.
  SECTION("process_location") {
    const auto data = t.setup_parameter_sweep();
    t.run_process_location_checks(data);
    const auto data_mixed_sign = t.setup_parameter_sweep(false, true);
    t.run_process_location_checks(data_mixed_sign);
  }

  SECTION("known_answer_tests") {
    t.run_known_answer_tests();
  }

  SECTION("runtime_option") {
    t.run_runtime_option_checks();
  }

  SECTION("context_mask") {
    t.run_context_mask_checks();
  }

  SECTION("edge_cases") {
    t.run_edge_case_checks();
  }

  SECTION("extreme_values") {
    t.run_extreme_value_checks();
  }

  SECTION("bfb") {
    t.run_bfb();
  }
}

} // anonymous namespace
