#include "eamxx_cld_fraction_process_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include <ekat/ekat_assert.hpp>
#include <ekat/util/ekat_units.hpp>
#include <ekat/util/ekat_units.hpp>

#include <array>

#ifdef EAMXX_ENABLE_CLDFRAC_LAPIS_SUPPORT
#include "cldfrac_ml_lapis.hpp"
#endif

#ifdef EAMXX_HAS_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

namespace scream
{
  using namespace cld_fraction;
// =========================================================================================
CldFraction::CldFraction (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void CldFraction::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto nondim = Units::nondimensional();

  m_grid = grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // Set of fields used strictly as input
  constexpr int ps = Pack::n;
  add_tracer<Required>("qi", m_grid, kg/kg, ps);
  add_field<Required>("cldfrac_liq", scalar3d_layout_mid, nondim, grid_name,ps);

  // Set of fields used strictly as output
  add_field<Computed>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name,ps);
  add_field<Computed>("cldfrac_ice", scalar3d_layout_mid, nondim, grid_name,ps);
  // Note, we track two versions of the cloud fraction.  The versions below have "_for_analysis"
  // attached to the name because they're meant for use with fields that are exclusively
  // related to writing output.  This is an important distinction here because the internal ice
  // cloud fraction needs to be 100% whenever any ice at all is present in the cell (in order
  // for the model's ice processes to act on that cell). Folks evaluating cloud, on the other hand,
  // expect cloud fraction to represent cloud visible to the human eye (which corresponds to
  // ~1e-5 kg/kg).
  add_field<Computed>("cldfrac_tot_for_analysis", scalar3d_layout_mid, nondim, grid_name,ps);
  add_field<Computed>("cldfrac_ice_for_analysis", scalar3d_layout_mid, nondim, grid_name,ps);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.

  // Gather parameters for ice cloud thresholds from parameter list:
  m_icecloud_threshold = m_params.get<double>("ice_cloud_threshold",1e-12);  // Default = 1e-12
  m_icecloud_for_analysis_threshold = m_params.get<double>("ice_cloud_for_analysis_threshold",1e-5); // Default = 1e-5
}

// =========================================================================================
void CldFraction::initialize_impl (const RunType /* run_type */)
{
  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  add_postcondition_check<Interval>(get_field_out("cldfrac_ice"),m_grid,0.0,1.0,false);
  add_postcondition_check<Interval>(get_field_out("cldfrac_tot"),m_grid,0.0,1.0,false);
  add_postcondition_check<Interval>(get_field_out("cldfrac_ice_for_analysis"),m_grid,0.0,1.0,false);
  add_postcondition_check<Interval>(get_field_out("cldfrac_tot_for_analysis"),m_grid,0.0,1.0,false);

#ifdef EAMXX_HAS_PYTHON
  if (m_py_module.has_value()) {
    const auto& py_module = std::any_cast<const py::module&>(m_py_module);
    py_module.attr("init")();
  }
#endif
#ifdef EAMXX_ENABLE_CLDFRAC_LAPIS_SUPPORT
  if (m_params.get<bool>("use_lapis_emulator")) {
    lapis_initialize();
  }
#endif
}

// =========================================================================================
void CldFraction::run_impl (const double /* dt */)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio.
  auto qi   = get_field_in("qi");
  auto liq_cld_frac = get_field_in("cldfrac_liq");
  auto ice_cld_frac = get_field_out("cldfrac_ice");
  auto tot_cld_frac = get_field_out("cldfrac_tot");
  auto ice_cld_frac_4out = get_field_out("cldfrac_ice_for_analysis");
  auto tot_cld_frac_4out = get_field_out("cldfrac_tot_for_analysis");

#ifdef EAMXX_ENABLE_CLDFRAC_LAPIS_SUPPORT
  if (m_params.get<bool>("use_lapis_emulator")) {
    printf("Running lapis code...\n");
    // LAPIS-generated emulator is compiled in and was requested. Use it
    using KT = KokkosTypes<DefaultDevice>;
    using ExeSpace = typename KT::ExeSpace;
    using ESU = ekat::ExeSpaceUtils<ExeSpace>;
    using MemberType = typename KT::MemberType;

    auto policy = ESU::get_default_team_policy(m_num_cols, m_num_levs);

    auto qi_v                = qi.get_view<const Real**>();
    auto liq_cld_frac_v      = liq_cld_frac.get_view<const Real**>();
    auto ice_cld_frac_v      = ice_cld_frac.get_view<Real**>();
    auto tot_cld_frac_v      = tot_cld_frac.get_view<Real**>();
    auto ice_cld_frac_4out_v = ice_cld_frac_4out.get_view<Pack**>();
    auto tot_cld_frac_4out_v = tot_cld_frac_4out.get_view<Pack**>();

    constexpr int max_shared = 4096;
    // This is a number that fits within all modern
    // Determine L0 and L1 per-team scratch size
    constexpr int scratch0_required = forward_L0_scratch_required(max_shared);
    constexpr int scratch1_required = forward_L1_scratch_required(max_shared);
    GlobalViews_forward globalViews;
    policy.set_scratch_size(0, Kokkos::PerTeam(scratch0_required));
    policy.set_scratch_size(1, Kokkos::PerTeam(scratch1_required));

    // Execute the appropriate specialization based on shared amount
    auto lambda = KOKKOS_LAMBDA (const MemberType& team) {
      int icol = team.league_rank();
      auto qi_col = ekat::subview(qi_v,icol);
      auto liq_col = ekat::subview(liq_cld_frac_v,icol);
      auto ice_col = ekat::subview(ice_cld_frac_v,icol);
      auto tot_col = ekat::subview(tot_cld_frac_v,icol);

      char* scratch0 = (char*)(team.team_scratch(0).get_shmem(scratch0_required));
      char* scratch1 = (char*)(team.team_scratch(1).get_shmem(scratch1_required));

      forward<ExeSpace, max_shared>(team, globalViews, ice_col, tot_col, qi_col, liq_col, scratch0, scratch1);
    };

    Kokkos::parallel_for(policy,lambda);
  }
  return;
#endif

#ifdef EAMXX_HAS_PYTHON
  if (m_py_module.has_value()) {
    printf("Running python module '%s'...\n",m_params.get<std::string>("py_module_name").c_str());
    // Python support is compiled in and was requested. Use it
    const auto& py_fields = m_py_fields_host.at(m_grid->name());

    const auto& py_qi                = std::any_cast<const py::array&>(py_fields.at("qi"));
    const auto& py_liq_cld_frac      = std::any_cast<const py::array&>(py_fields.at("cldfrac_liq"));
    const auto& py_ice_cld_frac      = std::any_cast<const py::array&>(py_fields.at("cldfrac_ice"));
    const auto& py_tot_cld_frac      = std::any_cast<const py::array&>(py_fields.at("cldfrac_tot"));
    const auto& py_ice_cld_frac_4out = std::any_cast<const py::array&>(py_fields.at("cldfrac_ice_for_analysis"));
    const auto& py_tot_cld_frac_4out = std::any_cast<const py::array&>(py_fields.at("cldfrac_tot_for_analysis"));

    // For now, we run Python code only on CPU, so sync inputs to host
    liq_cld_frac.sync_to_host();

    const auto& py_module = std::any_cast<const py::module&>(m_py_module);
    double ice_threshold      = m_params.get<double>("ice_cloud_threshold");
    double ice_4out_threshold = m_params.get<double>("ice_cloud_for_analysis_threshold");
    py_module.attr("main")(ice_threshold,ice_4out_threshold,py_qi,py_liq_cld_frac,py_ice_cld_frac,py_tot_cld_frac,py_ice_cld_frac_4out,py_tot_cld_frac_4out);

    // Sync outputs to dev
    qi.sync_to_dev();
    liq_cld_frac.sync_to_dev();
    ice_cld_frac.sync_to_dev();
    tot_cld_frac.sync_to_dev();
    ice_cld_frac_4out.sync_to_dev();
    tot_cld_frac_4out.sync_to_dev();

    return;
  }
#endif

  printf("Running regular c++ code...\n");

  // Either we don't have python/LAPIS support, or it was not requested. Run the regular C++ impl
  auto qi_v                = qi.get_view<const Pack**>();
  auto liq_cld_frac_v      = liq_cld_frac.get_view<const Pack**>();
  auto ice_cld_frac_v      = ice_cld_frac.get_view<Pack**>();
  auto tot_cld_frac_v      = tot_cld_frac.get_view<Pack**>();
  auto ice_cld_frac_4out_v = ice_cld_frac_4out.get_view<Pack**>();
  auto tot_cld_frac_4out_v = tot_cld_frac_4out.get_view<Pack**>();

  CldFractionFunc::main(m_num_cols,m_num_levs,m_icecloud_threshold,m_icecloud_for_analysis_threshold,
    qi_v,liq_cld_frac_v,ice_cld_frac_v,tot_cld_frac_v,ice_cld_frac_4out_v,tot_cld_frac_4out_v);
}

// =========================================================================================
void CldFraction::finalize_impl()
{
#ifdef EAMXX_ENABLE_CLDFRAC_LAPIS_SUPPORT
  if (m_params.get<bool>("use_lapis_emulator")) {
    printf("finalizing lapis...\n");
    lapis_finalize();
  }
#endif
}
// =========================================================================================

} // namespace scream
