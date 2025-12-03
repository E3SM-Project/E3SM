#include "eamxx_cld_frac_net_process_interface.hpp"

#ifdef EAMXX_HAS_PYTHON
#include "share/atm_process/atmosphere_process_pyhelpers.hpp"
#endif
#include "cld_frac_net.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream
{

CldFracNet::CldFracNet (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  if (m_params.get<std::string>("emulator")=="lapis") {
    lapis_initialize();
  }
#ifdef EAMXX_HAS_PYTHON
  else if (m_params.get<std::string>("emulator")=="pytorch") {
    EKAT_REQUIRE_MSG( has_py_module(),
        "[CldFracNet] Error! Something went wrong while initializing the python module.\n");
  }
#endif
  else {
    EKAT_ERROR_MSG("[CldFracNet] Error! No valid emulator type requested.\n");
  }
}

void CldFracNet::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  const auto nondim = Units::nondimensional();
  const auto grid = grids_manager->get_grid("physics");
  const auto grid_name = grid->name();
  const auto layout = grid->get_3d_scalar_layout(true);

  // Input fields
  add_tracer<Required>("qi", grid, kg/kg);
  add_field<Required>("cldfrac_liq", layout, nondim, grid_name);

  // Output fields
  add_field<Computed>("cldfrac_tot", layout, nondim, grid_name);
  add_field<Computed>("cldfrac_ice", layout, nondim, grid_name);
}

void CldFracNet::initialize_impl (const RunType /* run_type */)
{
#ifdef EAMXX_HAS_PYTHON
  if (m_params.get<std::string>("emulator")=="pytorch") {
    py_module_call("init");
  }
#endif
}

void CldFracNet::run_impl (const double /* dt */)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio.
  auto qi  = get_field_in("qi");
  auto liq = get_field_in("cldfrac_liq");
  auto ice = get_field_out("cldfrac_ice");
  auto tot = get_field_out("cldfrac_tot");

  if (m_params.get<std::string>("emulator")=="lapis") {
    // LAPIS-generated emulator is compiled in and was requested. Use it
    using KT = KokkosTypes<DefaultDevice>;
    using ExeSpace = typename KT::ExeSpace;
    using TPF = ekat::TeamPolicyFactory<ExeSpace>;
    using MemberType = typename KT::MemberType;

    const auto dims = qi.get_header().get_identifier().get_layout().dims();
    auto policy = TPF::get_default_team_policy(dims[0], dims[1]);

    auto qi_v  = qi.get_view<const Real**>();
    auto liq_v = liq.get_view<const Real**>();
    auto ice_v = ice.get_view<Real**>();
    auto tot_v = tot.get_view<Real**>();

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
      auto liq_col = ekat::subview(liq_v,icol);
      auto ice_col = ekat::subview(ice_v,icol);
      auto tot_col = ekat::subview(tot_v,icol);

      char* scratch0 = (char*)(team.team_scratch(0).get_shmem(scratch0_required));
      char* scratch1 = (char*)(team.team_scratch(1).get_shmem(scratch1_required));

      forward<ExeSpace, max_shared>(team, globalViews, ice_col, tot_col, qi_col, liq_col, scratch0, scratch1);
    };

    Kokkos::parallel_for(policy,lambda);
  }
#ifdef EAMXX_HAS_PYTHON
  else if (m_params.get<std::string>("emulator")=="pytorch") {
    auto py_qi  = get_py_field_dev("qi");
    auto py_liq = get_py_field_dev("cldfrac_liq");
    auto py_ice = get_py_field_dev("cldfrac_ice");
    auto py_tot = get_py_field_dev("cldfrac_tot");

    double ice_threshold = m_params.get<double>("ice_cloud_threshold");

    py_module_call("forward",ice_threshold,py_qi,py_liq,py_ice,py_tot);
  }
#endif
  else {
    EKAT_ERROR_MSG("[CldFracNet] Error! No valid emulator type requested.\n");
  }
}

void CldFracNet::finalize_impl()
{
  if (m_params.get<std::string>("emulator")=="lapis") {
    lapis_finalize();
  }
}

} // namespace scream
