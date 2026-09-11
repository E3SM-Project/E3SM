#include "eamxx_spa_process_interface.hpp"

#include "share/algorithm/eamxx_data_interpolation.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_assert.hpp>
#include <ekat_units.hpp>

namespace scream
{

SPA::SPA (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  EKAT_REQUIRE_MSG(m_params.isParameter("spa_data_file"),
      "ERROR: spa_data_file is missing from SPA parameter list.");
}

// =========================================================================================
void SPA::create_requests()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  constexpr int ps = SCREAM_PACK_SIZE;

  m_model_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_model_grid->name();

  ncol_ = m_model_grid->get_num_local_dofs();
  nlev_ = m_model_grid->get_num_vertical_levels();

  // Get bands info from file, and log it
  const auto spa_data_file = m_params.get<std::string>("spa_data_file");
  nswbands_ = scorpio::get_dimlen(spa_data_file,"swband");
  nlwbands_ = scorpio::get_dimlen(spa_data_file,"lwband");
  this->log(LogLevel::info,
      "SPA data file bands dimensions:\n"
      "  - num sw bands: " + std::to_string(nswbands_) + "\n"
      "  - num lw bands: " + std::to_string(nlwbands_) + "\n");

  // Define the different field layouts that will be used for this process
  auto scalar3d_mid    = m_model_grid->get_3d_scalar_layout(LEV);
  auto scalar2d        = m_model_grid->get_2d_scalar_layout();
  auto scalar1d_mid    = m_model_grid->get_vertical_layout(LEV);
  auto scalar3d_swband = m_model_grid->get_3d_vector_layout(LEV,nswbands_,"swband");
  auto scalar3d_lwband = m_model_grid->get_3d_vector_layout(LEV,nlwbands_,"lwband");

  // Set of fields used strictly as input
  add_field<Required>("p_mid",          scalar3d_mid, Pa,    grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa,    grid_name, ps);
  add_field<Required>("T_mid",          scalar3d_mid, K,     grid_name, ps);
  add_field<Required>("qv",             scalar3d_mid, kg/kg, grid_name, ps);

  // Set of fields used strictly as output
  add_field<Computed>("nccn",        scalar3d_mid,    1/kg, grid_name, ps);
  add_field<Computed>("aero_g_sw",   scalar3d_swband, none, grid_name, ps);
  add_field<Computed>("aero_ssa_sw", scalar3d_swband, none, grid_name, ps);
  add_field<Computed>("aero_tau_sw", scalar3d_swband, none, grid_name, ps);
  add_field<Computed>("aero_tau_lw", scalar3d_lwband, none, grid_name, ps);
}

// =========================================================================================
void SPA::initialize_impl (const RunType /* run_type */)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // NOTE: SPA does not have an internal persistent state, so run_type is irrelevant

  std::vector<Field> spa_fields = {
      get_field_out("nccn").alias("CCN3"),
      get_field_out("aero_g_sw").alias("AER_G_SW"),
      get_field_out("aero_ssa_sw").alias("AER_SSA_SW"),
      get_field_out("aero_tau_sw").alias("AER_EXT_SW"),
      get_field_out("aero_tau_lw").alias("AER_EXT_LW")
  };
  auto spa_data_file = m_params.get<std::string>("spa_data_file");
  auto spa_map_file  = m_params.get<std::string>("spa_remap_file","");
  auto time_interpolation_method = m_params.get<std::string>("time_interpolation_method","yearly_periodic");

  auto pmid = get_field_in("p_mid");

  m_data_interpolation = std::make_shared<DataInterpolation>(m_model_grid,spa_fields);
  if (time_interpolation_method=="yearly_periodic") {
    m_data_interpolation->setup_periodic_time_database ({spa_data_file});
  } else if (time_interpolation_method=="linear") {
    m_data_interpolation->setup_linear_time_database ({spa_data_file});
  } else {
    EKAT_ERROR_MSG("Error! Invalid time_interpolation_method: " +
                   time_interpolation_method +
                   ". Valid options are: yearly_periodic, linear.\n");
  }

  m_data_interpolation->create_horiz_remappers(spa_map_file, m_iop_data_manager);
  DataInterpolation::VertRemapData vremap_data;
  vremap_data.vr_type = DataInterpolation::Dynamic3DRef;
  vremap_data.pname = "PS";
  vremap_data.pmid = pmid;
  m_data_interpolation->create_vert_remapper (vremap_data);
  m_data_interpolation->init_time_interpolation (start_of_step_ts(),DataInterpolation::Linear);

  dz_ = view_2d("dz_", ncol_, nlev_);

  // Set property checks for fields in this process
  using FWI = FieldWithinIntervalCheck;
  const auto eps = std::numeric_limits<double>::epsilon();

  add_postcondition_check<FWI>(get_field_out("nccn"),m_model_grid,0,1e11,true,-1e11*eps);

  // TODO: add an epslon to max possible upper bound of aero_ssa_sw?
  add_postcondition_check<FWI>(get_field_out("aero_g_sw"),m_model_grid,0.0,1.0,true);
  add_postcondition_check<FWI>(get_field_out("aero_ssa_sw"),m_model_grid,0.0,1.0,true);
  add_postcondition_check<FWI>(get_field_out("aero_tau_sw"),m_model_grid,0.0,10.0,true);
  add_postcondition_check<FWI>(get_field_out("aero_tau_lw"),m_model_grid,0.0,10.0,true);
}

// =========================================================================================
void SPA::run_impl (const double /* dt */)
{
  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using ExeSpace = KT::ExeSpace;
  using MemberType = KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<ExeSpace>;

  m_data_interpolation->run(end_of_step_ts());

  // Convert aerosol extinction (m^-1) to layer optical depth (unitless) by multiplying layer depth (m)
  const int ncol = ncol_;
  const int nlev = nlev_;
  const int nswbands = nswbands_;
  const int nlwbands = nlwbands_;

  const auto p_mid          = get_field_in ("p_mid").get_view<const Real**>();
  const auto pseudo_density = get_field_in ("pseudo_density").get_view<const Real**>();
  const auto T_mid          = get_field_in ("T_mid").get_view<const Real**>();
  const auto qv             = get_field_in ("qv").get_view<const Real**>();
  auto aero_tau_sw    = get_field_out("aero_tau_sw").get_view<Real***>();
  auto aero_tau_lw    = get_field_out("aero_tau_lw").get_view<Real***>();
  auto dz = dz_;

  const auto policy = TPF::get_default_team_policy(ncol, nlev);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const MemberType& team) {
    const int icol = team.league_rank();

    const auto pseudo_density_icol = ekat::subview(pseudo_density, icol);
    const auto p_mid_icol = ekat::subview(p_mid, icol);
    const auto T_mid_icol = ekat::subview(T_mid, icol);
    const auto qv_icol = ekat::subview(qv, icol);
    auto dz_icol = ekat::subview(dz, icol);

    PF::calculate_dz(team, pseudo_density_icol, p_mid_icol,
                      T_mid_icol, qv_icol, dz_icol);
    team.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nswbands*nlev), [&] (const int&idx) {
      auto isw = idx / nlev;
      auto k = idx % nlev;
      const auto dz_ik = dz_icol(k);
      aero_tau_sw(icol,isw,k) *= dz_ik;
    });
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlwbands*nlev), [&] (const int&idx) {
      auto ilw = idx / nlev;
      auto k = idx % nlev;
      const auto dz_ik = dz_icol(k);
      aero_tau_lw(icol,ilw,k) *= dz_ik;
    });
  }); // icol parallel_for loop
}

} // namespace scream
