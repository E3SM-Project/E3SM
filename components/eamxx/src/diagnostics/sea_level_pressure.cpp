#include "diagnostics/sea_level_pressure.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream
{

// =========================================================================================
SeaLevelPressureDiagnostic::
SeaLevelPressureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void SeaLevelPressureDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  const auto m2 = pow(m,2);
  const auto s2 = pow(s,2);

  auto grid  = grids_manager->get_grid("Physics");
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar2d = grid->get_2d_scalar_layout();
  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid", scalar3d, K,     grid->name());
  add_field<Required>("p_mid", scalar3d, Pa,    grid->name());
  add_field<Required>("phis",  scalar2d, m2/s2, grid->name());

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d, Pa, grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void SeaLevelPressureDiagnostic::compute_diagnostic_impl()
{
  const auto& psl   = m_diagnostic_output.get_view<Real*>();
  const auto& T_mid = get_field_in("T_mid").get_view<const Real**>();
  const auto& p_mid = get_field_in("p_mid").get_view<const Real**>();
  const auto& phis  = get_field_in("phis").get_view<const Real*>();

  int surf_lev = m_num_levs - 1;
  using RP = typename KokkosTypes<DefaultDevice>::RangePolicy;
  using PF = scream::PhysicsFunctions<DefaultDevice>;

  auto lambda = KOKKOS_LAMBDA(const int& icol) {
    psl(icol) = PF::calculate_psl(T_mid(icol,surf_lev),p_mid(icol,surf_lev),phis(icol));
  };
  Kokkos::parallel_for("SeaLevelPressure", RP(0,m_num_cols), lambda);
}

} //namespace scream
