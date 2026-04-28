#include "sea_level_pressure.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

namespace scream
{

SeaLevelPressure::
SeaLevelPressure (const ekat::Comm& comm, const ekat::ParameterList& params,
                             const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;

  m_field_in_names.push_back("T_mid");
  m_field_in_names.push_back("p_mid");
  m_field_in_names.push_back("phis");

  auto diag_layout = m_grid->get_2d_scalar_layout();
  FieldIdentifier fid (name(), diag_layout, Pa, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void SeaLevelPressure::compute_diagnostic_impl()
{
  using RP = typename KokkosTypes<DefaultDevice>::RangePolicy;
  using PF = scream::PhysicsFunctions<DefaultDevice>;

  const auto& psl   = m_diagnostic_output.get_view<Real*>();
  const auto& T_mid = m_fields_in.at("T_mid").get_view<const Real**>();
  const auto& p_mid = m_fields_in.at("p_mid").get_view<const Real**>();
  const auto& phis  = m_fields_in.at("phis").get_view<const Real*>();

  const int ncols = m_grid->get_num_local_dofs();
  const int surf_lev = m_grid->get_num_vertical_levels() - 1;

  auto lambda = KOKKOS_LAMBDA(const int& icol) {
    psl(icol) = PF::calculate_psl(T_mid(icol,surf_lev),p_mid(icol,surf_lev),phis(icol));
  };
  Kokkos::parallel_for("SeaLevelPressure", RP(0,ncols), lambda);
}

} //namespace scream
