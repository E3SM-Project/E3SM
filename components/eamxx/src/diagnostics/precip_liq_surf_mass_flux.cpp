#include "diagnostics/precip_liq_surf_mass_flux.hpp"

namespace scream
{

// =========================================================================================
PrecipLiqSurfMassFluxDiagnostic::PrecipLiqSurfMassFluxDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void PrecipLiqSurfMassFluxDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto m2 = m*m;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank

  FieldLayout scalar2d_layout_mid { {COL}, {m_num_cols} };

  // The fields required for this diagnostic to be computed
  add_field<Required>("precip_liq_surf_mass", scalar2d_layout_mid, kg/m2, grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid(name(), scalar2d_layout_mid, m/s, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.get_header().get_alloc_properties().request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void PrecipLiqSurfMassFluxDiagnostic::compute_diagnostic_impl()
{
  const auto& precip_liq_surf_mass      = get_field_in("precip_liq_surf_mass").get_view<const Real*>();
  const auto& precip_liq_surf_mass_flux = m_diagnostic_output.get_view<Real*>();
  const auto dt = m_dt;
  
  Kokkos::parallel_for("PrecipLiqMassFluxDiagnostic",
                       KT::RangePolicy(0,m_num_cols),
                       KOKKOS_LAMBDA(const Int& icol) {
    precip_liq_surf_mass_flux(icol) = precip_liq_surf_mass(icol)/PC::RHO_H2O/dt;
  });
}
// =========================================================================================
} //namespace scream
