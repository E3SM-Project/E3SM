#include "diagnostics/precip_total_surf_mass_flux.hpp"

namespace scream
{

// =========================================================================================
PrecipTotalSurfMassFluxDiagnostic::PrecipTotalSurfMassFluxDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void PrecipTotalSurfMassFluxDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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
  add_field<Required>("precip_ice_surf_mass", scalar2d_layout_mid, kg/m2, grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid(name(), scalar2d_layout_mid, m/s, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.get_header().get_alloc_properties().request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void PrecipTotalSurfMassFluxDiagnostic::compute_diagnostic_impl()
{  
  const auto& ice_mass_field = get_field_in("precip_ice_surf_mass");
  const auto& ice_mass_view  = ice_mass_field.get_view<const Real*>();
  const auto& ice_mass_track = ice_mass_field.get_header().get_tracking();
  const auto& ice_t_start    = ice_mass_track.get_accum_start_time();
  const auto& ice_t_now      = ice_mass_track.get_time_stamp();
  const auto ice_dt          = ice_t_now - ice_t_start;
  auto ice_rhodt             = PC::RHO_H2O*ice_dt;

  const auto& liq_mass_field = get_field_in("precip_liq_surf_mass");
  const auto& liq_mass_view  = liq_mass_field.get_view<const Real*>();
  const auto& liq_mass_track = liq_mass_field.get_header().get_tracking();  
  const auto& liq_t_start    = liq_mass_track.get_accum_start_time();
  const auto& liq_t_now      = liq_mass_track.get_time_stamp();
  const auto liq_dt          = liq_t_now - liq_t_start;
  auto liq_rhodt             = PC::RHO_H2O*liq_dt;

  const auto& flux_view = m_diagnostic_output.get_view<Real*>();
  
  Kokkos::parallel_for("PrecipTotalSurfMassFluxDiagnostic",
                       KT::RangePolicy(0,m_num_cols),
                       KOKKOS_LAMBDA(const Int& icol) {
    flux_view(icol) = ice_mass_view(icol)/ice_rhodt + liq_mass_view(icol)/liq_rhodt;
  });
}
// =========================================================================================
} //namespace scream
