#include "diagnostics/surf_upward_latent_heat_flux.hpp"

#include "physics/share/physics_constants.hpp"

namespace scream {

// ==============================================================================
SurfaceUpwardLatentHeatFlux::
SurfaceUpwardLatentHeatFlux(const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm, params)
 , m_name("surface_upward_latent_heat_flux")
 , cf_long_name("surface_upward_latent_heat_flux_due_to_evaporation")
{
  // In the future we may add options to include latent heat fluxes due to other water species.
  // See precip_surf_mass_flux.hpp and *.cpp for an example.
  // We'll need to change the cf_long_name, too, when this happens.
}

// ==============================================================================
void SurfaceUpwardLatentHeatFlux::
set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  Units m2(m*m,"m2");

  const auto surf_evap_units = ekat::units::kg / m2 / ekat::units::s;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank

  FieldLayout scalar2d_layout_mid { {ShortFieldTagsNames::COL}, {m_num_cols} };

  // The fields required for this diagnostic to be computed
  // surf_evap is defined by SurfaceCouplingImporter
  add_field<Required>("surf_evap", scalar2d_layout_mid, surf_evap_units,  grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid(m_name, scalar2d_layout_mid, W/m2, grid_name);
  // handle parent class member variables
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.get_header().get_alloc_properties().request_allocation();
  m_diagnostic_output.allocate_view();
}

// ==============================================================================
void SurfaceUpwardLatentHeatFlux::compute_diagnostic_impl()
{
  using KT = ekat::KokkosTypes<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;

  constexpr auto latent_heat_evap = PC::LatVap; // [J/kg]

  Field::view_dev_t<const Real*> evap_view_d;

  auto evap = get_field_in("surf_evap");
  evap_view_d = evap.get_view<const Real*>();
  const auto& flux_view = m_diagnostic_output.get_view<Real*>();
  Kokkos::parallel_for("SurfaceUpwardLatentHeatFlux",
    KT::RangePolicy(0, m_num_cols),
    KOKKOS_LAMBDA (const Int& icol) {
      flux_view(icol) = evap_view_d(icol) * latent_heat_evap;
  });
}

} // namespace scream
