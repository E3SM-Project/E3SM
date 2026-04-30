#include "surf_upward_latent_heat_flux.hpp"

#include "share/physics/physics_constants.hpp"

namespace scream {

// ==============================================================================
SurfaceUpwardLatentHeatFlux::
SurfaceUpwardLatentHeatFlux(const ekat::Comm& comm, const ekat::ParameterList& params,
                            const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  using namespace ekat::units;

  m_field_in_names.push_back("surf_evap");

  auto m2 = (m*m).rename("m2");
  auto diag_layout = m_grid->get_2d_scalar_layout();
  auto name = "surface_upward_latent_heat_flux";
  FieldIdentifier fid(name, diag_layout, W/m2, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void SurfaceUpwardLatentHeatFlux::compute_impl()
{
  using KT = ekat::KokkosTypes<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;

  constexpr Real latent_heat_evap = PC::LatVap.value; // [J/kg]

  const int ncols = m_grid->get_num_local_dofs();
  const auto& surf_evap = m_fields_in.at("surf_evap").get_view<const Real*>();
  const auto& flux      = m_diagnostic_output.get_view<Real*>();

  typename KT::RangePolicy policy(0, ncols);
  auto lambda = KOKKOS_LAMBDA (const Int& icol) {
      flux(icol) = surf_evap(icol) * latent_heat_evap;
  };
  Kokkos::parallel_for("SurfaceUpwardLatentHeatFlux",policy,lambda);
}

} // namespace scream
