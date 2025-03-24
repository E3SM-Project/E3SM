#include "diagnostics/precip_surf_mass_flux.hpp"

#include "physics/share/physics_constants.hpp"

namespace scream
{

// ===============================================================================
PrecipSurfMassFlux::
PrecipSurfMassFlux (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  ci_string type = m_params.get<std::string>("precip_type");
  if (type=="ice") {
    m_type = s_ice;
  } else if (type=="liq") {
    m_type = s_liq;
  } else if (type=="total") {
    m_type = s_ice + s_liq;
  } else {
    EKAT_ERROR_MSG ("Error! Invalid choice for 'precip_type': " + type + "\n");
  }
  m_name = "precip_" + type + "_surf_mass_flux";
}

// ==============================================================================
void PrecipSurfMassFlux::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto m2 = m*m;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank

  FieldLayout scalar2d_layout_mid { {COL}, {m_num_cols} };

  // The fields required for this diagnostic to be computed
  if (m_type & s_ice) {
    add_field<Required>("precip_ice_surf_mass", scalar2d_layout_mid, kg/m2, grid_name);
  }
  if (m_type & s_liq) {
    add_field<Required>("precip_liq_surf_mass", scalar2d_layout_mid, kg/m2, grid_name);
  }

  // Construct and allocate the diagnostic field
  FieldIdentifier fid(m_name, scalar2d_layout_mid, m/s, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.get_header().get_alloc_properties().request_allocation();
  m_diagnostic_output.allocate_view();
}
// ===============================================================================
void PrecipSurfMassFlux::compute_diagnostic_impl()
{
  using KT = ekat::KokkosTypes<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;

  Field::view_dev_t<const Real*> mass_liq_d, mass_ice_d;

  const auto use_liq = m_type & s_liq;
  const auto use_ice = m_type & s_ice;

  std::int64_t dt;
  if (use_ice) {
    auto mass_ice = get_field_in("precip_ice_surf_mass");
    mass_ice_d = mass_ice.get_view<const Real*>();

    const auto& t_start = mass_ice.get_header().get_tracking().get_accum_start_time ();
    const auto& t_now   = mass_ice.get_header().get_tracking().get_time_stamp ();
    dt = t_now-t_start;
  }
  if (use_liq) {
    auto mass_liq = get_field_in("precip_liq_surf_mass");
    mass_liq_d = mass_liq.get_view<const Real*>();

    const auto& t_start = mass_liq.get_header().get_tracking().get_accum_start_time ();
    const auto& t_now   = mass_liq.get_header().get_tracking().get_time_stamp ();
    if (use_ice) {
      EKAT_REQUIRE_MSG (dt==(t_now-t_start),
          "Error! Liquid and ice precip mass fields have different accumulation time stamps!\n");
    } else {
      dt = t_now-t_start;
    }
  }

  if (dt==0) {
    // We can't compute a flux if no time has passed.
    // This may be happening b/c we're at t=0, and we have INSTANT output,
    // which always outputs fields at t=0. To avoid FPE's, we simply return,
    // setting the time stamp of the diag to invalid to signal that we did
    // not successfully compute it.
    m_diagnostic_output.get_header().get_tracking().invalidate_time_stamp();
    return;
  }

  auto rhodt = PC::RHO_H2O*dt;
  const auto& flux_view  = m_diagnostic_output.get_view<Real*>();
  Kokkos::parallel_for(m_name,
                       KT::RangePolicy(0,m_num_cols),
                       KOKKOS_LAMBDA(const Int& icol) {
    if (use_ice) {
      flux_view(icol) = mass_ice_d(icol) / rhodt;
      if (use_liq) {
        flux_view(icol) += mass_liq_d(icol) / rhodt;
      }
    } else {
      flux_view(icol) = mass_liq_d(icol) / rhodt;
    }
  });
}

} //namespace scream
