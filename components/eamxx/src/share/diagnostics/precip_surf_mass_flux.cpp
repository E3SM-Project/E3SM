#include "precip_surf_mass_flux.hpp"

#include "share/physics/physics_constants.hpp"

namespace scream
{

PrecipSurfMassFlux::
PrecipSurfMassFlux (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto type = m_params.get<std::string>("precip_type");
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

  if (m_type & s_ice) {
    m_field_in_names.push_back("precip_ice_surf_mass");
  }
  if (m_type & s_liq) {
    m_field_in_names.push_back("precip_liq_surf_mass");
  }

  auto diag_layout = m_grid->get_2d_scalar_layout();
  FieldIdentifier fid(m_name, diag_layout, m/s, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void PrecipSurfMassFlux::compute_impl()
{
  using KT = ekat::KokkosTypes<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;

  Field::view_dev_t<const Real*> mass_liq_d, mass_ice_d;

  const auto use_liq = m_type & s_liq;
  const auto use_ice = m_type & s_ice;

  double dt = 0;
  if (use_ice) {
    auto mass_ice = m_fields_in.at("precip_ice_surf_mass");
    mass_ice_d = mass_ice.get_view<const Real*>();

    const auto& t_start = mass_ice.get_header().get_tracking().get_accum_start_time ();
    const auto& t_now   = mass_ice.get_header().get_tracking().get_time_stamp ();
    dt = t_now.seconds_from(t_start);
    if (use_liq) {
      // Ensure liq/ice have same accumulation times
      auto mass_liq = m_fields_in.at("precip_liq_surf_mass");
      mass_liq_d = mass_liq.get_view<const Real*>();
      EKAT_REQUIRE_MSG (t_now==mass_liq.get_header().get_tracking().get_time_stamp() and
                        t_start==mass_liq.get_header().get_tracking().get_accum_start_time(),
          "Error! Liquid and ice precip mass fields have different accumulation time stamps!\n");
    }
  } else if (use_liq) {
    auto mass_liq = m_fields_in.at("precip_liq_surf_mass");
    mass_liq_d = mass_liq.get_view<const Real*>();

    const auto& t_start = mass_liq.get_header().get_tracking().get_accum_start_time ();
    const auto& t_now   = mass_liq.get_header().get_tracking().get_time_stamp ();
    dt = t_now.seconds_from(t_start);
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

  Real rhodt = PC::RHO_H2O.value*dt;
  const auto& flux_view  = m_diagnostic_output.get_view<Real*>();

  const int ncols = m_grid->get_num_local_dofs();
  auto lambda = KOKKOS_LAMBDA(const Int& icol) {
    if (use_ice) {
      flux_view(icol) = mass_ice_d(icol) / rhodt;
      if (use_liq) {
        flux_view(icol) += mass_liq_d(icol) / rhodt;
      }
    } else {
      flux_view(icol) = mass_liq_d(icol) / rhodt;
    }
  };
  Kokkos::parallel_for(m_name,KT::RangePolicy(0,ncols),lambda);
}

} //namespace scream
