#include "gw_functions.hpp"
#include "eamxx_gw_process_interface.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

namespace scream
{

/*------------------------------------------------------------------------------------------------*/
GWDrag::GWDrag(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params) {
  // Nothing to do here
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::create_requests()
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  constexpr int pack_size = Spack::n;

  // retrieve runtime options
  gw_opts.load_runtime_options(m_params);

  m_grid = grids_manager->get_grid("physics");

  const auto& grid_name = m_grid->name();
  const auto layout     = m_grid->get_3d_scalar_layout(true);

  // retrieve local grid parameters
  m_ncol = m_grid->get_num_local_dofs();
  m_nlev = m_grid->get_num_vertical_levels();

  const auto nondim = Units::nondimensional();
  const auto m2     = pow(m,2);
  const auto s2     = pow(s,2);
  const auto K2     = pow(K,2);

  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();        // 2D variables
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);    // 3D variables at mid-levels
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);   // 3D variables at interfaces
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);  // horiz_wind field

  // Input variables
  add_field<Required>("p_mid",                scalar3d_mid, Pa,     grid_name, pack_size);
  add_field<Required>("p_int",                scalar3d_int, Pa,     grid_name, pack_size);
  add_field<Required>("pseudo_density",       scalar3d_mid, Pa,     grid_name, pack_size);
  add_field<Required>("phis",                 scalar2d    , m2/s2,  grid_name);
  add_field<Required>("omega",                scalar3d_mid, Pa/s,   grid_name, pack_size);
  add_field<Required>("landfrac",             scalar2d    , nondim, grid_name);

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,      grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qc",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qi",                   m_grid,       kg/kg,             pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,    grid_name, pack_size);

  // // Output variables
  // add_field <Updated>("precip_liq_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");
  // add_field <Updated>("precip_ice_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");

  // Diagnostic Outputs
  add_field<Computed>("gw_activity",          scalar2d,     nondim, grid_name);
  add_field<Computed>("gw_T_mid_tend",        scalar3d_mid, K/s,    grid_name, pack_size);
  add_field<Computed>("gw_qv_tend",           scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("gw_u_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
  add_field<Computed>("gw_v_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
}

/*------------------------------------------------------------------------------------------------*/
void GWDrag::initialize_impl (const RunType) {
  // Set property checks for fields in this process
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("qv"),m_grid,1e-13,0.2,true);

  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_liq_surf_mass"),m_grid,0.0,false);
  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_ice_surf_mass"),m_grid,0.0,false);

}

/*------------------------------------------------------------------------------------------------*/
void GWDrag::run_impl (const double dt) {
  const int nlev_mid_packs   = ekat::npack<Spack>(m_nlev);
  //----------------------------------------------------------------------------
  // get fields

  // variables not updated by GWD
  const auto& phis        = get_field_in("phis")          .get_view<const Real*>();
  const auto& p_mid       = get_field_in("p_mid")         .get_view<const Spack**>();
  const auto& p_int       = get_field_in("p_int")         .get_view<const Spack**>();
  const auto& p_del       = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& omega       = get_field_in("omega")         .get_view<const Spack**>();
  const auto& landfrac    = get_field_in("landfrac")      .get_view<const Real*>();

  // variables updated by GWD
  const auto& T_mid       = get_field_out("T_mid")        .get_view<Spack**>();
  const auto& qv          = get_field_out("qv")           .get_view<Spack**>();
  const auto& qc          = get_field_out("qc")           .get_view<Spack**>();
  const auto& qi          = get_field_out("qi")           .get_view<Spack**>();
  const auto& hwinds_fld  = get_field_out("horiz_winds");
  const auto& uwind       = hwinds_fld.get_component(0)   .get_view<Spack**>();
  const auto& vwind       = hwinds_fld.get_component(1)   .get_view<Spack**>();
  //----------------------------------------------------------------------------

  // // Calculate local molecular diffusivity
  // if (do_molec_diff) {
  // }

  // Convective gravity waves (Beres scheme)
  if (gw_opts.use_gw_convect) {

    // // Determine wave sources
    // GWF::gw_beres_src();

    // // Solve for the drag profile
    // GWF::gw_drag_prof();

    // add the diffusion coefficients
    // do k = 0, pver
    //   egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
    // end do

    // ! Store constituents tendencies
    // do m=1, pcnst
    //    do k = 1, pver
    //       ptend%q(:ncol,k,m) = qtgw(:,k,m)
    //    end do
    // end do

    // ! add the momentum tendencies to the output tendency arrays
    // do k = 1, pver
    //    ptend%u(:ncol,k) = utgw(:,k)
    //    ptend%v(:ncol,k) = vtgw(:,k)
    //    ptend%s(:ncol,k) = ttgw(:,k)
    // end do

    // // Momentum & energy conservation
    // GWF::momentum_energy_conservation();

  }

  // Frontally generated gravity waves
  if (gw_opts.use_gw_frontal) {
    // GWF::gw_cm_src();
    // GWF::gw_drag_prof();
    // GWF::momentum_energy_conservation();
  }

  // Orographic stationary gravity waves
  if (gw_opts.use_gw_orogrph) {
    // GWF::gw_oro_src();
    // GWF::gw_drag_prof();
  }

  // Convert the tendencies for the dry constituents to dry air basis.
  // do m = 1, pcnst
  //    if (cnst_type(m).eq.'dry') then
  //       do k = 1, pver
  //          do i = 1, ncol
  //             ptend%q(i,k,m) = ptend%q(i,k,m)*state1%pdel(i,k)/state1%pdeldry(i,k)
  //          end do
  //       end do
  //    end if
  // end do

}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::finalize_impl ()
{
  // placeholder for final cleanup
}
/*------------------------------------------------------------------------------------------------*/
} // namespace scream
