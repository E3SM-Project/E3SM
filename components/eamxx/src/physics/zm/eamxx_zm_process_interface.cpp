#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "eamxx_zm_process_interface.hpp"
#include "share/physics/physics_constants.hpp"

#include "zm_eamxx_bridge.hpp"

#include <ekat_assert.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_reduction_utils.hpp>

#include <mpi.h> // Include the MPI header for special print statement diagnostics

namespace scream
{

/*------------------------------------------------------------------------------------------------*/
// Constructor for the ZMDeepConvection interface
ZMDeepConvection::ZMDeepConvection( const ekat::Comm& comm,
                                    const ekat::ParameterList& params)
 : AtmosphereProcess(comm,params)
{
  // Anything that can be initialized without grid information can be initialized here.
}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::create_requests ()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  constexpr int pack_size = Pack::n;

  // Gather runtime options from file
  zm_opts.load_runtime_options(m_params);

  m_grid = m_grids_manager->get_grid("physics");

  const auto& grid_name = m_grid->name();
  const auto layout = m_grid->get_3d_scalar_layout(LEV);

  // retrieve local grid parameters
  m_ncol = m_grid->get_num_local_dofs();
  m_nlev = m_grid->get_num_vertical_levels();

  const auto m2 = pow(m,2);
  const auto s2 = pow(s,2);
  const auto K2 = pow(K,2);

  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();       // 2D variables
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);    // 3D variables at mid-levels
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(ILEV);   // 3D variables at interfaces
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(LEV,2);  // horiz_wind field

  // Input variables
  add_field<Required>("p_mid",                scalar3d_mid, Pa,    grid_name, pack_size);
  add_field<Required>("p_int",                scalar3d_int, Pa,    grid_name, pack_size);
  add_field<Required>("pseudo_density",       scalar3d_mid, Pa,    grid_name, pack_size);
  add_field<Required>("phis",                 scalar2d    , m2/s2, grid_name);
  add_field<Required>("omega",                scalar3d_mid, Pa/s,  grid_name, pack_size);
  add_field<Required>("cldfrac_tot",          scalar3d_mid, none,  grid_name, pack_size);
  add_field<Required>("pbl_height",           scalar2d    , m,     grid_name);
  add_field<Required>("landfrac",             scalar2d    , none,  grid_name);
  add_field<Required>("thl_sec",              scalar3d_int, K2,    grid_name, pack_size); // thetal variance for PBL temperature perturbation
  add_tracer<Required>("qc",                  m_grid,       kg/kg,            pack_size);

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,      grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,             pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,    grid_name, pack_size);

  // Output variables
  add_field <Updated>("precip_liq_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");
  add_field <Updated>("precip_ice_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");

  // T/qv from previous time step for DCAPE
  add_field<Computed>("zm_t_prev",            scalar3d_mid, K,     grid_name);
  add_field<Computed>("zm_q_prev",            scalar3d_mid, kg/kg, grid_name);

  // Diagnostic Outputs
  add_field<Computed>("zm_prec",              scalar2d,     m/s,  grid_name);
  add_field<Computed>("zm_snow",              scalar2d,     m/s,  grid_name);
  add_field<Computed>("zm_cape",              scalar2d,     J/kg, grid_name);
  add_field<Computed>("zm_activity",          scalar2d,     none, grid_name);

  add_field<Computed>("zm_T_mid_tend",        scalar3d_mid, K/s,    grid_name, pack_size);
  add_field<Computed>("zm_qv_tend",           scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("zm_u_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
  add_field<Computed>("zm_v_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);

} // ZMDeepConvection::create_requests

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::initialize_impl (const RunType)
{
  // Set property checks for fields in this process
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("qv"),m_grid,1e-13,0.2,true);

  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_liq_surf_mass"),m_grid,0.0,false);
  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_ice_surf_mass"),m_grid,0.0,false);

  //----------------------------------------------------------------------------
  if (zm_opts.use_fortran_bridge) {
    // allocate host mirror input variables
    zm_input.h_phis       = ZMF::view_1dh<Scalar>("zm_input.h_phis",       m_ncol);
    zm_input.h_pblh       = ZMF::view_1dh<Scalar>("zm_input.h_pblh",       m_ncol);
    zm_input.h_tpert      = ZMF::view_1dh<Scalar>("zm_input.h_tpert",      m_ncol);
    zm_input.h_landfrac   = ZMF::view_1dh<Scalar>("zm_input.h_landfrac",   m_ncol);
    zm_input.h_z_mid      = ZMF::view_2dh<Real>  ("zm_input.h_z_mid",      m_ncol, m_nlev);
    zm_input.h_p_mid      = ZMF::view_2dh<Real>  ("zm_input.h_p_mid",      m_ncol, m_nlev);
    zm_input.h_p_del      = ZMF::view_2dh<Real>  ("zm_input.h_p_del",      m_ncol, m_nlev);
    zm_input.h_T_mid      = ZMF::view_2dh<Real>  ("zm_input.h_T_mid",      m_ncol, m_nlev);
    zm_input.h_qv         = ZMF::view_2dh<Real>  ("zm_input.h_qv",         m_ncol, m_nlev);
    zm_input.h_uwind      = ZMF::view_2dh<Real>  ("zm_input.h_uwind",      m_ncol, m_nlev);
    zm_input.h_vwind      = ZMF::view_2dh<Real>  ("zm_input.h_vwind",      m_ncol, m_nlev);
    zm_input.h_omega      = ZMF::view_2dh<Real>  ("zm_input.h_omega",      m_ncol, m_nlev);
    zm_input.h_cldfrac    = ZMF::view_2dh<Real>  ("zm_input.h_cldfrac",    m_ncol, m_nlev);
    zm_input.h_z_int      = ZMF::view_2dh<Real>  ("zm_input.h_z_int",      m_ncol, m_nlev+1);
    zm_input.h_p_int      = ZMF::view_2dh<Real>  ("zm_input.h_p_int",      m_ncol, m_nlev+1);
    // allocate host mirror output variables
    zm_output.h_activity  = ZMF::view_1dh<Int>   ("zm_output.h_activity",  m_ncol);
    zm_output.h_prec      = ZMF::view_1dh<Scalar>("zm_output.h_prec",      m_ncol);
    zm_output.h_snow      = ZMF::view_1dh<Scalar>("zm_output.h_snow",      m_ncol);
    zm_output.h_cape      = ZMF::view_1dh<Scalar>("zm_output.h_cape",      m_ncol);
    zm_output.h_tend_t    = ZMF::view_2dh<Real>  ("zm_output.h_tend_t",    m_ncol, m_nlev);
    zm_output.h_tend_qv   = ZMF::view_2dh<Real>  ("zm_output.h_tend_qv",   m_ncol, m_nlev);
    zm_output.h_tend_u    = ZMF::view_2dh<Real>  ("zm_output.h_tend_u",    m_ncol, m_nlev);
    zm_output.h_tend_v    = ZMF::view_2dh<Real>  ("zm_output.h_tend_v",    m_ncol, m_nlev);
    zm_output.h_rain_prod = ZMF::view_2dh<Real>  ("zm_output.h_rain_prod", m_ncol, m_nlev);
    zm_output.h_snow_prod = ZMF::view_2dh<Real>  ("zm_output.h_snow_prod", m_ncol, m_nlev);
    zm_output.h_prec_flux = ZMF::view_2dh<Real>  ("zm_output.h_prec_flux", m_ncol, m_nlev+1);
    zm_output.h_snow_flux = ZMF::view_2dh<Real>  ("zm_output.h_snow_flux", m_ncol, m_nlev+1);
    zm_output.h_mass_flux = ZMF::view_2dh<Real>  ("zm_output.h_mass_flux", m_ncol, m_nlev+1);
    // initialize variables on the fortran side
    zm::zm_eamxx_bridge_init(m_nlev);
  } // if zm_opts.use_fortran_bridge
  //----------------------------------------------------------------------------
  if (not zm_opts.use_fortran_bridge) {
    // Determine upper limit level index of deep convection based on the reference pressure profile
    const auto hyai_h = m_grid->get_geometry_data("hyai").get_view<const Real*, Host>();
    const auto hybi_h = m_grid->get_geometry_data("hybi").get_view<const Real*, Host>();
    const auto ps0 = PC::P0.value;
    zm_opts.limcnv = -1;
    if (ps0*hyai_h(0) + ps0*hybi_h(0) >= zm_opts.upper_limit_pref) {
      zm_opts.limcnv = 0;
    } else {
      for (int k = 0; k < m_nlev; ++k) {
        Real pk0 = ps0*hyai_h(k)   + ps0*hybi_h(k);
        Real pk1 = ps0*hyai_h(k+1) + ps0*hybi_h(k+1);
        if (pk0 < zm_opts.upper_limit_pref && pk1 >= zm_opts.upper_limit_pref) {
          zm_opts.limcnv = k;
          break;
        }
      }
      if (zm_opts.limcnv == -1) { zm_opts.limcnv = m_nlev+1; }
    }
  } // if not zm_opts.use_fortran_bridge
  //----------------------------------------------------------------------------
} // ZMDeepConvection::initialize_impl

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::run_impl (const double dt)
{
  const int nlev_mid = m_nlev;
  const int nlev_int = m_nlev+1;

  // calculate_z_int() contains a team-level parallel_scan, which requires a special policy
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_ncol, nlev_mid);

  auto ts_start      = start_of_step_ts();
  bool is_first_step = (ts_start.get_num_steps()==0);

  //----------------------------------------------------------------------------
  // get fields

  // variables not updated by ZM
  zm_input.phis        = get_field_in("phis")          .get_view<const Real*>();
  zm_input.landfrac    = get_field_in("landfrac")      .get_view<const Real*>();
  zm_input.pblh        = get_field_in("pbl_height")    .get_view<const Real*>();
  zm_input.p_mid       = get_field_in("p_mid")         .get_view<const Real**>();
  zm_input.p_int       = get_field_in("p_int")         .get_view<const Real**>();
  zm_input.p_del       = get_field_in("pseudo_density").get_view<const Real**>();
  zm_input.omega       = get_field_in("omega")         .get_view<const Real**>();
  zm_input.cldfrac     = get_field_in("cldfrac_tot")   .get_view<const Real**>();
  zm_input.thl_sec     = get_field_in("thl_sec")       .get_view<const Real**>();
  zm_input.qc          = get_field_in("qc")            .get_view<const Real**>();

  // variables updated by ZM
  const auto& T_mid    = get_field_out("T_mid")      .get_view<Real**>();
  const auto& qv       = get_field_out("qv")         .get_view<Real**>();
  const auto& uwind    = get_field_out("horiz_winds").get_component(0).get_view<Real**>();
  const auto& vwind    = get_field_out("horiz_winds").get_component(1).get_view<Real**>();
  const auto& t_prev   = get_field_out("zm_t_prev")  .get_view<Real**>();
  const auto& q_prev   = get_field_out("zm_q_prev")  .get_view<Real**>();

  zm_input.T_mid = T_mid;
  zm_input.qv    = qv;
  zm_input.uwind = uwind;
  zm_input.vwind = vwind;

  zm_input.t_prev = t_prev;
  zm_input.q_prev = q_prev;

  const auto& precip_liq_surf_mass = get_field_out("precip_liq_surf_mass").get_view<Real*>();
  const auto& precip_ice_surf_mass = get_field_out("precip_ice_surf_mass").get_view<Real*>();

  //----------------------------------------------------------------------------
  // initialize output buffer variables
  zm_output.init(m_ncol, nlev_mid);

  //----------------------------------------------------------------------------
  // calculate altitude on interfaces (z_int) and mid-points (z_mid)

  // create temporaries to avoid "Implicit capture" warning
  const auto loc_zm_input_p_mid = zm_input.p_mid;
  const auto loc_zm_input_p_del = zm_input.p_del;
  const auto loc_zm_input_T_mid = zm_input.T_mid;
  const auto loc_zm_input_qv    = zm_input.qv;
  auto loc_zm_input_z_mid = zm_input.z_mid;
  auto loc_zm_input_z_del = zm_input.z_del;
  auto loc_zm_input_z_int = zm_input.z_int;
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    const auto p_mid_i = ekat::subview(loc_zm_input_p_mid, i);
    const auto p_del_i = ekat::subview(loc_zm_input_p_del, i);
    const auto T_mid_i = ekat::subview(loc_zm_input_T_mid, i);
    const auto qv_i    = ekat::subview(loc_zm_input_qv,    i);
    auto z_mid_i = ekat::subview(loc_zm_input_z_mid, i);
    auto z_del_i = ekat::subview(loc_zm_input_z_del, i);
    auto z_int_i = ekat::subview(loc_zm_input_z_int, i);
    auto z_surf = 0.0; // ZM expects z_mid & z_int to be altitude above the surface
    PF::calculate_dz(team, p_del_i, p_mid_i, T_mid_i, qv_i, z_del_i);
    team.team_barrier();
    PF::calculate_z_int(team, nlev_mid, z_del_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, nlev_mid, z_int_i, z_mid_i);
    team.team_barrier();
  });

  //----------------------------------------------------------------------------
  // calculate temperature perturbation from SHOC thetal varance for ZM parcel/CAPE

  zm_input.calculate_tpert(m_ncol, nlev_mid, is_first_step);

  //----------------------------------------------------------------------------
  // run the ZM scheme

  if (zm_opts.use_fortran_bridge) {

    zm_eamxx_bridge_run(m_ncol, nlev_mid, dt, is_first_step, zm_input, zm_output, zm_opts);

  } else {

    ZMF::zm_conv_main(zm_opts, m_ncol, nlev_mid, nlev_int, is_first_step, dt,
                      // Inputs
                      zm_input.T_mid,
                      zm_input.qv,
                      zm_input.omega,
                      zm_input.p_mid,
                      zm_input.p_int,
                      zm_input.p_del,
                      zm_input.phis,
                      zm_input.z_mid,
                      zm_input.z_int,
                      zm_input.pblh,
                      zm_input.tpert,
                      zm_input.landfrac,
                      zm_input.t_prev,
                      zm_input.q_prev,
                      // Outputs
                      zm_output.msemax_klev,
                      zm_output.jctop,
                      zm_output.jcbot,
                      zm_output.jt,
                      zm_output.prec,
                      zm_output.tend_t,
                      zm_output.tend_qv,
                      zm_output.cape,
                      zm_output.dcape,
                      zm_output.mass_flux,
                      zm_output.prec_flux,
                      zm_output.zdu,
                      zm_output.mflx_up,
                      zm_output.entr_up,
                      zm_output.detr_up,
                      zm_output.mflx_dn,
                      zm_output.entr_dn,
                      zm_output.p_del_mb,
                      zm_output.dsubcld,
                      zm_output.ql,
                      zm_output.rliq,
                      zm_output.rain_prod,
                      zm_output.dlf );

    //--------------------------------------------------------------------------
    // MCSP modifies tendencies from zm_conv_main() prior to updating the state

    // initialize local output tendencies for MCSP
    // zm_tend_init( ncol, pver, local_tend_s, local_tend_q, local_tend_u, local_tend_v )

    // perform the MCSP calculations
    // call zm_conv_mcsp_tend( ncol, ncol, pver, pverp, &
    //                         dtime, jctop, zm_const, zm_param, &
    //                         state_p_mid, state_p_int, state_p_del, &
    //                         state_s, state_qv, state_u, state_v, &
    //                         output_tend_s, output_tend_q, &
    //                         local_tend_s, local_tend_q, &
    //                         local_tend_u, local_tend_v, &
    //                         mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
    //                         mcsp_freq, mcsp_shear, zm_depth )

    // add MCSP tendencies to ZM convective tendencies
    // do i = 1,ncol
    //   do k = 1,pver
    //     output_tend_s(i,k) = output_tend_s(i,k) + local_tend_s(i,k)
    //     output_tend_q(i,k) = output_tend_q(i,k) + local_tend_q(i,k)
    //     output_tend_u(i,k) = output_tend_u(i,k) + local_tend_u(i,k)
    //     output_tend_v(i,k) = output_tend_v(i,k) + local_tend_v(i,k)
    //   end do
    // end do

    // apply tendencies from zm_conv_main() & MCSP to local copy of state variables

    // call zm_physics_update( ncol, dtime, &
    //                         state_phis, local_state_zm, local_state_zi, &
    //                         state_p_mid, state_p_int, state_p_del, &
    //                         local_state_t, local_state_qv, &
    //                         output_tend_s, output_tend_q)

    //--------------------------------------------------------------------------
    // convective momentum transport

    // initialize local output tendencies for zm_conv_evap()
    // call zm_tend_init( ncol, pver, local_tend_s, local_tend_q, tx_wind_tend(:,:,1), tx_wind_tend(:,:,2) )

    // do i = 1,ncol
    //   do k = 1,pver
    //     tx_winds(i,k,1) = state_u(i,k)
    //     tx_winds(i,k,2) = state_v(i,k)
    //   end do
    // end do

    // call zm_transport_momentum( ncol, ncol, pver, pverp, tx_winds, 2, &
    //                             mu, md, du, eu, ed, dp, &
    //                             jt, maxg, ideep, 1, lengath, &
    //                             tx_wind_tend, tx_pguall, tx_pgdall, &
    //                             tx_icwu, tx_icwd, dtime, local_tend_s )

    // add tendencies from zm_transport_momentum() to output tendencies
    // do i = 1,ncol
    //   do k = 1,pver
    //     output_tend_s(i,k) = output_tend_s(i,k) + local_tend_s(i,k)
    //     output_tend_u(i,k) = output_tend_u(i,k) + tx_wind_tend(i,k,1)
    //     output_tend_v(i,k) = output_tend_v(i,k) + tx_wind_tend(i,k,2)
    //   end do
    // end do

    //--------------------------------------------------------------------------
    // convective tracer transport

    // this is just a placeholder for where to call zm_transport_tracer(...)

    //--------------------------------------------------------------------------
    // populate deep convection activity flag

    // if (lengath.gt.0) then
    //   do i=1,lengath
    //     output_activity(ideep(i)) = 1
    //   end do
    // end if

    //--------------------------------------------------------------------------
    // convert dry static energy tendency to temperature tendency

    Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol*nlev_mid), KOKKOS_LAMBDA (const int idx) {
      const int i = idx/nlev_mid;
      const int k = idx%nlev_mid;
      zm_output.tend_t(i,k) = zm_output.tend_t(i,k)/PC::CP.value;
    });

  }

  //----------------------------------------------------------------------------
  // create temporaries of output variables to avoid "Implicit capture" warning
  const auto loc_zm_output_prec     = zm_output.prec;
  const auto loc_zm_output_snow     = zm_output.snow;
  const auto loc_zm_output_cape     = zm_output.cape;
  const auto loc_zm_output_activity = zm_output.activity;
  const auto loc_zm_output_tend_t   = zm_output.tend_t;
  const auto loc_zm_output_tend_qv  = zm_output.tend_qv;
  const auto loc_zm_output_tend_u   = zm_output.tend_u;
  const auto loc_zm_output_tend_v   = zm_output.tend_v;

  //----------------------------------------------------------------------------
  // update prognostic fields

  if (zm_opts.apply_tendencies) {
    // accumulate surface precipitation fluxes
    Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol), KOKKOS_LAMBDA (const int i) {
      auto prec_liq = loc_zm_output_prec(i) - loc_zm_output_snow(i);
      precip_liq_surf_mass(i) += Kokkos::max(0.0,prec_liq) * PC::RHO_H2O.value * dt;
      precip_ice_surf_mass(i) += loc_zm_output_snow(i) * PC::RHO_H2O.value * dt;
    });
    // update prognostic variables
    Kokkos::parallel_for("zm_update_prognostic",KT::RangePolicy(0, m_ncol*nlev_mid), KOKKOS_LAMBDA (const int idx) {
      const int i = idx/nlev_mid;
      const int k = idx%nlev_mid;
      T_mid(i,k) += loc_zm_output_tend_t (i,k) * dt;
      qv   (i,k) += loc_zm_output_tend_qv(i,k) * dt;
      uwind(i,k) += loc_zm_output_tend_u (i,k) * dt;
      vwind(i,k) += loc_zm_output_tend_v (i,k) * dt;
      // update "previous" T/qv variabiables for DCAPE
      t_prev(i,k) = T_mid(i,k);
      q_prev(i,k) = qv   (i,k);
    });
  }

  //----------------------------------------------------------------------------
  // Update output fields

  // 2D output (no vertical dimension)
  const auto& zm_prec       = get_field_out("zm_prec")        .get_view<Real*>();
  const auto& zm_snow       = get_field_out("zm_snow")        .get_view<Real*>();
  const auto& zm_cape       = get_field_out("zm_cape")        .get_view<Real*>();
  const auto& zm_activity   = get_field_out("zm_activity")    .get_view<Real*>();
  Kokkos::parallel_for("zm_diag_outputs_2D",m_ncol, KOKKOS_LAMBDA (const int i) {
    zm_prec    (i) = loc_zm_output_prec    (i);
    zm_snow    (i) = loc_zm_output_snow    (i);
    zm_cape    (i) = loc_zm_output_cape    (i);
    zm_activity(i) = loc_zm_output_activity(i);
  });

  // 3D output (vertically resolved)
  const auto& zm_T_mid_tend = ekat::scalarize(get_field_out("zm_T_mid_tend")  .get_view<Pack**>());
  const auto& zm_qv_tend    = ekat::scalarize(get_field_out("zm_qv_tend")     .get_view<Pack**>());
  const auto& zm_u_tend     = ekat::scalarize(get_field_out("zm_u_tend")      .get_view<Pack**>());
  const auto& zm_v_tend     = ekat::scalarize(get_field_out("zm_v_tend")      .get_view<Pack**>());
  Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol*nlev_mid), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_mid;
    const int k = idx%nlev_mid;
    zm_T_mid_tend(i,k) = loc_zm_output_tend_t (i,k);
    zm_qv_tend   (i,k) = loc_zm_output_tend_qv(i,k);
    zm_u_tend    (i,k) = loc_zm_output_tend_u (i,k);
    zm_v_tend    (i,k) = loc_zm_output_tend_v (i,k);
  });

} // ZMDeepConvection::run_impl

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::finalize_impl ()
{
  // placeholder for final cleanup
}

/*------------------------------------------------------------------------------------------------*/

size_t ZMDeepConvection::requested_buffer_size_in_bytes() const
{
  const int nlev_mid = m_nlev;
  const int nlev_int = m_nlev+1;
  size_t zm_buffer_size = 0;

  zm_buffer_size+= ZMF::ZmInputState::num_1d_intgr * sizeof(Int)   * m_ncol;
  zm_buffer_size+= ZMF::ZmInputState::num_1d_scalr * sizeof(Scalar)* m_ncol;
  zm_buffer_size+= ZMF::ZmInputState::num_2d_midlv * sizeof(Real)  * m_ncol * nlev_mid;
  zm_buffer_size+= ZMF::ZmInputState::num_2d_intfc * sizeof(Real)  * m_ncol * nlev_int;

  zm_buffer_size+= ZMF::ZmOutputTend::num_1d_intgr * sizeof(Int)   * m_ncol;
  zm_buffer_size+= ZMF::ZmOutputTend::num_1d_scalr * sizeof(Scalar)* m_ncol;
  zm_buffer_size+= ZMF::ZmOutputTend::num_2d_midlv * sizeof(Real)  * m_ncol * nlev_mid;
  zm_buffer_size+= ZMF::ZmOutputTend::num_2d_intfc * sizeof(Real)  * m_ncol * nlev_int;

  int num_f_mid  = (9+6);
  int num_f_int  = (2+3);
  zm_buffer_size+= num_f_mid * sizeof(Real) * m_ncol * m_nlev;
  zm_buffer_size+= num_f_int * sizeof(Real) * m_ncol * (m_nlev+1);

  return zm_buffer_size;
}

/*------------------------------------------------------------------------------------------------*/

void ZMDeepConvection::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");

  const int nlev_mid = m_nlev;
  const int nlev_int = m_nlev+1;

  constexpr auto num_1d_intgr = ZMF::ZmInputState::num_1d_intgr + ZMF::ZmOutputTend::num_1d_intgr;
  constexpr auto num_1d_scalr = ZMF::ZmInputState::num_1d_scalr + ZMF::ZmOutputTend::num_1d_scalr;
  constexpr auto num_2d_midlv = ZMF::ZmInputState::num_2d_midlv + ZMF::ZmOutputTend::num_2d_midlv;
  constexpr auto num_2d_intfc = ZMF::ZmInputState::num_2d_intfc + ZMF::ZmOutputTend::num_2d_intfc;

  constexpr int num_f_mid  = (9+6);
  constexpr int num_f_int  = (2+3);

  //----------------------------------------------------------------------------
  Int* i_mem = reinterpret_cast<Int*>(buffer_manager.get_memory());
  //----------------------------------------------------------------------------
  // device 1D integer variables
  ZMF::uview_1d<Int>* ptrs_1d_intgr[num_1d_intgr]             = { &zm_output.activity,
                                                                  &zm_output.msemax_klev,
                                                                  &zm_output.jctop,
                                                                  &zm_output.jcbot,
                                                                  &zm_output.jt,
                                                                };
  for (auto& v : ptrs_1d_intgr) {
    *v = ZMF::uview_1d<Int>(i_mem, m_ncol);
    i_mem += v->size();
  }
  //----------------------------------------------------------------------------
  Scalar* scl_mem = reinterpret_cast<Scalar*>(i_mem);
  //----------------------------------------------------------------------------
  // device 1D scalar scalars
  ZMF::uview_1d<Scalar>* ptrs_1d_scalr[num_1d_scalr]          = { &zm_input.tpert,
                                                                  &zm_output.prec,
                                                                  &zm_output.snow,
                                                                  &zm_output.cape,
                                                                  &zm_output.dcape,
                                                                  &zm_output.dsubcld,
                                                                  &zm_output.rliq,
                                                                };
  for (auto& v : ptrs_1d_scalr) {
    *v = ZMF::uview_1d<Scalar>(scl_mem, m_ncol);
    scl_mem += v->size();
  }
  //----------------------------------------------------------------------------

  // ***************************************************************************
  // TEMPORARY
  // ***************************************************************************
  Real* r_mem = reinterpret_cast<Real*>(scl_mem);
  //----------------------------------------------------------------------------
  // device 2D views on mid-point levels
  ZMF::uview_2dl<Real>* ptrs_f_midlv[num_f_mid]               = { &zm_input.f_z_mid,
                                                                  &zm_input.f_p_mid,
                                                                  &zm_input.f_p_del,
                                                                  &zm_input.f_T_mid,
                                                                  &zm_input.f_qv,
                                                                  &zm_input.f_uwind,
                                                                  &zm_input.f_vwind,
                                                                  &zm_input.f_omega,
                                                                  &zm_input.f_cldfrac,
                                                                  &zm_output.f_tend_t,
                                                                  &zm_output.f_tend_qv,
                                                                  &zm_output.f_tend_u,
                                                                  &zm_output.f_tend_v,
                                                                  &zm_output.f_rain_prod,
                                                                  &zm_output.f_snow_prod,
                                                                };
  for (auto& v : ptrs_f_midlv) {
    *v = ZMF::uview_2dl<Real>(r_mem, m_ncol, m_nlev);
    r_mem += v->size();
  }
  //----------------------------------------------------------------------------
  // device 2D views on interface levels
  ZMF::uview_2dl<Real>* ptrs_f_intfc[num_f_int]               = { &zm_input.f_z_int,
                                                                  &zm_input.f_p_int,
                                                                  &zm_output.f_prec_flux,
                                                                  &zm_output.f_snow_flux,
                                                                  &zm_output.f_mass_flux,
                                                                };
  for (auto& v : ptrs_f_intfc) {
    *v = ZMF::uview_2dl<Real>(r_mem, m_ncol, (m_nlev+1));
    r_mem += v->size();
  }
  //----------------------------------------------------------------------------
  // ***************************************************************************
  // TEMPORARY
  // ***************************************************************************
  // device 2D views on mid-point levels
  ZMF::uview_2d<Real>* ptrs_2d_midlv[num_2d_midlv]            = { &zm_input.z_mid,
                                                                  &zm_input.z_del,
                                                                  &zm_output.tend_t,
                                                                  &zm_output.tend_qv,
                                                                  &zm_output.tend_u,
                                                                  &zm_output.tend_v,
                                                                  &zm_output.rain_prod,
                                                                  &zm_output.snow_prod,
                                                                  &zm_output.zdu,
                                                                  &zm_output.mflx_up,
                                                                  &zm_output.entr_up,
                                                                  &zm_output.detr_up,
                                                                  &zm_output.mflx_dn,
                                                                  &zm_output.entr_dn,
                                                                  &zm_output.p_del_mb,
                                                                  &zm_output.ql,
                                                                  &zm_output.dlf,
                                                                };
  for (auto& v : ptrs_2d_midlv) {
    *v = ZMF::uview_2d<Real>(r_mem, m_ncol, nlev_mid);
    r_mem += v->size();
  }
  //----------------------------------------------------------------------------
  // device 2D views on interface levels
  ZMF::uview_2d<Real>* ptrs_2d_intfc[num_2d_intfc]            = { &zm_input.z_int,
                                                                  &zm_output.prec_flux,
                                                                  &zm_output.snow_flux,
                                                                  &zm_output.mass_flux,
                                                                };
  for (auto& v : ptrs_2d_intfc) {
    *v = ZMF::uview_2d<Real>(r_mem, m_ncol, nlev_int);
    r_mem += v->size();
  }
  //----------------------------------------------------------------------------
  size_t used_mem = (r_mem - buffer_manager.get_memory())*sizeof(Real);
  auto mem_chk = ( used_mem == requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory != requested memory for ZMDeepConvection.");
  //----------------------------------------------------------------------------
}

/*------------------------------------------------------------------------------------------------*/


} // namespace scream
