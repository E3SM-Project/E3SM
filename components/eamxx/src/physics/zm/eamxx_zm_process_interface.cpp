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

  // set runtime options
  ZMF::s_zm_opts.load_runtime_options(m_params);

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

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,      grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qc",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qi",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("nc",                   m_grid,       1/kg,              pack_size);
  add_tracer<Updated>("ni",                   m_grid,       1/kg,              pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,    grid_name, pack_size);

  // Output variables
  add_field <Updated>("precip_liq_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");
  add_field <Updated>("precip_ice_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");

  // T/qv from previous time step for DCAPE
  add_field<Computed>("zm_t_prev",            scalar3d_mid, K,      grid_name);
  add_field<Computed>("zm_q_prev",            scalar3d_mid, kg/kg,  grid_name);

  // Diagnostic Outputs
  add_field<Computed>("zm_prec",              scalar2d,     m/s,    grid_name);
  add_field<Computed>("zm_snow",              scalar2d,     m/s,    grid_name);
  add_field<Computed>("zm_cape",              scalar2d,     J/kg,   grid_name);
  add_field<Computed>("zm_dcape",             scalar2d,   J/kg/s,   grid_name);
  add_field<Computed>("zm_activity",          scalar2d,     none,   grid_name);

  add_field<Computed>("zm_detr_qc",           scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("zm_detr_qi",           scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("zm_detr_nc",           scalar3d_mid, 1/kg/s, grid_name, pack_size);
  add_field<Computed>("zm_detr_ni",           scalar3d_mid, 1/kg/s, grid_name, pack_size);

  add_field<Computed>("zm_mass_flux_int",     scalar3d_int, kg/m2/s,grid_name, pack_size);
  add_field<Computed>("zm_mflx_up",           scalar3d_mid, kg/m2/s,grid_name, pack_size);
  add_field<Computed>("zm_entr_up",           scalar3d_mid, 1/s,    grid_name, pack_size);
  add_field<Computed>("zm_detr_up",           scalar3d_mid, 1/s,    grid_name, pack_size);
  add_field<Computed>("zm_mflx_dn",           scalar3d_mid, kg/m2/s,grid_name, pack_size);
  add_field<Computed>("zm_entr_dn",           scalar3d_mid, 1/s,    grid_name, pack_size);

  add_field<Computed>("mcsp_ds_out",          scalar3d_mid, K/s,    grid_name, pack_size);
  add_field<Computed>("mcsp_dq_out",          scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("mcsp_du_out",          scalar3d_mid, m/s/s,  grid_name, pack_size);
  add_field<Computed>("mcsp_dv_out",          scalar3d_mid, m/s/s,  grid_name, pack_size);
  add_field<Computed>("mcsp_freq",            scalar2d,     none,   grid_name);
  add_field<Computed>("mcsp_shear",           scalar2d,     m/s,    grid_name);
  add_field<Computed>("zm_depth",             scalar2d,     Pa,     grid_name);

  add_field<Computed>("evap_ds_out",          scalar3d_mid, K/s,    grid_name, pack_size);
  add_field<Computed>("evap_dq_out",          scalar3d_mid, kg/kg/s,grid_name, pack_size);

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
  ZMF::s_zm_opts.set_limcnv(m_grid);
  if (this->get_comm().am_i_root()) ZMF::s_zm_opts.print();
  //----------------------------------------------------------------------------
  if (ZMF::s_zm_opts.use_fortran_bridge) {
    // allocate host mirror variables for fortran bridging
    zm_input.init_host_mirrors (m_ncol, m_nlev);
    zm_output.init_host_mirrors(m_ncol, m_nlev);
    // initialize variables on the fortran side
    zm::zm_eamxx_bridge_init( m_nlev, ZMF::s_zm_opts.limcnv+1,
                              ZMF::s_zm_opts.trig_dcape,
                              ZMF::s_zm_opts.trig_ull,
                              ZMF::s_zm_opts.clos_dyn_adj,
                              ZMF::s_zm_opts.mcsp_enabled);
  } // if use_fortran_bridge
  //----------------------------------------------------------------------------
} // ZMDeepConvection::initialize_impl

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::run_impl (const double dt)
{
  const int nlev_mid = m_nlev;
  const int nlev_int = m_nlev+1;

  auto ts_start      = start_of_step_ts();
  bool is_first_step = (ts_start.get_num_steps()==0);

  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const auto team_policy = TPF::get_default_team_policy(m_ncol, nlev_mid);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_ncol, nlev_mid);

  auto zm_opts = ZMF::s_zm_opts;
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

  // variables updated by ZM
  const auto& T_mid    = get_field_out("T_mid")      .get_view<Real**>();
  const auto& qv       = get_field_out("qv")         .get_view<Real**>();
  const auto& qc       = get_field_out("qc")         .get_view<Real**>();
  const auto& qi       = get_field_out("qi")         .get_view<Real**>();
  const auto& nc       = get_field_out("nc")         .get_view<Real**>();
  const auto& ni       = get_field_out("ni")         .get_view<Real**>();
  const auto& winds_f  = get_field_out("horiz_winds"); // (ncol, nwind, nlev)
  const auto  winds_v  = winds_f.get_view<Real***>();
  const auto& uwind    = winds_f.get_component(0).get_view<Real**>();
  const auto& vwind    = winds_f.get_component(1).get_view<Real**>();
  const auto& t_prev   = get_field_out("zm_t_prev")  .get_view<Real**>();
  const auto& q_prev   = get_field_out("zm_q_prev")  .get_view<Real**>();

  zm_input.T_mid = T_mid;
  zm_input.qv    = qv;
  zm_input.qc    = qc;
  zm_input.uwind = uwind;
  zm_input.vwind = vwind;

  zm_input.t_prev = t_prev;
  zm_input.q_prev = q_prev;

  if (is_first_step) {
    Kokkos::deep_copy(zm_input.t_prev, zm_input.T_mid);
    Kokkos::deep_copy(zm_input.q_prev, zm_input.qv);
  }

  const auto& precip_liq_surf_mass = get_field_out("precip_liq_surf_mass").get_view<Real*>();
  const auto& precip_ice_surf_mass = get_field_out("precip_ice_surf_mass").get_view<Real*>();

  zm_output.mcsp_ds_out = get_field_out("mcsp_ds_out").get_view<Real**>();
  zm_output.mcsp_dq_out = get_field_out("mcsp_dq_out").get_view<Real**>();
  zm_output.mcsp_du_out = get_field_out("mcsp_du_out").get_view<Real**>();
  zm_output.mcsp_dv_out = get_field_out("mcsp_dv_out").get_view<Real**>();
  zm_output.mcsp_freq   = get_field_out("mcsp_freq")  .get_view<Real*>();
  zm_output.mcsp_shear  = get_field_out("mcsp_shear") .get_view<Real*>();
  zm_output.zm_depth    = get_field_out("zm_depth")   .get_view<Real*>();

  zm_output.evap_ds_out = get_field_out("evap_ds_out").get_view<Real**>();
  zm_output.evap_dq_out = get_field_out("evap_dq_out").get_view<Real**>();

  //----------------------------------------------------------------------------
  // initialize output buffer variables
  zm_output.init_all(m_ncol, nlev_mid);

  //----------------------------------------------------------------------------
  // create local temporaries of input variables to avoid "Implicit capture" warning
  const auto loc_zm_input_phis      = zm_input.phis;
  const auto loc_zm_input_z_mid     = zm_input.z_mid;
  const auto loc_zm_input_z_del     = zm_input.z_del;
  const auto loc_zm_input_z_int     = zm_input.z_int;
  const auto loc_zm_input_p_mid     = zm_input.p_mid;
  const auto loc_zm_input_p_int     = zm_input.p_int;
  const auto loc_zm_input_p_del     = zm_input.p_del;
  const auto loc_zm_input_T_mid     = zm_input.T_mid;
  const auto loc_zm_input_qv        = zm_input.qv;
  const auto loc_zm_input_cldfrac   = zm_input.cldfrac;
  const auto loc_zm_input_uwind     = zm_input.uwind;
  const auto loc_zm_input_vwind     = zm_input.vwind;
  const auto loc_zm_input_tmp_s_mid = zm_input.tmp_s_mid;
  const auto loc_zm_input_tmp_T_mid = zm_input.tmp_T_mid;
  const auto loc_zm_input_tmp_qv    = zm_input.tmp_qv;
  const auto loc_zm_input_tmp_winds = zm_input.tmp_winds;

  //----------------------------------------------------------------------------
  // create local temporaries of final output variables to avoid "Implicit capture" warning
  const auto loc_zm_output_prec             = zm_output.prec;
  const auto loc_zm_output_snow             = zm_output.snow;
  const auto loc_zm_output_cape             = zm_output.cape;
  const auto loc_zm_output_dcape            = zm_output.dcape;
  const auto loc_zm_output_msemax_klev      = zm_output.msemax_klev;
  const auto loc_zm_output_jctop            = zm_output.jctop;
  const auto loc_zm_output_jcbot            = zm_output.jcbot;
  const auto loc_zm_output_jt               = zm_output.jt;
  const auto loc_zm_output_activity         = zm_output.activity;
  const auto loc_zm_output_tend_out_t       = zm_output.tend_out_t;
  const auto loc_zm_output_tend_out_s       = zm_output.tend_out_s;
  const auto loc_zm_output_tend_out_qv      = zm_output.tend_out_qv;
  const auto loc_zm_output_tend_out_u       = zm_output.tend_out_u;
  const auto loc_zm_output_tend_out_v       = zm_output.tend_out_v;
  const auto loc_zm_output_tend_tmp_s       = zm_output.tend_tmp_s;
  const auto loc_zm_output_tend_tmp_qv      = zm_output.tend_tmp_qv;
  const auto loc_zm_output_tend_tmp_winds   = zm_output.tend_tmp_winds;

  const auto loc_zm_output_rain_prod        = zm_output.rain_prod;
  const auto loc_zm_output_tend_s_snwprd    = zm_output.tend_s_snwprd;
  const auto loc_zm_output_tend_s_snwevmlt  = zm_output.tend_s_snwevmlt;
  const auto loc_zm_output_ntprprd          = zm_output.ntprprd;
  const auto loc_zm_output_ntsnprd          = zm_output.ntsnprd;
  const auto loc_zm_output_flxprec          = zm_output.flxprec;
  const auto loc_zm_output_flxsnow          = zm_output.flxsnow;
  const auto loc_zm_output_mcsp_ds_out      = zm_output.mcsp_ds_out;
  const auto loc_zm_output_mcsp_dq_out      = zm_output.mcsp_dq_out;
  const auto loc_zm_output_mcsp_du_out      = zm_output.mcsp_du_out;
  const auto loc_zm_output_mcsp_dv_out      = zm_output.mcsp_dv_out;
  const auto loc_zm_output_mcsp_freq        = zm_output.mcsp_freq;
  const auto loc_zm_output_mcsp_shear       = zm_output.mcsp_shear;
  const auto loc_zm_output_zm_depth         = zm_output.zm_depth;

  const auto loc_zm_output_evap_ds_out      = zm_output.evap_ds_out;
  const auto loc_zm_output_evap_dq_out      = zm_output.evap_dq_out;

  const auto loc_zm_output_mflx_up          = zm_output.mflx_up;
  const auto loc_zm_output_entr_up          = zm_output.entr_up;
  const auto loc_zm_output_detr_up          = zm_output.detr_up;
  const auto loc_zm_output_mflx_dn          = zm_output.mflx_dn;
  const auto loc_zm_output_entr_dn          = zm_output.entr_dn;
  const auto loc_zm_output_p_del_mb         = zm_output.p_del_mb;
  const auto loc_zm_output_dsubcld          = zm_output.dsubcld;
  const auto loc_zm_output_ql               = zm_output.ql;
  const auto loc_zm_output_rliq             = zm_output.rliq;
  const auto loc_zm_output_dlf              = zm_output.dlf;

  //----------------------------------------------------------------------------
  // calculate altitude on interfaces (z_int) and mid-points (z_mid)
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    const auto p_mid_i = ekat::subview(loc_zm_input_p_mid, i);
    const auto p_del_i = ekat::subview(loc_zm_input_p_del, i);
    const auto T_mid_i = ekat::subview(loc_zm_input_T_mid, i);
    const auto qv_i    = ekat::subview(loc_zm_input_qv,    i);
    const auto z_mid_i = ekat::subview(loc_zm_input_z_mid, i);
    const auto z_del_i = ekat::subview(loc_zm_input_z_del, i);
    const auto z_int_i = ekat::subview(loc_zm_input_z_int, i);
    const auto z_surf = 0.0; // ZM expects z_mid & z_int to be altitude above the surface
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
    // Allocate the Workspace for the MCSP / evap / momentum kernels below.
    // Slot length is nlev_int*nwind to fit zm_transport_momentum's (nwind, nlev_int) buffers.
    WSM wsm( nlev_int*nwind, 16, team_policy);

    //--------------------------------------------------------------------------
    // run the main ZM scheme
    ZMF::zm_conv_main(zm_opts, m_ncol, nlev_mid, nlev_int, is_first_step, dt,
                      zm_input.T_mid, zm_input.qv, zm_input.omega,
                      zm_input.p_mid, zm_input.p_int, zm_input.p_del,
                      zm_input.phis, zm_input.z_mid, zm_input.z_int,
                      zm_input.pblh, zm_input.tpert, zm_input.landfrac,
                      zm_input.t_prev, zm_input.q_prev,
                      zm_output.msemax_klev, zm_output.jctop, zm_output.jcbot, zm_output.jt,
                      zm_output.activity,
                      zm_output.prec, zm_output.tend_out_s, zm_output.tend_out_qv,
                      zm_output.cape, zm_output.dcape,
                      zm_output.mass_flux, zm_output.prec_flux,
                      zm_output.zdu,
                      zm_output.mflx_up, zm_output.entr_up, zm_output.detr_up,
                      zm_output.mflx_dn, zm_output.entr_dn,
                      zm_output.p_del_mb, zm_output.dsubcld,
                      zm_output.ql, zm_output.rliq, zm_output.rain_prod, zm_output.dlf,
                      zm_output.ktm, zm_output.kbm );

    //--------------------------------------------------------------------------
    // MCSP modifies tendencies from zm_conv_main() prior to updating the state

    if (zm_opts.mcsp_enabled) {

      // initialize intermediate output tendencies for MCSP
      zm_output.init_tmp(m_ncol, nlev_mid);

      // perform the MCSP calculations
      Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
        const Int i = team.league_rank();
        const auto phis_i         = loc_zm_input_phis(i);
        const auto jctop_i        = loc_zm_output_jctop(i);
        const auto z_mid_i        = ekat::subview(loc_zm_input_z_mid, i);
        const auto p_mid_i        = ekat::subview(loc_zm_input_p_mid, i);
        const auto p_int_i        = ekat::subview(loc_zm_input_p_int, i);
        const auto p_del_i        = ekat::subview(loc_zm_input_p_del, i);
        const auto T_mid_i        = ekat::subview(loc_zm_input_T_mid, i);
        const auto qv_i           = ekat::subview(loc_zm_input_qv, i);
        const auto uwind_i        = ekat::subview(loc_zm_input_uwind, i);
        const auto vwind_i        = ekat::subview(loc_zm_input_vwind, i);
        const auto tmp_s_mid_i    = ekat::subview(loc_zm_input_tmp_s_mid, i);
        const auto tend_out_s_i   = ekat::subview(loc_zm_output_tend_out_s, i);
        const auto tend_out_qv_i  = ekat::subview(loc_zm_output_tend_out_qv, i);
        const auto tend_out_u_i   = ekat::subview(loc_zm_output_tend_out_u, i);
        const auto tend_out_v_i   = ekat::subview(loc_zm_output_tend_out_v, i);
        const auto tend_tmp_s_i   = ekat::subview(loc_zm_output_tend_tmp_s, i);
        const auto tend_tmp_qv_i  = ekat::subview(loc_zm_output_tend_tmp_qv, i);
        const auto tend_tmp_winds_i = ekat::subview(loc_zm_output_tend_tmp_winds, i); // (nwind, nlev)
        const auto tend_tmp_u_i   = ekat::subview(tend_tmp_winds_i, 0);               // (nlev)
        const auto tend_tmp_v_i   = ekat::subview(tend_tmp_winds_i, 1);               // (nlev)
        const auto mcsp_ds_out_i  = ekat::subview(loc_zm_output_mcsp_ds_out, i);
        const auto mcsp_dq_out_i  = ekat::subview(loc_zm_output_mcsp_dq_out, i);
        const auto mcsp_du_out_i  = ekat::subview(loc_zm_output_mcsp_du_out, i);
        const auto mcsp_dv_out_i  = ekat::subview(loc_zm_output_mcsp_dv_out, i);
        auto& mcsp_freq_i         = loc_zm_output_mcsp_freq(i);
        auto& mcsp_shear_i        = loc_zm_output_mcsp_shear(i);
        auto& zm_depth_i          = loc_zm_output_zm_depth(i);

        // calculate DSE for MCSP
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_mid), [&] (const int k) {
          tmp_s_mid_i(k) = T_mid_i(k)*PC::CP.value + z_mid_i(k)*PC::gravit.value + phis_i;
        });
        team.team_barrier();

        // call the MCSP scheme
        ZMF::zm_conv_mcsp_tend( team, wsm.get_workspace(team), zm_opts, nlev_mid, nlev_int, dt,
          jctop_i, p_mid_i, p_int_i, p_del_i, tmp_s_mid_i, qv_i, uwind_i, vwind_i,
          tend_out_s_i, tend_out_qv_i,
          tend_tmp_s_i, tend_tmp_qv_i, tend_tmp_u_i, tend_tmp_v_i,
          mcsp_ds_out_i, mcsp_dq_out_i, mcsp_du_out_i, mcsp_dv_out_i,
          mcsp_freq_i, mcsp_shear_i, zm_depth_i );

        // add MCSP tendencies to output tendencies
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_mid), [&] (const int k) {
          tend_out_s_i(k)  += tend_tmp_s_i(k);
          tend_out_qv_i(k) += tend_tmp_qv_i(k);
          tend_out_u_i(k)  += tend_tmp_u_i(k);
          tend_out_v_i(k)  += tend_tmp_v_i(k);
        });
      });

    }

    //--------------------------------------------------------------------------
    // apply tendencies from zm_conv_main() (& MCSP) to temporary state variables for zm_conv_evap()

    Kokkos::parallel_for(KT::RangePolicy(0, m_ncol*nlev_mid), KOKKOS_LAMBDA (const int idx) {
      const int i = idx/nlev_mid;
      const int k = idx%nlev_mid;
      // calculate temporary T/qv for zm_conv_evap()
      loc_zm_input_tmp_T_mid(i,k) = loc_zm_input_T_mid(i,k) + loc_zm_output_tend_out_s(i,k)/PC::CP.value * dt;
      loc_zm_input_tmp_qv(i,k)    = loc_zm_input_qv(i,k)    + loc_zm_output_tend_out_qv(i,k) * dt;
      // calculate temporary winds for zm_transport_momentum()
      loc_zm_input_tmp_winds(i,0,k) = loc_zm_input_uwind(i,k) + loc_zm_output_tend_out_u(i,k) * dt;
      loc_zm_input_tmp_winds(i,1,k) = loc_zm_input_vwind(i,k) + loc_zm_output_tend_out_v(i,k) * dt;
    });

    //--------------------------------------------------------------------------
    // Compute the precipitation, rain evaporation, and snow formation/melting
    // Note - this routine expects an updated state following zm_conv_main() (+MCSP)

    // initialize intermediate output tendencies for zm_conv_evap()
    zm_output.init_tmp(m_ncol, nlev_mid);

    // perform the convective evaporation calculations
    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
      const Int i = team.league_rank();
      // skip inactive columns: zm_conv_main leaves p_del_mb (=dp) zero for them,
      // and zm_transport_momentum divides by dp, which would produce NaNs
      if (!loc_zm_output_activity(i)) return;
      const auto p_mid_i            = ekat::subview(loc_zm_input_p_mid, i);
      const auto p_del_i            = ekat::subview(loc_zm_input_p_del, i);
      const auto T_mid_i            = ekat::subview(loc_zm_input_T_mid, i);
      const auto qv_i               = ekat::subview(loc_zm_input_qv, i);
      const auto tmp_T_mid_i        = ekat::subview(loc_zm_input_tmp_T_mid, i);
      const auto tmp_qv_i           = ekat::subview(loc_zm_input_tmp_qv, i);
      const auto rain_prod_i        = ekat::subview(loc_zm_output_rain_prod, i);
      const auto cldfrac_i          = ekat::subview(loc_zm_input_cldfrac, i);
      const auto tend_out_s_i       = ekat::subview(loc_zm_output_tend_out_s, i);
      const auto tend_out_qv_i      = ekat::subview(loc_zm_output_tend_out_qv, i);
      const auto tend_tmp_s_i       = ekat::subview(loc_zm_output_tend_tmp_s, i);
      const auto tend_tmp_qv_i      = ekat::subview(loc_zm_output_tend_tmp_qv, i);
      const auto tend_s_snwprd_i    = ekat::subview(loc_zm_output_tend_s_snwprd, i);
      const auto tend_s_snwevmlt_i  = ekat::subview(loc_zm_output_tend_s_snwevmlt, i);
      const auto ntprprd_i          = ekat::subview(loc_zm_output_ntprprd, i);
      const auto ntsnprd_i          = ekat::subview(loc_zm_output_ntsnprd, i);
      const auto flxprec_i          = ekat::subview(loc_zm_output_flxprec, i);
      const auto flxsnow_i          = ekat::subview(loc_zm_output_flxsnow, i);
      const auto evap_ds_out_i      = ekat::subview(loc_zm_output_evap_ds_out, i);
      const auto evap_dq_out_i      = ekat::subview(loc_zm_output_evap_dq_out, i);

      // call the ZM evap scheme
      ZMF::zm_conv_evap( team, zm_opts, nlev_mid, nlev_int, dt,
        p_mid_i, p_del_i, tmp_T_mid_i, tmp_qv_i, rain_prod_i, cldfrac_i,
        tend_tmp_s_i, tend_tmp_qv_i, tend_s_snwprd_i, tend_s_snwevmlt_i,
        loc_zm_output_prec(i), loc_zm_output_snow(i),
        ntprprd_i, ntsnprd_i, flxprec_i, flxsnow_i );

      // add zm_conv_evap() tendencies to output tendencies and update temporary state variables
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_mid), [&] (const int k) {
        tend_out_s_i(k)  += tend_tmp_s_i(k);
        tend_out_qv_i(k) += tend_tmp_qv_i(k);
        evap_ds_out_i(k) = tend_tmp_s_i(k);
        evap_dq_out_i(k) = tend_tmp_qv_i(k);
      });

    });

    //--------------------------------------------------------------------------
    // convective momentum transport

    // initialize intermediate output tendencies for zm_transport_momentum()
    zm_output.init_tmp(m_ncol, nlev_mid);

    // domain-wide convection level bounds from zm_conv_main (host scalars captured by value)
    const Int ktm = zm_output.ktm;
    const Int kbm = zm_output.kbm;

    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
      const Int i = team.league_rank();
      // skip inactive columns: zm_conv_main leaves p_del_mb (=dp) zero for them,
      // and zm_transport_momentum divides by dp, which would produce NaNs
      if (!loc_zm_output_activity(i)) return;
      // MCSP-updated winds in (nwind, nlev) == (m,k) layout for wind_mid
      const auto wind_mid_i     = ekat::subview(loc_zm_input_tmp_winds, i);
      const auto mflx_up_i      = ekat::subview(loc_zm_output_mflx_up, i);
      const auto mflx_dn_i      = ekat::subview(loc_zm_output_mflx_dn, i);
      const auto detr_up_i      = ekat::subview(loc_zm_output_detr_up, i);
      const auto entr_up_i      = ekat::subview(loc_zm_output_entr_up, i);
      const auto entr_dn_i      = ekat::subview(loc_zm_output_entr_dn, i);
      const auto p_del_mb_i     = ekat::subview(loc_zm_output_p_del_mb, i);

      const auto tend_out_s_i   = ekat::subview(loc_zm_output_tend_out_s, i);
      const auto tend_out_u_i   = ekat::subview(loc_zm_output_tend_out_u, i);
      const auto tend_out_v_i   = ekat::subview(loc_zm_output_tend_out_v, i);

      // seten output buffer, and combined (nwind, nlev) momentum tendency output
      const auto tend_tmp_s_i   = ekat::subview(loc_zm_output_tend_tmp_s, i);
      const auto wind_tend_i    = ekat::subview(loc_zm_output_tend_tmp_winds, i);

      // call the ZM momentum transport scheme
      ZMF::zm_transport_momentum( team, wsm.get_workspace(team), nlev_mid, nlev_int, dt,
        wind_mid_i, 2,
        mflx_up_i, mflx_dn_i, detr_up_i, entr_up_i, entr_dn_i, p_del_mb_i,
        loc_zm_output_jt(i), loc_zm_output_msemax_klev(i), ktm, kbm,
        wind_tend_i,    // momentum tendency (nwind, nlev)
        tend_tmp_s_i ); // seten (KE-dissipation dry static energy tendency)

      // add zm_transport_momentum tendencies to output tendencies
      // (mirrors EAM: output_tend_{u,v} += tx_wind_tend; output_tend_s += seten)
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_mid), [&] (const int k) {
        tend_out_s_i(k)  += tend_tmp_s_i(k);
        tend_out_u_i(k)  += wind_tend_i(0,k);
        tend_out_v_i(k)  += wind_tend_i(1,k);
      });

    });

    //--------------------------------------------------------------------------
  }

  //----------------------------------------------------------------------------
  // update prognostic fields

  // accumulate surface precipitation fluxes
  Kokkos::parallel_for("zm_update_precip",
    KT::RangePolicy(0, m_ncol), KOKKOS_LAMBDA (const int i) {
    auto prec_liq = loc_zm_output_prec(i) - loc_zm_output_snow(i);
    precip_liq_surf_mass(i) += Kokkos::max(0.0,prec_liq) * PC::RHO_H2O.value * dt;
    precip_ice_surf_mass(i) += loc_zm_output_snow(i) * PC::RHO_H2O.value * dt;
  });

  // update 3D prognostic variables
  const auto& zm_detr_qc = get_field_out("zm_detr_qc").get_view<Real**>();
  const auto& zm_detr_qi = get_field_out("zm_detr_qi").get_view<Real**>();
  const auto& zm_detr_nc = get_field_out("zm_detr_nc").get_view<Real**>();
  const auto& zm_detr_ni = get_field_out("zm_detr_ni").get_view<Real**>();
  Kokkos::parallel_for("zm_update_prognostic",
    KT::RangePolicy(0, m_ncol*nlev_mid), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_mid;
    const int k = idx%nlev_mid;
    // partition detrained cloud water (dlf) into liquid and ice by temperature
    // (this method is equivalent to how ZM dlf was applied in CLUBB for EAMv2)
    Real ice_frac;
    if      (T_mid(i,k) > ZMF::ZMC::dlf_tk1) ice_frac = 0;
    else if (T_mid(i,k) < ZMF::ZMC::dlf_tk2) ice_frac = 1;
    else               ice_frac = (ZMF::ZMC::dlf_tk1 - T_mid(i,k))
                                 /(ZMF::ZMC::dlf_tk1 - ZMF::ZMC::dlf_tk2);
    // liq/ice mass detrainment tendencies
    zm_detr_qc(i,k) = (1 - ice_frac) * Kokkos::max(loc_zm_output_dlf(i,k), 0.0);
    zm_detr_qi(i,k) = ice_frac       * Kokkos::max(loc_zm_output_dlf(i,k), 0.0);
    // liq/ice number detrainment tendencies
    zm_detr_nc(i,k) = zm_detr_qc(i,k) / ( 4.0/3.0*PC::Pi
                                          * Kokkos::pow(ZMF::ZMC::cld_liq_radius,3.0)
                                          * ZMF::ZMC::cld_liq_density);
    zm_detr_ni(i,k) = zm_detr_qi(i,k) / ( 4.0/3.0*PC::Pi
                                          * Kokkos::pow(ZMF::ZMC::cld_ice_radius,3.0)
                                          * ZMF::ZMC::cld_ice_density);
    // convert dry static energy tendency to temperature tendency (if not using fortran bridge)
    if (not zm_opts.use_fortran_bridge) {
      loc_zm_output_tend_out_t(i,k) = loc_zm_output_tend_out_s(i,k)/PC::CP.value;
    }
    // apply tendencies (winds via the combined (ncol,nwind,nlev) field view)
    T_mid(i,k)     += loc_zm_output_tend_out_t (i,k) * dt;
    qv   (i,k)     += loc_zm_output_tend_out_qv(i,k) * dt;
    if (zm_opts.apply_detr_tend) {
      qc (i,k)     += zm_detr_qc(i,k) * dt;
      qi (i,k)     += zm_detr_qi(i,k) * dt;
      nc (i,k)     += zm_detr_nc(i,k) * dt;
      ni (i,k)     += zm_detr_ni(i,k) * dt;
    }
    winds_v(i,0,k) += loc_zm_output_tend_out_u (i,k) * dt;
    winds_v(i,1,k) += loc_zm_output_tend_out_v (i,k) * dt;
    // update "previous" T/qv variabiables for DCAPE
    t_prev(i,k) = T_mid(i,k);
    q_prev(i,k) = qv   (i,k);
  });

  //----------------------------------------------------------------------------
  // Update output fields

  // 2D diagnostic output (no vertical dimension)
  const auto& zm_prec       = get_field_out("zm_prec")    .get_view<Real*>();
  const auto& zm_snow       = get_field_out("zm_snow")    .get_view<Real*>();
  const auto& zm_cape       = get_field_out("zm_cape")    .get_view<Real*>();
  const auto& zm_dcape      = get_field_out("zm_dcape")   .get_view<Real*>();
  const auto& zm_activity   = get_field_out("zm_activity").get_view<Real*>();
  Kokkos::parallel_for("zm_diag_outputs_2D",m_ncol, KOKKOS_LAMBDA (const int i) {
    zm_prec    (i) = loc_zm_output_prec    (i);
    zm_snow    (i) = loc_zm_output_snow    (i);
    zm_cape    (i) = loc_zm_output_cape    (i);
    zm_dcape   (i) = loc_zm_output_dcape   (i);
    zm_activity(i) = loc_zm_output_activity(i);
  });

  // 3D mid-level output
  const auto& mflx_up = get_field_out("zm_mflx_up").get_view<Real**>();
  const auto& entr_up = get_field_out("zm_entr_up").get_view<Real**>();
  const auto& detr_up = get_field_out("zm_detr_up").get_view<Real**>();
  const auto& mflx_dn = get_field_out("zm_mflx_dn").get_view<Real**>();
  const auto& entr_dn = get_field_out("zm_entr_dn").get_view<Real**>();
  Kokkos::parallel_for("zm_diag_outputs_3D_mid",
    KT::RangePolicy(0, m_ncol*nlev_mid), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_mid;
    const int k = idx%nlev_mid;
    mflx_up(i,k) = loc_zm_output_mflx_up(i,k);
    entr_up(i,k) = loc_zm_output_entr_up(i,k);
    detr_up(i,k) = loc_zm_output_detr_up(i,k);
    mflx_dn(i,k) = loc_zm_output_mflx_dn(i,k);
    entr_dn(i,k) = loc_zm_output_entr_dn(i,k);
  });

  // 3D interface output
  const auto loc_zm_output_mass_flux = zm_output.mass_flux;
  const auto& mass_flux_int = get_field_out("zm_mass_flux_int").get_view<Real**>();
  Kokkos::parallel_for("zm_diag_outputs_3D_int",
    KT::RangePolicy(0, m_ncol*nlev_int), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_int;
    const int k = idx%nlev_int;
    mass_flux_int(i,k) = loc_zm_output_mass_flux(i,k);
  });

} // ZMDeepConvection::run_impl

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::finalize_impl ()
{
  ZMF::zm_finalize();
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
  zm_buffer_size+= ZMF::ZmInputState::num_3d_midlv * sizeof(Real)  * m_ncol * nwind * nlev_mid;

  zm_buffer_size+= ZMF::ZmOutputTend::num_1d_intgr * sizeof(Int)   * m_ncol;
  zm_buffer_size+= ZMF::ZmOutputTend::num_1d_scalr * sizeof(Scalar)* m_ncol;
  zm_buffer_size+= ZMF::ZmOutputTend::num_2d_midlv * sizeof(Real)  * m_ncol * nlev_mid;
  zm_buffer_size+= ZMF::ZmOutputTend::num_2d_intfc * sizeof(Real)  * m_ncol * nlev_int;
  zm_buffer_size+= ZMF::ZmOutputTend::num_3d_midlv * sizeof(Real)  * m_ncol * nwind * nlev_mid;

  // fortran-bridge (LayoutLeft) transpose buffers are only used when running
  // the fortran bridge, so only reserve space for them in that case
  if (ZMF::s_zm_opts.use_fortran_bridge) {
    constexpr int num_f_mid = ZMF::ZmInputState::num_f_midlv + ZMF::ZmOutputTend::num_f_midlv;
    constexpr int num_f_int = ZMF::ZmInputState::num_f_intfc + ZMF::ZmOutputTend::num_f_intfc;
    zm_buffer_size+= num_f_mid * sizeof(Real) * m_ncol * m_nlev;
    zm_buffer_size+= num_f_int * sizeof(Real) * m_ncol * (m_nlev+1);
  }

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

  constexpr int num_f_mid  = ZMF::ZmInputState::num_f_midlv + ZMF::ZmOutputTend::num_f_midlv;
  constexpr int num_f_int  = ZMF::ZmInputState::num_f_intfc + ZMF::ZmOutputTend::num_f_intfc;

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
  // fortran-bridge (LayoutLeft) transpose buffers are only carved out when
  // running the fortran bridge (must match requested_buffer_size_in_bytes())
  if (ZMF::s_zm_opts.use_fortran_bridge) {
    // device 2D views on mid-point levels
    ZMF::uview_2dl<Real>* ptrs_f_midlv[num_f_mid]             = { &zm_input.f_z_mid,
                                                                  &zm_input.f_p_mid,
                                                                  &zm_input.f_p_del,
                                                                  &zm_input.f_T_mid,
                                                                  &zm_input.f_qv,
                                                                  &zm_input.f_uwind,
                                                                  &zm_input.f_vwind,
                                                                  &zm_input.f_omega,
                                                                  &zm_input.f_cldfrac,
                                                                  &zm_input.f_t_prev,
                                                                  &zm_input.f_q_prev,
                                                                  &zm_output.f_tend_t,
                                                                  &zm_output.f_tend_qv,
                                                                  &zm_output.f_tend_u,
                                                                  &zm_output.f_tend_v,
                                                                  &zm_output.f_rain_prod,
                                                                  &zm_output.f_snow_prod,
                                                                  &zm_output.f_dlf,
                                                                  &zm_output.f_mcsp_ds_out,
                                                                  &zm_output.f_mcsp_dq_out,
                                                                  &zm_output.f_mcsp_du_out,
                                                                  &zm_output.f_mcsp_dv_out,
                                                                  &zm_output.f_evap_ds_out,
                                                                  &zm_output.f_evap_dq_out,
                                                                };
    for (auto& v : ptrs_f_midlv) {
      *v = ZMF::uview_2dl<Real>(r_mem, m_ncol, m_nlev);
      r_mem += v->size();
    }
    //--------------------------------------------------------------------------
    // device 2D views on interface levels
    ZMF::uview_2dl<Real>* ptrs_f_intfc[num_f_int]             = { &zm_input.f_z_int,
                                                                  &zm_input.f_p_int,
                                                                  &zm_output.f_prec_flux,
                                                                  &zm_output.f_snow_flux,
                                                                  &zm_output.f_mass_flux,
                                                                };
    for (auto& v : ptrs_f_intfc) {
      *v = ZMF::uview_2dl<Real>(r_mem, m_ncol, (m_nlev+1));
      r_mem += v->size();
    }
  }
  //----------------------------------------------------------------------------
  // ***************************************************************************
  // TEMPORARY
  // ***************************************************************************
  // device 2D views on mid-point levels
  ZMF::uview_2d<Real>* ptrs_2d_midlv[num_2d_midlv]            = { &zm_input.z_mid,
                                                                  &zm_input.z_del,
                                                                  &zm_input.tmp_s_mid,
                                                                  &zm_input.tmp_T_mid,
                                                                  &zm_input.tmp_qv,
                                                                  &zm_output.tend_out_t,
                                                                  &zm_output.tend_out_s,
                                                                  &zm_output.tend_out_qv,
                                                                  &zm_output.tend_out_u,
                                                                  &zm_output.tend_out_v,
                                                                  &zm_output.tend_tmp_s,
                                                                  &zm_output.tend_tmp_qv,
                                                                  &zm_output.tend_s_snwprd,
                                                                  &zm_output.tend_s_snwevmlt,
                                                                  &zm_output.rain_prod,
                                                                  &zm_output.snow_prod,
                                                                  &zm_output.ntprprd,
                                                                  &zm_output.ntsnprd,
                                                                  &zm_output.flxprec,
                                                                  &zm_output.flxsnow,
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
  // device 3D views on mid-point levels (ncol, nwind, nlev)
  zm_input.tmp_winds = ZMF::uview_3d<Real>(r_mem, m_ncol, nwind, nlev_mid);
  r_mem += zm_input.tmp_winds.size();
  zm_output.tend_tmp_winds = ZMF::uview_3d<Real>(r_mem, m_ncol, nwind, nlev_mid);
  r_mem += zm_output.tend_tmp_winds.size();
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
