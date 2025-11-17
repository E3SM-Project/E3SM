#include "eamxx_config.h" // for SCREAM_CIME_BUILD

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
void ZMDeepConvection::set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  constexpr int pack_size = Spack::n;

  // Gather runtime options from file
  zm_opts.load_runtime_options(m_params);

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
  add_field<Required>("cldfrac_tot",          scalar3d_mid, nondim, grid_name, pack_size);
  add_field<Required>("pbl_height",           scalar2d    , m,      grid_name);
  add_field<Required>("landfrac",             scalar2d    , nondim, grid_name);
  add_field<Required>("thl_sec",              scalar3d_int, K2,     grid_name, pack_size); // thetal variance for PBL temperature perturbation
  add_tracer<Required>("qc",                  m_grid,       kg/kg,             pack_size);

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,      grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,             pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,    grid_name, pack_size);

  // Output variables
  add_field <Updated>("precip_liq_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");
  add_field <Updated>("precip_ice_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");

  // Diagnostic Outputs
  add_field<Computed>("zm_prec",              scalar2d,     m/s,    grid_name);
  add_field<Computed>("zm_snow",              scalar2d,     m/s,    grid_name);
  add_field<Computed>("zm_cape",              scalar2d,     J/kg,   grid_name);
  add_field<Computed>("zm_activity",          scalar2d,     nondim, grid_name);

  add_field<Computed>("zm_T_mid_tend",        scalar3d_mid, K/s,    grid_name, pack_size);
  add_field<Computed>("zm_qv_tend",           scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("zm_u_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
  add_field<Computed>("zm_v_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);

}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::initialize_impl (const RunType)
{
  // Set property checks for fields in this process
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("qv"),m_grid,1e-13,0.2,true);

  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_liq_surf_mass"),m_grid,0.0,false);
  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_ice_surf_mass"),m_grid,0.0,false);

  //----------------------------------------------------------------------------
  // allocate host mirror variables

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

  //----------------------------------------------------------------------------
  // initialize variables on the fortran side
  zm::zm_eamxx_bridge_init(m_nlev);

}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::run_impl (const double dt)
{
  const int nlev_mid_packs   = ekat::npack<Spack>(m_nlev);

  // calculate_z_int() contains a team-level parallel_scan, which requires a special policy
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_ncol, nlev_mid_packs);

  auto ts_start      = start_of_step_ts();
  bool is_first_step = (ts_start.get_num_steps()==0);

  //----------------------------------------------------------------------------
  // get fields

  // variables not updated by ZM
  const auto& phis        = get_field_in("phis")          .get_view<const Real*>();
  const auto& p_mid       = get_field_in("p_mid")         .get_view<const Spack**>();
  const auto& p_int       = get_field_in("p_int")         .get_view<const Spack**>();
  const auto& p_del       = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& omega       = get_field_in("omega")         .get_view<const Spack**>();
  const auto& cldfrac     = get_field_in("cldfrac_tot")   .get_view<const Spack**>();
  const auto& pblh        = get_field_in("pbl_height")    .get_view<const Real*>();
  const auto& landfrac    = get_field_in("landfrac")      .get_view<const Real*>();
  const auto& thl_sec     = get_field_in("thl_sec")       .get_view<const Spack**>();
  const auto& qc          = get_field_in("qc")            .get_view<const Spack**>();

  // variables updated by ZM
  const auto& T_mid       = get_field_out("T_mid")        .get_view<Spack**>();
  const auto& qv          = get_field_out("qv")           .get_view<Spack**>();
  const auto& hwinds_fld  = get_field_out("horiz_winds");
  const auto& uwind       = hwinds_fld.get_component(0)   .get_view<Spack**>();
  const auto& vwind       = hwinds_fld.get_component(1)   .get_view<Spack**>();

  const auto& precip_liq_surf_mass = get_field_out("precip_liq_surf_mass").get_view<Real*>();
  const auto& precip_ice_surf_mass = get_field_out("precip_ice_surf_mass").get_view<Real*>();

  //----------------------------------------------------------------------------
  // prepare input struct

  zm_input.dtime          = dt;
  zm_input.is_first_step  = is_first_step;
  zm_input.phis           = phis;
  zm_input.p_mid          = p_mid;
  zm_input.p_int          = p_int;
  zm_input.p_del          = p_del;
  zm_input.T_mid          = T_mid;
  zm_input.qv             = qv;
  zm_input.uwind          = uwind;
  zm_input.vwind          = vwind;
  zm_input.omega          = omega;
  zm_input.cldfrac        = cldfrac;
  zm_input.pblh           = pblh;
  zm_input.landfrac       = landfrac;
  zm_input.thl_sec        = thl_sec;
  zm_input.qc             = qc;

  // initialize output buffer variables
  zm_output.init(m_ncol, m_nlev);

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
  auto loc_nlev = m_nlev;

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
    PF::calculate_z_int(team, loc_nlev, z_del_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, loc_nlev, z_int_i, z_mid_i);
    team.team_barrier();
  });

  //----------------------------------------------------------------------------
  // calculate temperature perturbation from SHOC thetal varance for ZM parcel/CAPE

  zm_input.calculate_tpert(m_ncol,m_nlev,is_first_step);

  //----------------------------------------------------------------------------
  // run the ZM scheme

  zm_eamxx_bridge_run(m_ncol, m_nlev, zm_input, zm_output, zm_opts);

  //----------------------------------------------------------------------------
  // create temporaries of output variables to avoid "Implicit capture" warning

  const auto loc_zm_output_prec    = zm_output.prec;
  const auto loc_zm_output_snow    = zm_output.snow;
  const auto loc_zm_output_cape        = zm_output.cape;
  const auto loc_zm_output_activity    = zm_output.activity;
  const auto loc_zm_output_tend_t  = zm_output.tend_t;
  const auto loc_zm_output_tend_qv = zm_output.tend_qv;
  const auto loc_zm_output_tend_u  = zm_output.tend_u;
  const auto loc_zm_output_tend_v  = zm_output.tend_v;

  //----------------------------------------------------------------------------
  // update prognostic fields

  if (zm_opts.apply_tendencies) {
    // accumulate surface precipitation fluxes
    Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol), KOKKOS_LAMBDA (const int i) {
      auto prec_liq = loc_zm_output_prec(i) - loc_zm_output_snow(i);
      precip_liq_surf_mass(i) += ekat::impl::max(0.0,prec_liq) * PC::RHO_H2O * dt;
      precip_ice_surf_mass(i) += loc_zm_output_snow(i) * PC::RHO_H2O * dt;
    });

    Kokkos::parallel_for("zm_update_prognostic",KT::RangePolicy(0, m_ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int idx) {
      const int i = idx/nlev_mid_packs;
      const int k = idx%nlev_mid_packs;
      T_mid(i,k) += loc_zm_output_tend_t (i,k) * dt;
      qv   (i,k) += loc_zm_output_tend_qv(i,k) * dt;
      uwind(i,k) += loc_zm_output_tend_u (i,k) * dt;
      vwind(i,k) += loc_zm_output_tend_v (i,k) * dt;
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
  const auto& zm_T_mid_tend = get_field_out("zm_T_mid_tend")  .get_view<Spack**>();
  const auto& zm_qv_tend    = get_field_out("zm_qv_tend")     .get_view<Spack**>();
  const auto& zm_u_tend     = get_field_out("zm_u_tend")      .get_view<Spack**>();
  const auto& zm_v_tend     = get_field_out("zm_v_tend")      .get_view<Spack**>();
  Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_mid_packs;
    const int k = idx%nlev_mid_packs;
    zm_T_mid_tend(i,k) = loc_zm_output_tend_t (i,k);
    zm_qv_tend   (i,k) = loc_zm_output_tend_qv(i,k);
    zm_u_tend    (i,k) = loc_zm_output_tend_u (i,k);
    zm_v_tend    (i,k) = loc_zm_output_tend_v (i,k);
  });

}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::finalize_impl ()
{
  // placeholder for final cleanup
}

/*------------------------------------------------------------------------------------------------*/

size_t ZMDeepConvection::requested_buffer_size_in_bytes() const
{
  const int nlev_mid_packs = ekat::npack<Spack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Spack>(m_nlev+1);
  size_t zm_buffer_size = 0;

  zm_buffer_size+= ZMF::zm_input_state::num_1d_intgr * sizeof(Int)   * m_ncol;
  zm_buffer_size+= ZMF::zm_input_state::num_1d_scalr * sizeof(Scalar)* m_ncol;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_midlv * sizeof(Spack) * m_ncol * nlev_mid_packs;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_intfc * sizeof(Spack) * m_ncol * nlev_int_packs;

  zm_buffer_size+= ZMF::zm_output_tend::num_1d_intgr * sizeof(Int)   * m_ncol;
  zm_buffer_size+= ZMF::zm_output_tend::num_1d_scalr * sizeof(Scalar)* m_ncol;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_midlv * sizeof(Spack) * m_ncol * nlev_mid_packs;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_intfc * sizeof(Spack) * m_ncol * nlev_int_packs;

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

  const int nlev_mid_packs = ekat::npack<Spack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Spack>(m_nlev+1);

  constexpr auto num_1d_intgr = ZMF::zm_input_state::num_1d_intgr + ZMF::zm_output_tend::num_1d_intgr;
  constexpr auto num_1d_scalr = ZMF::zm_input_state::num_1d_scalr + ZMF::zm_output_tend::num_1d_scalr;
  constexpr auto num_2d_midlv = ZMF::zm_input_state::num_2d_midlv + ZMF::zm_output_tend::num_2d_midlv;
  constexpr auto num_2d_intfc = ZMF::zm_input_state::num_2d_intfc + ZMF::zm_output_tend::num_2d_intfc;

  constexpr int num_f_mid  = (9+6);
  constexpr int num_f_int  = (2+3);

  //----------------------------------------------------------------------------
  Int* i_mem = reinterpret_cast<Int*>(buffer_manager.get_memory());
  //----------------------------------------------------------------------------
  // device 1D integer variables
  ZMF::uview_1d<Int>* ptrs_1d_intgr[num_1d_intgr]             = { &zm_output.activity };
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
  Spack* spk_mem = reinterpret_cast<Spack*>(r_mem);
  // ***************************************************************************
  // TEMPORARY
  // ***************************************************************************
  // Spack* spk_mem = reinterpret_cast<Spack*>(scl_mem);
  //----------------------------------------------------------------------------
  // device 2D views on mid-point levels
  ZMF::uview_2d<Spack>* ptrs_2d_midlv[num_2d_midlv]           = { &zm_input.z_mid,
                                                                  &zm_input.z_del,
                                                                  &zm_output.tend_t,
                                                                  &zm_output.tend_qv,
                                                                  &zm_output.tend_u,
                                                                  &zm_output.tend_v,
                                                                  &zm_output.rain_prod,
                                                                  &zm_output.snow_prod,
                                                                };
  for (auto& v : ptrs_2d_midlv) {
    *v = ZMF::uview_2d<Spack>(spk_mem, m_ncol, nlev_mid_packs);
    spk_mem += v->size();
  }
  //----------------------------------------------------------------------------
  // device 2D views on interface levels
  ZMF::uview_2d<Spack>* ptrs_2d_intfc[num_2d_intfc]           = { &zm_input.z_int,
                                                                  &zm_output.prec_flux,
                                                                  &zm_output.snow_flux,
                                                                  &zm_output.mass_flux,
                                                                };
  for (auto& v : ptrs_2d_intfc) {
    *v = ZMF::uview_2d<Spack>(spk_mem, m_ncol, nlev_int_packs);
    spk_mem += v->size();
  }
  //----------------------------------------------------------------------------
  Real* total_mem = reinterpret_cast<Real*>(spk_mem);
  size_t used_mem = (reinterpret_cast<Real*>(total_mem) - buffer_manager.get_memory())*sizeof(Real);
  auto mem_chk = ( used_mem == requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory != requested memory for ZMDeepConvection.");
  //----------------------------------------------------------------------------
}

/*------------------------------------------------------------------------------------------------*/


} // namespace scream
