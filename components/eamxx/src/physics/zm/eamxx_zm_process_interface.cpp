#include "eamxx_config.h" // for SCREAM_CIME_BUILD

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "eamxx_zm_process_interface.hpp"
#include "physics/share/physics_constants.hpp"

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
  const auto comm       = m_grid->get_comm();

  // retrieve local grid parameters
  m_ncol = m_grid->get_num_local_dofs();
  m_nlev = m_grid->get_num_vertical_levels();

  // get max ncol value across ranks to mimic how pcols is used on the fortran side
  m_pcol = m_ncol;
  comm.all_reduce(&m_pcol, 1, MPI_MAX);

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

  // initialize variables on the fortran side
  zm::zm_eamxx_bridge_init(m_pcol, m_nlev);
}

/*------------------------------------------------------------------------------------------------*/
void ZMDeepConvection::run_impl (const double dt)
{
  constexpr int pack_size = Spack::n;
  const int nlevm_packs   = ekat::npack<Spack>(m_nlev);

  // calculate_z_int() contains a team-level parallel_scan, which requires a special policy
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_ncol, nlevm_packs);
  const auto team_policy = TPF::get_default_team_policy(m_ncol, nlevm_packs);

  auto ts_start      = start_of_step_ts();
  bool is_first_step = (ts_start.get_num_steps()==0);

  //----------------------------------------------------------------------------
  // get constants

  const Real cpair  = PC::Cpair;
  const Real latvap = PC::LatVap;

  //----------------------------------------------------------------------------
  // get fields

  // variables not updated by ZM
  const auto& phis        = get_field_in("phis")          .get_view<const Real*>();
  const auto& p_mid       = get_field_in("p_mid")         .get_view<const Spack**, Host>();
  const auto& p_int       = get_field_in("p_int")         .get_view<const Spack**, Host>();
  const auto& p_del       = get_field_in("pseudo_density").get_view<const Spack**, Host>();
  const auto& omega       = get_field_in("omega")         .get_view<const Spack**, Host>();
  const auto& cldfrac     = get_field_in("cldfrac_tot")   .get_view<const Spack**, Host>();
  const auto& pblh        = get_field_in("pbl_height")    .get_view<const Real*>();
  const auto& landfrac    = get_field_in("landfrac")      .get_view<const Real*>();
  const auto& thl_sec     = get_field_in("thl_sec")       .get_view<const Spack**, Host>();
  const auto& qc          = get_field_in("qc")            .get_view<const Spack**, Host>();

  // variables updated by ZM
  const auto& T_mid       = get_field_out("T_mid")        .get_view<Spack**, Host>();
  const auto& qv          = get_field_out("qv")           .get_view<Spack**, Host>();
  const auto& hwinds_fld  = get_field_out("horiz_winds");
  const auto& uwind       = hwinds_fld.get_component(0)   .get_view<Spack**, Host>();
  const auto& vwind       = hwinds_fld.get_component(1)   .get_view<Spack**, Host>();

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

  // initialize output buffer variables
  zm_output.init(m_pcol, m_nlev);

  //----------------------------------------------------------------------------
  // calculate altitude on interfaces (z_int) and mid-points (z_mid)

  // const auto zm_input_loc = zm_input;
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    const auto z_mid_i = ekat::subview(zm_input.z_mid, i);
    const auto z_del_i = ekat::subview(zm_input.z_del, i);
    const auto z_int_i = ekat::subview(zm_input.z_int, i);
    const auto p_mid_i = ekat::subview(zm_input.p_mid, i);
    const auto p_del_i = ekat::subview(zm_input.p_del, i);
    const auto T_mid_i = ekat::subview(zm_input.T_mid, i);
    const auto qv_i    = ekat::subview(zm_input.qv,    i);
    auto z_surf = 0.0; // ZM expects z_mid & z_int to be altitude above the surface
    PF::calculate_dz(team, p_del_i, p_mid_i, T_mid_i, qv_i, z_del_i);
    team.team_barrier();
    PF::calculate_z_int(team, m_nlev, z_del_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, m_nlev, z_int_i, z_mid_i);
    team.team_barrier();
  });

  //----------------------------------------------------------------------------
  // calculate temperature perturbation from SHOC thetal varance for ZM parcel/CAPE

  Kokkos::parallel_for("zm_calculate_tpert",m_ncol, KOKKOS_LAMBDA (const int i) {
    if (is_first_step) {
      zm_input.tpert(i) = 0.0;
    } else {
      // identify interface index for top of PBL
      int pblh_k_ind = -1;
      for (int k=0; k<m_nlev; ++k) {
        auto z_int_tmp_k   = zm_input.z_int(i,k/pack_size)[k%pack_size];
        auto z_int_tmp_kp1 = zm_input.z_int(i,k/pack_size)[k%pack_size];
        if ( z_int_tmp_k>pblh(i) && z_int_tmp_kp1<=pblh(i) ) {
          pblh_k_ind = k;
        }
      }
      if (pblh_k_ind==-1) {
        // PBL top index not found, so just set the perturbation to zero
        zm_input.tpert(i) = 0.0;
      } else {
        // calculate tpert as std deviation of temperature from SHOC's theta-l variance
        auto exner_pbl    = PF::exner_function( p_mid(i,pblh_k_ind/pack_size)[pblh_k_ind%pack_size] );
        auto qc_pbl       = qc(i,pblh_k_ind/pack_size)[pblh_k_ind%pack_size];
        auto thl_sec_pbl  = thl_sec(i,pblh_k_ind/pack_size)[pblh_k_ind%pack_size];
        auto thl_std_pbl  = sqrt( thl_sec_pbl ); // std deviation of thetal;
        zm_input.tpert(i) = ( thl_std_pbl + (latvap/cpair)*qc_pbl ) / exner_pbl;
        zm_input.tpert(i) = std::min(2.0,zm_input.tpert(i)); // apply limiter
      }
    }
  });

  //----------------------------------------------------------------------------
  // run the ZM scheme

  zm_eamxx_bridge_run(m_ncol, m_nlev, zm_input, zm_output, zm_opts);

  //----------------------------------------------------------------------------
  // update prognostic fields

  const auto zm_output_loc = zm_output; // create temp to avoid "Implicit capture" warning

  if (zm_opts.apply_tendencies) {
    // accumulate surface precipitation fluxes
    Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol), KOKKOS_LAMBDA (const int i) {
      auto prec_liq = zm_output_loc.prec(i) - zm_output_loc.snow(i);
      precip_liq_surf_mass(i) += std::max(prec_liq,0.0) * PC::RHO_H2O * dt;
      precip_ice_surf_mass(i) += zm_output_loc.snow(i)  * PC::RHO_H2O * dt;
    });
    // update 3D prognostic variables
    Kokkos::parallel_for("zm_update_prognostics",team_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevm_packs), [&](const int k) {
        T_mid(i,k) += zm_output_loc.tend_s (i,k)/cpair * dt;
        qv   (i,k) += zm_output_loc.tend_qv(i,k)       * dt;
        uwind(i,k) += zm_output_loc.tend_u (i,k)       * dt;
        vwind(i,k) += zm_output_loc.tend_v (i,k)       * dt;
      });
    });
  }

  //----------------------------------------------------------------------------
  // Update output fields

  // NOTE - in the future we might want to clean this up using Kokkos::deep_copy(),
  // but this is currently not possible due to the pcol/ncol thing for the fortran bridge

  // 2D output (no vertical dimension)
  const auto& zm_prec       = get_field_out("zm_prec")        .get_view<Real*>();
  const auto& zm_snow       = get_field_out("zm_snow")        .get_view<Real*>();
  const auto& zm_cape       = get_field_out("zm_cape")        .get_view<Real*>();
  const auto& zm_activity   = get_field_out("zm_activity")    .get_view<Real*>();
  Kokkos::parallel_for("zm_diag_outputs_2D",m_ncol, KOKKOS_LAMBDA (const int i) {
    zm_prec    (i) = zm_output_loc.prec    (i);
    zm_snow    (i) = zm_output_loc.snow    (i);
    zm_cape    (i) = zm_output_loc.cape    (i);
    zm_activity(i) = zm_output_loc.activity(i);
  });

  // 3D output (vertically resolved)
  const auto& zm_T_mid_tend = get_field_out("zm_T_mid_tend")  .get_view<Spack**, Host>();
  const auto& zm_qv_tend    = get_field_out("zm_qv_tend")     .get_view<Spack**, Host>();
  const auto& zm_u_tend     = get_field_out("zm_u_tend")      .get_view<Spack**, Host>();
  const auto& zm_v_tend     = get_field_out("zm_v_tend")      .get_view<Spack**, Host>();
  Kokkos::parallel_for("zm_diag_outputs",team_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const auto i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevm_packs), [&](const int k) {
      zm_T_mid_tend(i,k) = zm_output_loc.tend_s (i,k)/cpair;
      zm_qv_tend   (i,k) = zm_output_loc.tend_qv(i,k);
      zm_u_tend    (i,k) = zm_output_loc.tend_u (i,k);
      zm_v_tend    (i,k) = zm_output_loc.tend_v (i,k);
    });
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
  const int nlevm_packs = ekat::npack<Spack>(m_nlev);
  const int nlevi_packs = ekat::npack<Spack>(m_nlev+1);
  size_t zm_buffer_size = 0;

  zm_buffer_size+= ZMF::zm_input_state::num_1d_intgr_views   * sizeof(Int)   * m_pcol;
  zm_buffer_size+= ZMF::zm_input_state::num_1d_scalr_views   * sizeof(Scalar)* m_pcol;

  zm_buffer_size+= ZMF::zm_input_state::num_2d_midlv_c_views * sizeof(Spack) * m_pcol * nlevm_packs;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_intfc_c_views * sizeof(Spack) * m_pcol * nlevi_packs;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_midlv_f_views * sizeof(Real)  * m_pcol * m_nlev;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_intfc_f_views * sizeof(Real)  * m_pcol * (m_nlev+1);

  zm_buffer_size+= ZMF::zm_output_tend::num_1d_scalr_views   * sizeof(Scalar)* m_pcol;
  zm_buffer_size+= ZMF::zm_output_tend::num_1d_intgr_views   * sizeof(Int)   * m_pcol;

  zm_buffer_size+= ZMF::zm_output_tend::num_2d_midlv_c_views * sizeof(Spack) * m_pcol * nlevm_packs;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_intfc_c_views * sizeof(Spack) * m_pcol * nlevi_packs;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_midlv_f_views * sizeof(Real)  * m_pcol * m_nlev;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_intfc_f_views * sizeof(Real)  * m_pcol * (m_nlev+1);

  return zm_buffer_size;
}

/*------------------------------------------------------------------------------------------------*/

void ZMDeepConvection::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");

  const int nlevm_packs = ekat::npack<Spack>(m_nlev);
  const int nlevi_packs = ekat::npack<Spack>(m_nlev+1);

  constexpr auto num_1d_intgr_views   = ZMF::zm_input_state::num_1d_intgr_views   + ZMF::zm_output_tend::num_1d_intgr_views;
  constexpr auto num_1d_scalr_views   = ZMF::zm_input_state::num_1d_scalr_views   + ZMF::zm_output_tend::num_1d_scalr_views;
  constexpr auto num_2d_midlv_c_views = ZMF::zm_input_state::num_2d_midlv_c_views + ZMF::zm_output_tend::num_2d_midlv_c_views;
  constexpr auto num_2d_intfc_c_views = ZMF::zm_input_state::num_2d_intfc_c_views + ZMF::zm_output_tend::num_2d_intfc_c_views;
  constexpr auto num_2d_midlv_f_views = ZMF::zm_input_state::num_2d_midlv_f_views + ZMF::zm_output_tend::num_2d_midlv_f_views;
  constexpr auto num_2d_intfc_f_views = ZMF::zm_input_state::num_2d_intfc_f_views + ZMF::zm_output_tend::num_2d_intfc_f_views;

  //----------------------------------------------------------------------------
  Int* i_mem = reinterpret_cast<Int*>(buffer_manager.get_memory());
  //----------------------------------------------------------------------------
  // 1D integer variables
  ZMF::uview_1d<Int>* int_ptrs[num_1d_intgr_views]    = { &zm_output.activity
                                                        };
  for (int i=0; i<num_1d_intgr_views; ++i) {
    *int_ptrs[i] = ZMF::uview_1d<Int>(i_mem, m_pcol);
    i_mem += int_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Scalar* scl_mem = reinterpret_cast<Scalar*>(i_mem);
  //----------------------------------------------------------------------------
  // 1D scalar variables
  ZMF::uview_1d<Scalar>* scl_ptrs[num_1d_scalr_views] = { &zm_input.tpert,
                                                          &zm_output.prec,
                                                          &zm_output.snow,
                                                          &zm_output.cape
                                                        };
  for (int i=0; i<num_1d_scalr_views; ++i) {
    *scl_ptrs[i] = ZMF::uview_1d<Scalar>(scl_mem, m_pcol);
    scl_mem += scl_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Real* r_mem = reinterpret_cast<Real*>(scl_mem);
  //----------------------------------------------------------------------------
  // 2D "f_" views on mid-point levels
  ZMF::uview_2dl<Real>* midlv_f_ptrs[num_2d_midlv_f_views]  = { &zm_input.f_z_mid,
                                                                &zm_input.f_p_mid,
                                                                &zm_input.f_p_del,
                                                                &zm_input.f_T_mid,
                                                                &zm_input.f_qv,
                                                                &zm_input.f_uwind,
                                                                &zm_input.f_vwind,
                                                                &zm_input.f_omega,
                                                                &zm_input.f_cldfrac,
                                                                &zm_output.f_tend_s,
                                                                &zm_output.f_tend_qv,
                                                                &zm_output.f_tend_u,
                                                                &zm_output.f_tend_v,
                                                                &zm_output.f_rain_prod,
                                                                &zm_output.f_snow_prod
                                                              };
  for (int i=0; i<num_2d_midlv_f_views; ++i) {
    *midlv_f_ptrs[i] = ZMF::uview_2dl<Real>(r_mem, m_pcol, m_nlev);
    r_mem += midlv_f_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  // 2D "f_" views on interface levels
  ZMF::uview_2dl<Real>* intfc_f_ptrs[num_2d_intfc_f_views]  = { &zm_input.f_z_int,
                                                                &zm_input.f_p_int,
                                                                &zm_output.f_prec_flux,
                                                                &zm_output.f_snow_flux,
                                                                &zm_output.f_mass_flux
                                                              };
  for (int i=0; i<num_2d_intfc_f_views; ++i) {
    *intfc_f_ptrs[i] = ZMF::uview_2dl<Real>(r_mem, m_pcol, (m_nlev+1));
    r_mem += intfc_f_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Spack* spk_mem = reinterpret_cast<Spack*>(r_mem);
  //----------------------------------------------------------------------------
  // 2D views on mid-point levels
  ZMF::uview_2d<Spack>* midlv_c_ptrs[num_2d_midlv_c_views]  = { &zm_input.z_mid,
                                                                &zm_input.z_del,
                                                                &zm_output.tend_s,
                                                                &zm_output.tend_qv,
                                                                &zm_output.tend_u,
                                                                &zm_output.tend_v,
                                                                &zm_output.rain_prod,
                                                                &zm_output.snow_prod
                                                              };
  for (int i=0; i<num_2d_midlv_c_views; ++i) {
    *midlv_c_ptrs[i] = ZMF::uview_2d<Spack>(spk_mem, m_pcol, nlevm_packs);
    spk_mem += midlv_c_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  // 2D variables on interface levels
  ZMF::uview_2d<Spack>* intfc_c_ptrs[num_2d_intfc_c_views]  = { &zm_input.z_int,
                                                                &zm_output.prec_flux,
                                                                &zm_output.snow_flux,
                                                                &zm_output.mass_flux
                                                              };
  for (int i=0; i<num_2d_intfc_c_views; ++i) {
    *intfc_c_ptrs[i] = ZMF::uview_2d<Spack>(spk_mem, m_pcol, nlevi_packs);
    spk_mem += intfc_c_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Real* total_mem = reinterpret_cast<Real*>(spk_mem);
  size_t used_mem = (reinterpret_cast<Real*>(total_mem) - buffer_manager.get_memory())*sizeof(Real);
  auto mem_chk = ( used_mem == requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory != requested memory for ZMDeepConvection.");
}

/*------------------------------------------------------------------------------------------------*/


} // namespace scream
