#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "radiation.h"
#include "pam_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "pam_feedback.h"
#include "pam_state.h"
#include "pam_radiation.h"
#include "pam_statistics.h"
#include "pam_output.h"
#include "pam_accelerate.h"
#include "pam_variance_transport.h"
#include "pam_hyperdiffusion.h"
#include "sponge_layer.h"
#include "surface_friction.h"
#include "scream_cxx_interface_finalize.h"

// Needed for p3_init
#include "p3_functions.hpp"

#include "pam_debug.h"
bool constexpr enable_check_state = false;


inline int pam_driver_set_subcycle_timestep( pam::PamCoupler &coupler, real crm_dt_fixed ) {
  // calculate the CFL condition and adjust the PAM time loop subcylcing
  //------------------------------------------------------------------------------------------------
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicMax;
  //------------------------------------------------------------------------------------------------
  auto nens          = coupler.get_option<int>("ncrms");
  auto gcm_nlev      = coupler.get_option<int>("gcm_nlev");
  auto crm_nz        = coupler.get_option<int>("crm_nz");
  auto crm_nx        = coupler.get_option<int>("crm_nx");
  auto crm_ny        = coupler.get_option<int>("crm_ny");
  auto crm_dx        = coupler.get_option<real>("crm_dx");
  auto crm_dy        = coupler.get_option<real>("crm_dy");
  auto &dm_device = coupler.get_data_manager_device_readonly();
  auto &dm_host   = coupler.get_data_manager_host_readonly();
  auto uvel       = dm_device.get<real const,4>("uvel");
  auto wvel       = dm_device.get<real const,4>("wvel");
  auto input_zint = dm_host.get<real const,2>("input_zint").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  yakl::ParallelMax<real,yakl::memDevice> pmax( crm_nz*nens );
  real cfl = 0;
  int num_subcycle = 1;
  int constexpr max_num_subcycle = 10;
  real2d wvel_max("wvel_max",crm_nz,nens);
  real2d uvel_max("uvel_max",crm_nz,nens);
  real2d cfl_max("cfl_max",  crm_nz,nens);
  //------------------------------------------------------------------------------------------------
  // initialize max U and W arrays
  parallel_for( SimpleBounds<2>(crm_nz,nens) , YAKL_LAMBDA (int k, int n) {
    wvel_max(k,n) = 0.0;
    uvel_max(k,n) = 0.0;
  });
  // calculate max U and W
  parallel_for( SimpleBounds<4>(crm_nz,crm_ny,crm_nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
    yakl::atomicMax(uvel_max(k,n), std::abs(uvel(k,j,i,n)) );
    yakl::atomicMax(wvel_max(k,n), std::abs(wvel(k,j,i,n)) );
  });
  // find max CFL between horizontal and vertical CFL values
  parallel_for( SimpleBounds<2>(crm_nz,nens) , YAKL_LAMBDA (int k, int n) {
    int k_gcm = gcm_nlev-1-k;
    real crm_dz = input_zint(k_gcm,n) - input_zint(k_gcm+1,n);
    real cfl_u = uvel_max(k,n)*crm_dt_fixed/crm_dx;
    real cfl_w = wvel_max(k,n)*crm_dt_fixed/crm_dz;
    cfl_max(k,n) = std::max(cfl_u,cfl_w);
  });
  // calculate final CFL across ensemble
  real cfl_loc = yakl::intrinsics::maxval(cfl_max);
  cfl = std::max(cfl,cfl_loc);
  // update number of subcycles and time step
  num_subcycle = std::max(num_subcycle,std::max(1,static_cast<int>(ceil(cfl/0.7))));
  real crm_dt_subcycle = crm_dt_fixed / num_subcycle;
  coupler.set_option<real>("crm_dt",crm_dt_subcycle);
  // check for excessive subcylcing - don't exit, just print
  if(num_subcycle > max_num_subcycle) {
    real umax = pmax(uvel_max.data());
    real wmax = pmax(wvel_max.data());
    printf("PAM_DRIVER - WARNING: excessive subcycling!"
          " - num_subcycle: %3.3d  dt: %8.4f  cfl: %8.4f  umax: %8.2f  wmax: %8.2f \n",
          num_subcycle,crm_dt_subcycle,cfl,umax,wmax);
    // exit(-1);
  }

  return num_subcycle;
  //------------------------------------------------------------------------------------------------
}


extern "C" void pam_driver() {
  // This is the primary method for running the PAM CRM in E3SM-MMF.
  //------------------------------------------------------------------------------------------------
  auto &coupler = pam_interface::get_coupler();
  //------------------------------------------------------------------------------------------------
  // retreive coupler options
  auto nens          = coupler.get_option<int>("ncrms");
  auto gcm_nlev      = coupler.get_option<int>("gcm_nlev");
  auto crm_nz        = coupler.get_option<int>("crm_nz");
  auto crm_nx        = coupler.get_option<int>("crm_nx");
  auto crm_ny        = coupler.get_option<int>("crm_ny");
  // auto crm_dx        = coupler.get_option<real>("crm_dx");
  // auto crm_dy        = coupler.get_option<real>("crm_dy");
  auto gcm_dt        = coupler.get_option<real>("gcm_dt");
  auto crm_dt_fixed  = coupler.get_option<real>("crm_dt");
  auto is_first_step = coupler.get_option<bool>("is_first_step");
  auto is_restart    = coupler.get_option<bool>("is_restart");
  bool use_crm_accel = coupler.get_option<bool>("use_crm_accel");
  bool use_MMF_VT    = coupler.get_option<bool>("use_MMF_VT");
  bool enable_physics_tend_stats = coupler.get_option<bool>("enable_physics_tend_stats");
  //------------------------------------------------------------------------------------------------
  // set various coupler options
  coupler.set_option<real>("gcm_physics_dt",gcm_dt);
  coupler.set_option<int>("sponge_num_layers",crm_nz*0.3); // depth of sponge layer
  coupler.set_option<real>("sponge_time_scale",60);        // minimum damping timescale at top
  coupler.set_option<bool>("crm_acceleration_ceaseflag",false);
  //------------------------------------------------------------------------------------------------
  // coupler options for SPAM dycor
  coupler.set_option<int>("crm_dyn_per_phys",1);
  coupler.set_option<bool>("spam_couple_wind_exact_inverse",true);
  coupler.set_option<bool>("spam_clip_negative_densities",true);
  coupler.set_option<bool>("spam_clip_vertical_velocities",true);
  coupler.set_option<bool>("spam_adjust_crm_per_phys_using_vert_cfl",true);
  coupler.set_option<real>("spam_target_cfl",0.7);
  coupler.set_option<real>("spam_max_w",30.0);
  //------------------------------------------------------------------------------------------------
  // Allocate the coupler state and retrieve host/device data managers
  coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

  // set up the grid - this needs to happen before initializing coupler objects
  pam_state_set_grid(coupler);
  //------------------------------------------------------------------------------------------------
  // get seperate data manager objects for host and device
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  // Create objects for dycor, microphysics, and turbulence and initialize them
  bool verbose = is_first_step || is_restart;
  Microphysics micro;
  SGS          sgs;
  Dycore       dycore;
  Radiation    rad;
  micro .init(coupler);
  sgs   .init(coupler);
  dycore.init(coupler,verbose); // pass is_first_step to control verbosity in PAM-C
  rad   .init(coupler);
  //------------------------------------------------------------------------------------------------
  // update coupler GCM state with input GCM state
  pam_state_update_gcm_state(coupler);

  // Copy input CRM state (saved by the GCM) to coupler
  pam_state_copy_input_to_coupler(coupler);

  #ifdef MMF_DISABLE_DENSITY_FORCING
    // update CRM dry density to match GCM and disable dry density forcing
    pam_state_update_dry_density(coupler);
  #endif

  // if debugging - initialize saved state variables and check initial CRM state
  if (enable_check_state) {
    pam_debug_init(coupler);
    pam_debug_check_state(coupler, 0, 0);
  }

  // Compute CRM forcing tendencies
  modules::compute_gcm_forcing_tendencies(coupler);

  // Copy input radiation tendencies to coupler
  pam_radiation_copy_input_to_coupler(coupler);

  // initialize aggregated variables needed for radiation
  pam_radiation_init(coupler);

  // initialize aggregated variables for output statistics
  pam_statistics_init(coupler);

  // initialize variables for CRM mean-state acceleration
  if (use_crm_accel) { pam_accelerate_init(coupler); }

  if (use_MMF_VT) {
    pam_variance_transport_init(coupler);
    pam_variance_transport_compute_forcing(coupler);
  }

  // initilize surface "psuedo-friction" (psuedo => doesn't match "real" GCM friction)
  auto input_tau  = dm_host.get<real const,1>("input_tau00").createDeviceCopy();
  auto input_bflx = dm_host.get<real const,1>("input_bflxls").createDeviceCopy();
  modules::surface_friction_init(coupler, input_tau, input_bflx);

  // Perturb the CRM at the only on first CRM call
  if (is_first_step) {
    auto global_column_id = dm_host.get<int const,1>("global_column_id").createDeviceCopy();
    modules::perturb_temperature( coupler , global_column_id );
  }

  // Microphysics initialization - load lookup tables
  #if defined(P3_CXX)
    if (is_first_step || is_restart) {
      auto am_i_root = coupler.get_option<bool>("am_i_root");
      using P3F = scream::p3::Functions<scream::Real, scream::DefaultDevice>;
      P3F::p3_init(/*write_tables=*/false, am_i_root);
      pam::p3_init_lookup_tables(); // Load P3 lookup table data - avoid re-loading every CRM call
    }
  #endif

  // dycor initialization
  bool do_density_save_recall = false;
  #if defined(MMF_PAM_DYCOR_SPAM)
    // The PAM-C anelastic dycor can dramatically change the dry density due to the assumption that
    // the total density does not change, so we will save it here and recall after the dycor
    do_density_save_recall = true;
    pam_state_set_reference_state(coupler);
    dycore.pre_time_loop(coupler);
  #elif defined(MMF_PAM_DYCOR_AWFL)
    dycore.declare_current_profile_as_hydrostatic(coupler,/*use_gcm_data=*/true);
  #endif

  #ifdef MMF_DISABLE_DENSITY_RECALL
    do_density_save_recall = false;
  #endif
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------

  // set number of CRM steps
  int nstop = int(gcm_dt/crm_dt_fixed);

  // for mean-state acceleration adjust nstop and diagnose horizontal means
  if (use_crm_accel) { 
    pam_accelerate_nstop(coupler,nstop);
    pam_accelerate_diagnose(coupler);
  };

  // Run the CRM
  int nstep = 0;
  while (nstep < nstop) {
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    auto num_subcycle = pam_driver_set_subcycle_timestep(coupler,crm_dt_fixed);
    #if defined(MMF_PAM_DYCOR_SPAM)
      dycore.update_dt(coupler);
    #endif

    if (enable_check_state) { pam_debug_check_state(coupler, 1, nstep); }

    // loop for adaptive subcyling based on CFL
    for(int icycle=1; icycle<=num_subcycle; icycle++) {

      // Apply forcing tendencies
      if (use_MMF_VT) { pam_variance_transport_apply_forcing(coupler); }
      if (enable_check_state) { pam_debug_check_state(coupler, 2, nstep); }
      coupler.run_module( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
      if (enable_check_state) { pam_debug_check_state(coupler, 3, nstep); }
      coupler.run_module( "radiation"                    , [&] (pam::PamCoupler &coupler) {rad   .timeStep(coupler);} );
      if (enable_check_state) { pam_debug_check_state(coupler, 4, nstep); }

      // Dynamics
      if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
      if (do_density_save_recall)    { pam_state_save_dry_density(coupler); }
      coupler.run_module( "dycore", [&] (pam::PamCoupler &coupler) {dycore.timeStep(coupler);} );
      if (do_density_save_recall)    { pam_state_recall_dry_density(coupler); }
      if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"dycor"); }
      if (enable_check_state)        { pam_debug_check_state(coupler, 5, nstep); }

      // Sponge layer damping
      if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
      coupler.run_module( "sponge_layer", modules::sponge_layer );
      if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"sponge"); }

      // Apply hyperdiffusion to account for lack of horizontal mixing in SHOC
      pam_hyperdiffusion(coupler);

      // Turbulence - SHOC
      coupler.run_module( "compute_surface_friction", modules::compute_surface_friction );
      if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
      coupler.run_module( "sgs", [&] (pam::PamCoupler &coupler) {sgs   .timeStep(coupler);} );
      if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"sgs"); }
      if (enable_check_state)        { pam_debug_check_state(coupler, 6, nstep); }

      // Microphysics - P3
      if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
      coupler.run_module( "micro", [&] (pam::PamCoupler &coupler) {micro .timeStep(coupler);} );
      if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"micro"); }
      if (enable_check_state)        { pam_debug_check_state(coupler, 7, nstep); }

    } // num_subcycle

    // CRM mean state acceleration
    if (use_crm_accel && !coupler.get_option<bool>("crm_acceleration_ceaseflag")) {
      pam_accelerate(coupler, nstep, nstop);
      pam_accelerate_diagnose(coupler);
    }

    // Diagnostic aggregation
    pam_radiation_timestep_aggregation(coupler);
    pam_statistics_timestep_aggregation(coupler);

    nstep += 1;
  }
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------

  // Compute CRM feedback tendencies and copy to host
  pam_feedback_compute_tendencies(coupler,gcm_dt);
  pam_feedback_copy_to_host(coupler);
  
  // Copy the final CRM state to the host to be saved for next time step
  pam_state_copy_to_host(coupler);

  // Compute horizontal means of CRM state variables and copy to host
  pam_output_compute_means(coupler);
  pam_output_copy_to_host(coupler);

  // convert aggregated radiation quantities to means and copy to host
  pam_radiation_compute_means(coupler);
  pam_radiation_copy_to_host(coupler);

  // convert aggregated diagnostic quantities to means and copy to host
  pam_statistics_compute_means(coupler);
  pam_statistics_copy_to_host(coupler);

  if (use_MMF_VT) {
    pam_variance_transport_compute_feedback(coupler);
    pam_variance_transport_copy_to_host(coupler);
  }

  //------------------------------------------------------------------------------------------------
  // Finalize and clean up
  micro .finalize(coupler);
  sgs   .finalize(coupler);
  dycore.finalize(coupler);
  rad   .finalize(coupler);
  pam_interface::finalize();
  //------------------------------------------------------------------------------------------------
}

extern "C" void pam_finalize() {
  #if defined(P3_CXX) || defined(SHOC_CXX)
  pam::deallocate_scream_cxx_globals();
  // if using SL tracer advection then COMPOSE will call Kokkos::finalize(), otherwise, call it here
  // pam::call_kokkos_finalize();
  #endif
}
