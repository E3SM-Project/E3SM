#include "pam_coupler.h"
// #include "params.h"
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
#include "sponge_layer.h"
#include "surface_friction.h"
#include "scream_cxx_interface_finalize.h"

#include "pam_hyperdiffusion.h"

// Needed for p3_init
#include "p3_functions.hpp"
#include "p3_f90.hpp"

#include "pam_debug.h"
bool constexpr enable_check_state = false;

extern "C" void pam_driver() {
  //------------------------------------------------------------------------------------------------
  using yakl::intrinsics::abs;
  using yakl::intrinsics::maxval;
  using yakl::atomicAdd;
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &coupler = pam_interface::get_coupler();
  //------------------------------------------------------------------------------------------------
  // retreive coupler options
  auto nens          = coupler.get_option<int>("ncrms");
  auto gcm_nlev      = coupler.get_option<int>("gcm_nlev");
  auto crm_nz        = coupler.get_option<int>("crm_nz");
  auto crm_nx        = coupler.get_option<int>("crm_nx");
  auto crm_ny        = coupler.get_option<int>("crm_ny");
  auto gcm_dt        = coupler.get_option<real>("gcm_dt");
  auto crm_dt        = coupler.get_option<real>("crm_dt");
  auto is_first_step = coupler.get_option<bool>("is_first_step");
  auto is_restart    = coupler.get_option<bool>("is_restart");
  bool use_crm_accel = coupler.get_option<bool>("use_crm_accel");
  bool use_MMF_VT    = coupler.get_option<bool>("use_MMF_VT");
  bool enable_physics_tend_stats = coupler.get_option<bool>("enable_physics_tend_stats");
  //------------------------------------------------------------------------------------------------
  // set various coupler options
  coupler.set_option<real>("gcm_physics_dt",gcm_dt);
  #ifdef MMF_PAM_DPP
  // this is leftover from debugging, but it might still be useful for testing values of crm_per_phys
  coupler.set_option<int>("crm_per_phys",MMF_PAM_DPP);
  #else
  coupler.set_option<int>("crm_per_phys",2);               // # of PAM-C dynamics steps per physics
  #endif
  coupler.set_option<int>("sponge_num_layers",crm_nz*0.3); // depth of sponge layer
  coupler.set_option<real>("sponge_time_scale",60);        // minimum damping timescale at top
  coupler.set_option<bool>("crm_acceleration_ceaseflag",false);
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

  // // update CRM dry density to match GCM and disable dry density forcing
  // pam_state_update_dry_density(coupler);

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
      scream::p3::p3_init(/*write_tables=*/false, am_i_root);
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

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------

  // set number of CRM steps
  int nstop = int(gcm_dt/crm_dt);

  // for mean-state acceleration adjust nstop and diagnose horizontal means
  if (use_crm_accel) { 
    pam_accelerate_nstop(coupler,nstop);
    pam_accelerate_diagnose(coupler);
  };

  // Run the CRM
  real etime_crm = 0;
  int nstep = 0;
  // while (etime_crm < gcm_dt) {
  while (nstep < nstop) {
    if (crm_dt == 0.) { crm_dt = dycore.compute_time_step(coupler); }
    if (etime_crm + crm_dt > gcm_dt) { crm_dt = gcm_dt - etime_crm; }

    if (enable_check_state) { pam_debug_check_state(coupler, 1, nstep); }

    // Apply forcing tendencies
    if (use_MMF_VT) { pam_variance_transport_apply_forcing(coupler); }
    coupler.run_module( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.run_module( "radiation"                    , [&] (pam::PamCoupler &coupler) {rad   .timeStep(coupler);} );
    if (enable_check_state) { pam_debug_check_state(coupler, 2, nstep); }

    // Dynamics
    if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
    if (do_density_save_recall)    { pam_state_save_dry_density(coupler); }
    coupler.run_module( "dycore", [&] (pam::PamCoupler &coupler) {dycore.timeStep(coupler);} );
    if (do_density_save_recall)    { pam_state_recall_dry_density(coupler); }
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"dycor"); }
    if (enable_check_state)        { pam_debug_check_state(coupler, 3, nstep); }

    // Sponge layer damping
    if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
    coupler.run_module( "sponge_layer", modules::sponge_layer );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"sponge"); }
    if (enable_check_state)        { pam_debug_check_state(coupler, 4, nstep); }

    // Apply hyperdiffusion to account for lack of horizontal mixing in SHOC
    pam_hyperdiffusion(coupler);

    // Turbulence - SHOC
    coupler.run_module( "compute_surface_friction", modules::compute_surface_friction );
    if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
    coupler.run_module( "sgs", [&] (pam::PamCoupler &coupler) {sgs   .timeStep(coupler);} );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"sgs"); }
    if (enable_check_state)        { pam_debug_check_state(coupler, 5, nstep); }

    // Microphysics - P3
    if (enable_physics_tend_stats) { pam_statistics_save_state(coupler); }
    coupler.run_module( "micro", [&] (pam::PamCoupler &coupler) {micro .timeStep(coupler);} );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"micro"); }
    if (enable_check_state)        { pam_debug_check_state(coupler, 6, nstep); }

    // CRM mean state acceleration
    if (use_crm_accel && !coupler.get_option<bool>("crm_acceleration_ceaseflag")) {
      pam_accelerate(coupler, nstep, nstop);
      pam_accelerate_diagnose(coupler);
    }

    // Diagnostic aggregation
    pam_radiation_timestep_aggregation(coupler);
    pam_statistics_timestep_aggregation(coupler);

    etime_crm += crm_dt;
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
