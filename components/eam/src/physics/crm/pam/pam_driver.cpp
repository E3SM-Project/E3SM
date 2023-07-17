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
#include "sponge_layer.h"
#include "surface_friction.h"
#include "scream_cxx_interface_finalize.h"

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
  bool enable_physics_tend_stats = coupler.get_option<bool>("enable_physics_tend_stats");
  //------------------------------------------------------------------------------------------------
  // set various coupler options
  coupler.set_option<real>("gcm_physics_dt",gcm_dt);
  coupler.set_option<real>("sponge_time_scale",60);
  coupler.set_option<int>("sponge_num_layers",crm_nz*0.3);
#ifdef CRM_DYN_PER_PHYS
  int cdpp = CRM_DYN_PER_PHYS;
  coupler.set_option<int>("crm_per_phys",cdpp);
#else
  coupler.set_option<int>("crm_per_phys",4);
  if (crm_dt== 1) { coupler.set_option<int>("crm_per_phys",1); }
  if (crm_dt== 5) { coupler.set_option<int>("crm_per_phys",2); }
  if (crm_dt==10) { coupler.set_option<int>("crm_per_phys",4); }
  if (crm_dt==20) { coupler.set_option<int>("crm_per_phys",2); }
  if (crm_dt==30) { coupler.set_option<int>("crm_per_phys",3); }
  if (crm_dt==60) { coupler.set_option<int>("crm_per_phys",6); }
#endif

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
  Microphysics micro;
  SGS          sgs;
  Dycore       dycore;
  Radiation    rad;
  micro .init(coupler);
  sgs   .init(coupler);
  dycore.init(coupler,is_first_step);
  rad   .init(coupler);
  //------------------------------------------------------------------------------------------------
  // update coupler GCM state with input GCM state
  pam_state_update_gcm_state(coupler);

  // Copy input CRM state (saved by the GCM) to coupler
  pam_state_copy_input_to_coupler(coupler);

  // update horizontal mean of CRM dry density to match GCM (also disables dry density forcing)
  // pam_state_update_dry_density(coupler);

  if (enable_check_state) { pam_debug_check_state(coupler, 0, 0); }

  // now that initial state is set, more dycor initialization
  coupler.update_hydrostasis();

  // Compute CRM forcing tendencies
  modules::compute_gcm_forcing_tendencies(coupler);

  // Copy input rad tendencies to coupler
  pam_radiation_copy_input_to_coupler(coupler);

  // initialize rad output variables
  pam_radiation_init(coupler);

  // initialize stat variables
  pam_statistics_init(coupler);

  // initilize surface "psuedo-friction" (psuedo => doesn't match "real" GCM friction)
  auto input_tau  = dm_host.get<real const,1>("input_tau00").createDeviceCopy();
  auto input_bflx = dm_host.get<real const,1>("input_bflxls").createDeviceCopy();
  modules::surface_friction_init(coupler, input_tau, input_bflx);

  if (is_first_step) {
    // Perturb the CRM at the beginning of the run
    auto global_column_id = dm_host.get<int const,1>("global_column_id").createDeviceCopy();
    modules::perturb_temperature( coupler , global_column_id );

    #if defined(P3_CXX)
      auto am_i_root = coupler.get_option<bool>("am_i_root");
      scream::p3::p3_init(/*write_tables=*/false, am_i_root);
      // Load P3 lookup table data to avoid re-loading it every CRM call
      pam::p3_init_lookup_tables();
    #endif
  }

  #if defined(MMF_PAM_DYCOR_SPAM)
    pam_state_set_reference_state(coupler);
    dycore.pre_time_loop(coupler);
  #endif

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------

  // Run the CRM
  real etime_crm = 0;
  int nstep = 0;
  while (etime_crm < gcm_dt) {
    if (crm_dt == 0.) { crm_dt = dycore.compute_time_step(coupler); }
    if (etime_crm + crm_dt > gcm_dt) { crm_dt = gcm_dt - etime_crm; }

    if (enable_check_state) { pam_debug_check_state(coupler, 1, nstep); }

    // run a PAM time step
    coupler.run_module( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.run_module( "radiation"                    , [&] (pam::PamCoupler &coupler) {rad   .timeStep(coupler);} );
    if (enable_check_state) { pam_debug_check_state(coupler, 2, nstep); }

    // anelastic dycor can dramatically change the dry density due to the assumption that the total 
    // density does not change, so we will save it here and recall after the dycor
    pam_state_save_dry_density(coupler);
    pam_statistics_save_state(coupler);
    coupler.run_module( "dycore"                       , [&] (pam::PamCoupler &coupler) {dycore.timeStep(coupler);} );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"dycor"); }
    pam_state_recall_dry_density(coupler);
    if (enable_check_state) { pam_debug_check_state(coupler, 3, nstep); }

    pam_statistics_save_state(coupler);
    coupler.run_module( "sponge_layer"                 , modules::sponge_layer );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"sponge"); }
    if (enable_check_state) { pam_debug_check_state(coupler, 4, nstep); }

    // calculate psuedo friction which will be an input to SHOC
    coupler.run_module( "compute_surface_friction"     , modules::compute_surface_friction );

    pam_statistics_save_state(coupler);
    coupler.run_module( "sgs"                          , [&] (pam::PamCoupler &coupler) {sgs   .timeStep(coupler);} );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"sgs"); }
    if (enable_check_state) { pam_debug_check_state(coupler, 5, nstep); }

    pam_statistics_save_state(coupler);
    coupler.run_module( "micro"                        , [&] (pam::PamCoupler &coupler) {micro .timeStep(coupler);} );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"micro"); }
    if (enable_check_state) { pam_debug_check_state(coupler, 6, nstep); }

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
  pam::call_kokkos_finalize();
  #endif
}
