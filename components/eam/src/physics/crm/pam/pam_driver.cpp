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
  coupler.set_option<int>("crm_per_phys",4);
  if (crm_dt==10) { coupler.set_option<int>("crm_per_phys",10); }
  if (crm_dt==20) { coupler.set_option<int>("crm_per_phys",2); }
  if (crm_dt==30) { coupler.set_option<int>("crm_per_phys",3); }
  if (crm_dt==60) { coupler.set_option<int>("crm_per_phys",6); }
  //------------------------------------------------------------------------------------------------
  // Allocate the coupler state and retrieve host/device data managers
  coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

  // set up the grid - this needs to happen before initializing coupler objects
  pam_state_set_grid(coupler);
  //------------------------------------------------------------------------------------------------
  // get seperate data manager objects for host and device
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();

  // temporary variable for saving and recalling dry density in pam_state
  dm_device.register_and_allocate<real>("density_dry_save", "temporary CRM dry density", {crm_nz,crm_ny,crm_nx,nens}, {"z","y","x","nens"} );
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
    auto gcolp = dm_host.get<int const,1>("gcolp").createDeviceCopy();
    modules::perturb_temperature( coupler , gcolp );

    #if defined(P3_CXX)
      auto am_i_root = coupler.get_option<bool>("am_i_root");
      scream::p3::p3_init(/*write_tables=*/false, am_i_root);

      // Load P3 lookup table data to avoid re-loading it every CRM call
      pam::p3_init_lookup_tables();
    #endif
  }

  #if defined(MMF_PAM_DYCOR_SPAM)
    // set anelastic reference state using GCM state
    pam_state_set_reference_state(coupler);
    dycore.pre_time_loop(coupler);
  #endif

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------

  real4d save_temp("save_temp",crm_nz,crm_ny,crm_nx,nens);
  real4d save_rhod("save_rhod",crm_nz,crm_ny,crm_nx,nens);
  real4d save_rhov("save_rhov",crm_nz,crm_ny,crm_nx,nens);
  real4d save_rhoc("save_rhoc",crm_nz,crm_ny,crm_nx,nens);
  real4d save_rhoi("save_rhoi",crm_nz,crm_ny,crm_nx,nens);
  real4d save_uvel("save_uvel",crm_nz,crm_ny,crm_nx,nens);
  real4d save_wvel("save_wvel",crm_nz,crm_ny,crm_nx,nens);

  auto ref_rho_d  = dm_device.get<real const,2>("ref_density_dry");
  auto ref_rho_v  = dm_device.get<real const,2>("ref_density_vapor");
  auto ref_rho_c  = dm_device.get<real const,2>("ref_density_liq");
  auto ref_rho_i  = dm_device.get<real const,2>("ref_density_ice");
  auto ref_temp   = dm_device.get<real const,2>("ref_temp");
  auto gcm_rho_d  = dm_device.get<real const,2>("gcm_density_dry");
  auto gcm_temp   = dm_device.get<real const,2>("gcm_temp"       );
  auto gcm_rho_v  = dm_device.get<real const,2>("gcm_water_vapor");
  auto gcm_rho_c  = dm_device.get<real const,2>("gcm_cloud_water");
  auto gcm_rho_i  = dm_device.get<real const,2>("gcm_cloud_ice"  );

  auto temp = dm_device.get<real,4>("temp");
  auto rhod = dm_device.get<real,4>("density_dry");
  auto rhov = dm_device.get<real,4>("water_vapor");
  auto rhoc = dm_device.get<real,4>("cloud_water");
  auto rhoi = dm_device.get<real,4>("ice");
  auto zmid = dm_device.get<real,2>("vertical_midpoint_height" );
  auto zint = dm_device.get<real,2>("vertical_interface_height");
  auto uvel = dm_device.get<real,4>("uvel");
  auto wvel = dm_device.get<real,4>("wvel");


  // real grav = 9.80616;
  // auto input_phis = dm_host.get<real const,1>("input_phis").createDeviceCopy();
  // auto lat        = dm_host.get<real const,1>("latitude"  ).createDeviceCopy();
  // auto lon        = dm_host.get<real const,1>("longitude" ).createDeviceCopy();

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

    // if (enable_check_state) { 
      // save stuff for after-dycor check
      parallel_for("", SimpleBounds<4>(crm_nz,crm_ny,crm_nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
        save_temp(k,j,i,iens) = temp(k,j,i,iens);
        save_rhod(k,j,i,iens) = rhod(k,j,i,iens);
        save_rhov(k,j,i,iens) = rhov(k,j,i,iens);
        save_rhoc(k,j,i,iens) = rhoc(k,j,i,iens);
        save_rhoi(k,j,i,iens) = rhoi(k,j,i,iens);
        save_uvel(k,j,i,iens) = uvel(k,j,i,iens);
        save_wvel(k,j,i,iens) = wvel(k,j,i,iens);
      });
    // }

    // anelastic dycor messes up the dry density ddue to how the GCM forcing changes total density
    // so we need to save it here and recall after the dycor
    pam_state_save_dry_density(coupler);

    pam_statistics_save_state(coupler);
    coupler.run_module( "dycore"                       , [&] (pam::PamCoupler &coupler) {dycore.timeStep(coupler);} );
    if (enable_physics_tend_stats) { pam_statistics_aggregate_tendency(coupler,"dycor"); }

    // recall the dry density to the values before the dycor
    pam_state_recall_dry_density(coupler);

    // Check stuff
    // if (enable_check_state) { 
      parallel_for("", SimpleBounds<3>(crm_ny,crm_nx,nens), YAKL_LAMBDA (int j, int i, int iens) {
        int  tmp_i = -1;
        int  tmp_j = -1;
        int  tmp_n = -1;
        // bool problem_found = false;
        for (int k=0; k<crm_nz; k++) {
          if ( temp(k,j,i,iens)<0 || isnan(temp(k,j,i,iens)) ) {
            // problem_found = true;
            tmp_i = i;
            tmp_j = j;
            tmp_n = iens;
            real dz = (zint(k+1,iens) - zint(k,iens));
            // printf("WHDEBUG - st:%3.3d n:%3.3d i:%3.3d k:%3.3d z:%8.2g dz:%8.2g =curr= t: %8.2g rv: %8.2g rc: %8.2g ri: %8.2g rd: %8.2g u: %8.2g w: %8.2g =prev= t: %8.2g rv: %8.2g rc: %8.2g ri: %8.2g u: %8.2g w: %8.2g \n",
            //   nstep,iens,i,k,
            //   zmid(k,iens),
            //   dz,
            //   temp(k,j,i,iens),
            //   rhov(k,j,i,iens),
            //   rhoc(k,j,i,iens),
            //   rhoi(k,j,i,iens),
            //   rhod(k,j,i,iens),
            //   uvel(k,j,i,iens),
            //   wvel(k,j,i,iens),
            //   save_temp(k,j,i,iens),
            //   save_rhov(k,j,i,iens),
            //   save_rhoc(k,j,i,iens),
            //   save_rhoi(k,j,i,iens),
            //   save_uvel(k,j,i,iens),
            //   save_wvel(k,j,i,iens)
            // );
            printf("WHDEBUG - st:%3.3d n:%3.3d i:%3.3d k:%3.3d z:%8.2g dz:%8.2g =curr= t: %8.2g u: %8.2g w: %8.2g =prev= t: %8.2g u: %8.2g w: %8.2g \n",
              nstep,iens,i,k,
              zmid(k,iens),
              dz,
              temp(k,j,i,iens),
              uvel(k,j,i,iens),
              wvel(k,j,i,iens),
              save_temp(k,j,i,iens),
              save_uvel(k,j,i,iens),
              save_wvel(k,j,i,iens)
            );
          }
        }
      });
    // }

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
