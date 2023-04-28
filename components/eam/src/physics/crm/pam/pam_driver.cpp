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
#include "sponge_layer.h"
#include "surface_friction.h"
#include "scream_cxx_interface_finalize.h"

void nan_chk(pam::PamCoupler &coupler, std::string id, std::string var_name) {
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nx   = coupler.get_option<int>("crm_nx");
  auto ny   = coupler.get_option<int>("crm_ny");
  auto nz   = coupler.get_option<int>("crm_nz");
  auto nens = coupler.get_option<int>("ncrms");
  auto tvar = dm_device.get<real,4>(var_name);
  parallel_for("Horz mean of CRM dry density", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    if ( isnan(tvar(k,j,i,iens)) ) {
      printf("  WHDEBUG - NaN detected in driver - k:%d  i:%d  e:%d  %s:%d \n",k,i,iens,var_name,tvar(k,j,i,iens));
    }
  });
}
void neg_chk(pam::PamCoupler &coupler, std::string id, std::string var_name) {
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nx   = coupler.get_option<int>("crm_nx");
  auto ny   = coupler.get_option<int>("crm_ny");
  auto nz   = coupler.get_option<int>("crm_nz");
  auto nens = coupler.get_option<int>("ncrms");
  auto tvar = dm_device.get<real,4>(var_name);
  parallel_for("Horz mean of CRM dry density", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    if ( tvar(k,j,i,iens)<0 ) {
      printf("  WHDEBUG - negative value detected in driver - k:%d  i:%d  e:%d  %s:%d \n",k,i,iens,var_name,tvar(k,j,i,iens));
    }
  });
}
void max_chk(pam::PamCoupler &coupler, std::string id, std::string var_name, real max_val) {
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nx   = coupler.get_option<int>("crm_nx");
  auto ny   = coupler.get_option<int>("crm_ny");
  auto nz   = coupler.get_option<int>("crm_nz");
  auto nens = coupler.get_option<int>("ncrms");
  auto tvar = dm_device.get<real,4>(var_name);
  parallel_for("Horz mean of CRM dry density", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    if ( tvar(k,j,i,iens)>max_val ) {
      printf("  WHDEBUG - variable exceeds max threshold in driver - k:%d  i:%d  e:%d  %s:%d \n",k,i,iens,var_name,tvar(k,j,i,iens));
    }
  });
}
void chk_state( pam::PamCoupler &coupler, std::string id ) {
  // nan_chk(coupler, id, "temp");
  // nan_chk(coupler, id, "density_dry");
  // nan_chk(coupler, id, "water_vapor");
  // nan_chk(coupler, id, "cloud_water");
  // nan_chk(coupler, id, "ice");
  // nan_chk(coupler, id, "uvel");
  // nan_chk(coupler, id, "wvel");

  // neg_chk(coupler, id, "temp");
  // neg_chk(coupler, id, "density_dry");

  // max_chk(coupler, id, "temp", 400);
  // max_chk(coupler, id, "density_dry", 100);

  // fflush(stdout);
}


void print_state( pam::PamCoupler &coupler, std::string id ) {
  // auto &dm_device = coupler.get_data_manager_device_readwrite();
  // auto nz         = coupler.get_option<int>("crm_nz");
  // auto nx         = coupler.get_option<int>("crm_nx");
  // auto ny         = coupler.get_option<int>("crm_ny");
  // auto nens       = coupler.get_option<int>("ncrms");
  // // auto gcm_nlev   = coupler.get_option<int>("gcm_nlev");
  // auto pmid       = coupler.compute_pressure_array();
  // auto zmid       = dm_device.get<real,2>("vertical_midpoint_height" );
  // // auto gcm_rho_d  = dm_device.get<real,2>("gcm_density_dry");
  // auto temp       = dm_device.get<real,4>("temp");
  // auto crm_rho_v  = dm_device.get<real,4>("water_vapor");
  // auto crm_rho_c  = dm_device.get<real,4>("cloud_water");
  // auto crm_rho_d  = dm_device.get<real,4>("density_dry");
  // for (int k=0; k<nz; k++) { 
  //   // int k_gcm = (gcm_nlev+1)-1-k;
  //   real max_temp = -1e20;
  //   real max_rhod = -1e20;
  //   real max_rhov = -1e20;
  //   real max_rhoc = -1e20;
  //   real min_temp =  1e20;
  //   real min_rhod =  1e20;
  //   real min_rhov =  1e20;
  //   real min_rhoc =  1e20;
  //   for (int j=0; j<ny; j++) { 
  //     for (int i=0; i<nx; i++) { 
  //       for (int n=0; n<nens; n++) { 
  //         max_temp = std::max(max_temp,     temp(k,j,i,n));
  //         max_rhod = std::max(max_rhod,crm_rho_d(k,j,i,n));
  //         max_rhov = std::max(max_rhov,crm_rho_v(k,j,i,n));
  //         max_rhoc = std::max(max_rhoc,crm_rho_c(k,j,i,n));
  //         min_temp = std::min(min_temp,     temp(k,j,i,n));
  //         min_rhod = std::min(min_rhod,crm_rho_d(k,j,i,n));
  //         min_rhov = std::min(min_rhov,crm_rho_v(k,j,i,n));
  //         min_rhoc = std::min(min_rhoc,crm_rho_c(k,j,i,n));
  //       }
  //     }
  //   }
  //   std::cout<<"  WHDEBUG print_state - "<<id
  //   <<"  k:"<<k
  //   <<"  z:"<<zmid(k,0)
  //   <<"  p:"<<pmid(k,0,0,0)
  //   // <<"  t:"<<temp(k,0,0,0)
  //   // <<"  rhod:"<<crm_rho_d(k,0,0,0)
  //   // <<"  rhov:"<<crm_rho_v(k,0,0,0)
  //   // <<"  rhoc:"<<crm_rho_c(k,0,0,0)
  //   <<"  max_temp:"<<max_temp
  //   <<"  max_rhod:"<<max_rhod
  //   <<"  max_rhov:"<<max_rhov
  //   <<"  max_rhoc:"<<max_rhoc
  //   <<"  min_temp:"<<min_temp
  //   <<"  min_rhod:"<<min_rhod
  //   <<"  min_rhov:"<<min_rhov
  //   <<"  min_rhoc:"<<min_rhoc
  //   <<std::endl;
  // }
}

extern "C" void pam_driver() {
  //------------------------------------------------------------------------------------------------
  using yakl::intrinsics::abs;
  using yakl::intrinsics::maxval;
  // using yakl::atomicAdd; // temporary - only for debugging
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
  //------------------------------------------------------------------------------------------------
  // set various coupler options
  coupler.set_option<real>("gcm_physics_dt",gcm_dt);
  coupler.set_option<std::string>("p3_lookup_data_path","./p3_data");
  coupler.set_option<int>("sponge_num_layers",crm_nz*0.25);
  coupler.set_option<real>("sponge_time_scale",60*6);
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
  dycore.init(coupler);
  rad   .init(coupler);
  //------------------------------------------------------------------------------------------------
  // update coupler GCM state with input GCM state
  pam_state_update_gcm_state(coupler);

  // Copy input CRM state (saved by the GCM) to coupler
  pam_state_copy_input_to_coupler(coupler);

  chk_state(coupler, "0");
  print_state(coupler, "0");

  // #if defined(MMF_PAM_DYCOR_SPAM)
  // pam_state_update_dry_density(coupler); // this update effectively disables dry density forcing
  // #endif

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

  // // initilize quantities for surface "psuedo-friction"
  // auto input_tau  = dm_host.get<real const,1>("input_tau00").createDeviceCopy();
  // auto input_bflx = dm_host.get<real const,1>("input_bflxls").createDeviceCopy();
  // modules::surface_friction_init(coupler, input_tau, input_bflx);

  // Perturb the CRM at the beginning of the run
  if (is_first_step) {
    auto gcolp = dm_host.get<int const,1>("gcolp").createDeviceCopy();
    modules::perturb_temperature( coupler , gcolp );
  }

  #if defined(MMF_PAM_DYCOR_SPAM)
    pam_state_update_reference_state(coupler, dycore);
    dycore.pre_time_loop(coupler);
  #endif

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------

  // Run the CRM
  real etime_crm = 0;
  while (etime_crm < gcm_dt) {
    if (crm_dt == 0.) { crm_dt = dycore.compute_time_step(coupler); }
    if (etime_crm + crm_dt > gcm_dt) { crm_dt = gcm_dt - etime_crm; }

    // printf("  WHDEBUG - etime_crm:%d \n",etime_crm);
    // fflush(stdout);
    std::cout<<"  WHDEBUG - etime_crm:"<<etime_crm<<std::endl;

    chk_state(coupler, "1");
    print_state(coupler, "1");

    // run a PAM time step
    coupler.run_module( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.run_module( "radiation"                    , [&] (pam::PamCoupler &coupler) {rad   .timeStep(coupler);} );

    chk_state(coupler, "2");
    print_state(coupler, "2");

    #if defined(MMF_PAM_DYCOR_SPAM)
    // pam_state_update_dry_density(coupler); // redo this update to keep domain mean dry density matching GCM
    pam_state_update_reference_state(coupler, dycore);
    #endif

    // chk_state(coupler, "2.5");
    // print_state(coupler, "2.5");

    pam_statistics_save_state(coupler);
    coupler.run_module( "dycore"                       , [&] (pam::PamCoupler &coupler) {dycore.timeStep(coupler);} );
    pam_statistics_aggregate_tendency(coupler,"dycor");

    chk_state(coupler, "3");
    print_state(coupler, "3");

    pam_statistics_save_state(coupler);
    coupler.run_module( "sponge_layer"                 , modules::sponge_layer );
    pam_statistics_aggregate_tendency(coupler,"sponge");

    chk_state(coupler, "4");

    // coupler.run_module( "compute_surface_friction"     , modules::compute_surface_friction );

    pam_statistics_save_state(coupler);
    coupler.run_module( "sgs"                          , [&] (pam::PamCoupler &coupler) {sgs   .timeStep(coupler);} );
    pam_statistics_aggregate_tendency(coupler,"sgs");

    chk_state(coupler, "5");

    pam_statistics_save_state(coupler);
    coupler.run_module( "micro"                        , [&] (pam::PamCoupler &coupler) {micro .timeStep(coupler);} );
    pam_statistics_aggregate_tendency(coupler,"micro");

    chk_state(coupler, "6");

    pam_radiation_timestep_aggregation(coupler);
    pam_statistics_timestep_aggregation(coupler);

    etime_crm += crm_dt;
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
