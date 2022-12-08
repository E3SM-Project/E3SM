#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "pam_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "pam_feedback.h"
#include "pam_state_update.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"
#include "broadcast_initial_gcm_column.h"
#include "output.h"

extern "C" void pam_driver() {
  //------------------------------------------------------------------------------------------------
  using yakl::intrinsics::abs;
  using yakl::intrinsics::maxval;

  //------------------------------------------------------------------------------------------------
  // retreive coupler options
  auto nens        = coupler.get_option<int>('ncrms');
  auto gcm_nlev    = coupler.get_option<int>('gcm_nlev');
  auto crm_nz      = coupler.get_option<int>('crm_nz');
  auto crm_nx      = coupler.get_option<int>('crm_nx');
  auto crm_ny      = coupler.get_option<int>('crm_ny');
  auto crm_dx      = coupler.get_option<double>('crm_dx');
  auto crm_dy      = coupler.get_option<double>('crm_dy');
  // time step increments
  auto gcm_dt      = coupler.get_option<double>('gcm_dt');
  auto crm_dt      = coupler.get_option<double>('crm_dt');
  // flag to determine whether to perturb CRM state
  auto is_first_step = coupler.get_option<double>('is_first_step');
  //------------------------------------------------------------------------------------------------
  // Allocate the coupler state and retrieve host/device data managers
  auto &coupler = pam_interface::get_coupler();
  coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  // wrap host data in YAKL arrays
  auto gcolp = dm_host.get<real const,1>("gcolp").createDeviceCopy();

  // auto input_bflxls        = dm_host.get<real const,2>("input_bflxls").createDeviceCopy();
  // auto input_wndls         = dm_host.get<real const,2>("input_wndls").createDeviceCopy();
  // auto input_zmid          = dm_host.get<real const,2>("input_zmid").createDeviceCopy();
  // auto input_zint          = dm_host.get<real const,2>("input_zint").createDeviceCopy();
  // auto input_pmid          = dm_host.get<real const,2>("input_pmid").createDeviceCopy();
  // auto input_pint          = dm_host.get<real const,2>("input_pint").createDeviceCopy();
  // auto input_pdel          = dm_host.get<real const,2>("input_pdel").createDeviceCopy();
  // auto input_tau00         = dm_host.get<real const,2>("input_tau00").createDeviceCopy();

  //------------------------------------------------------------------------------------------------
  // Create objects for dycor, microphysics, and turbulence
  Dycore       dycore;
  Microphysics micro;
  SGS          sgs;
  //------------------------------------------------------------------------------------------------
  // Initialize dycor, microphysics, and turbulence
  dycore.init( coupler );
  micro .init( coupler );
  sgs   .init( coupler );
  //------------------------------------------------------------------------------------------------
  // Set physical constants for coupler at thread 0 using microphysics data
  coupler.set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );
  //------------------------------------------------------------------------------------------------
  // Set the vertical grid in the coupler
  coupler.set_grid( crm_dx , crm_dy , input_zint );
  //------------------------------------------------------------------------------------------------
  // Define hydrostasis (only for PAM-A / AWFL dycor)
  coupler.update_hydrostasis();
  //------------------------------------------------------------------------------------------------
  // update the coupler GCM state variables using the input GCM state
  update_gcm_state( coupler )
  //------------------------------------------------------------------------------------------------
  // Copy the input CRM state saved by the GCM to the PAM coupler state
  copy_input_state_to_coupler( coupler );
  // Copy input radiation tendencies to coupler
  copy_input_rad_to_coupler( coupler );
  //------------------------------------------------------------------------------------------------
  // Compute CRM forcing tendencies
  modules::compute_gcm_forcing_tendencies( coupler, gcm_dt );
  //------------------------------------------------------------------------------------------------
  // Perturb the CRM only at the beginning of the run
  if (is_first_step) { modules::perturb_temperature( coupler , gcolp ); }
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  // Run the CRM
  real etime_crm = 0;
  while (etime_crm < gcm_dt) {
    if (crm_dt == 0.) { crm_dt = dycore.compute_time_step(coupler); }
    if (etime_crm + crm_dt > simTime_crm) { crm_dt = simTime_crm - etime_crm; }

    coupler.run_module( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.run_module( "apply_rad_forcing_tendencies" , modules::apply_rad_forcing_tendencies );
    coupler.run_module( "sponge_layer"                 , modules::sponge_layer );
    coupler.run_module( "dycore"                       , [&] (pam::PamCoupler &coupler) {dycore.timeStep(coupler);} );
    coupler.run_module( "sgs"                          , [&] (pam::PamCoupler &coupler) {sgs   .timeStep(coupler);} );
    coupler.run_module( "micro"                        , [&] (pam::PamCoupler &coupler) {micro .timeStep(coupler);} );

    etime_crm += crm_dt;
  }
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  // Compute primary feedback tendencies for GCM
  compute_crm_feedback_tendencies( coupler, gcm_dt )
  
  // Compute horizontal means of CRM state variables that are not forced
  compute_crm_mean_state( coupler )

  // Copy the output CRM state from the PAM coupler to the GCM
  copy_output_state_to_gcm( coupler );

  // Copy input radiation tendencies to coupler
  copy_output_rad_to_gcm( coupler );

  //------------------------------------------------------------------------------------------------
  // Output aggregated surface flux of water specices
  // precip_liq = dm_device.get<real const,2>("p3_output_prec_liq")
  // precip_ice = dm_device.get<real const,2>("p3_output_prec_ice")
  // precip_tot = precip_liq + precip_ice
  // auto precip_tot_host = dm_host.get<real,1>("output_precc");
  // auto precip_liq_host = dm_host.get<real,1>("output_precl");
  // auto precip_ice_host = dm_host.get<real,1>("output_precsc");
  // precip_tot.deep_copy_to(precip_tot_host);
  // precip_liq.deep_copy_to(precip_liq_host);
  // precip_ice.deep_copy_to(precip_ice_host);
  //------------------------------------------------------------------------------------------------
  // Finalize the coupler and clean up
  dycore.finalize( coupler );
  pam_interface::finalize();
}