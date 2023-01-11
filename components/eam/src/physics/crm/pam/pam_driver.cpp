#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "radiation.h"
#include "pam_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
// #include "gcm_rad_coupling.h"
#include "pam_feedback.h"
#include "pam_state_update.h"
#include "pam_radiation.h"
#include "pam_statistics.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"
#include "broadcast_initial_gcm_column.h"
// #include "output.h"

extern "C" void pam_driver() {
  //------------------------------------------------------------------------------------------------
  using yakl::intrinsics::abs;
  using yakl::intrinsics::maxval;
  auto &coupler = pam_interface::get_coupler();
  //------------------------------------------------------------------------------------------------
  // retreive coupler options
  auto nens        = coupler.get_option<int>("ncrms");
  auto gcm_nlev    = coupler.get_option<int>("gcm_nlev");
  auto crm_nz      = coupler.get_option<int>("crm_nz");
  auto crm_nx      = coupler.get_option<int>("crm_nx");
  auto crm_ny      = coupler.get_option<int>("crm_ny");
  auto crm_dx      = coupler.get_option<double>("crm_dx");
  auto crm_dy      = coupler.get_option<double>("crm_dy");
  auto gcm_dt      = coupler.get_option<double>("gcm_dt");
  auto crm_dt      = coupler.get_option<double>("crm_dt");
  coupler.set_option<real>("gcm_physics_dt",gcm_dt);
  // auto rad_ny      = coupler.get_option<int>("rad_ny");
  // auto rad_nx      = coupler.get_option<int>("rad_nx");
  //------------------------------------------------------------------------------------------------
  // Allocate the coupler state and retrieve host/device data managers
  coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  // wrap host data in YAKL arrays
  auto input_zint          = dm_host.get<real const,2>("input_zint").createDeviceCopy();
  // auto input_zmid          = dm_host.get<real const,2>("input_zmid").createDeviceCopy();
  // auto input_pint          = dm_host.get<real const,2>("input_pint").createDeviceCopy();
  // auto input_pmid          = dm_host.get<real const,2>("input_pmid").createDeviceCopy();
  // auto input_pdel          = dm_host.get<real const,2>("input_pdel").createDeviceCopy();
  // auto input_tau00         = dm_host.get<real const,2>("input_tau00").createDeviceCopy();
  // auto input_bflxls        = dm_host.get<real const,2>("input_bflxls").createDeviceCopy();
  // auto input_wndls         = dm_host.get<real const,2>("input_wndls").createDeviceCopy();
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

  // Set the vertical grid in the coupler
  coupler.set_grid( crm_dx , crm_dy , input_zint );

  // update coupler GCM state with input GCM state
  update_gcm_state(coupler);

  // Copy input CRM state saved by the GCM to coupler
  copy_input_state_to_coupler(coupler);

  // Copy input rad tendencies to coupler
  pam_radiation_copy_input_to_coupler(coupler);

  // initialize rad output variables
  pam_radiation_init(coupler);

  // initialize stat variables
  pam_statistics_init(coupler);

  // set CRM dry density using gcm_density_dry (set in update_gcm_state)
  modules::broadcast_initial_gcm_column_dry_density(coupler); 

  // Compute CRM forcing tendencies
  modules::compute_gcm_forcing_tendencies(coupler);

  // Define hydrostasis (only for PAM-A/AWFL)
  coupler.update_hydrostasis();

  // Perturb the CRM at the beginning of the run
  auto is_first_step = coupler.get_option<bool>("is_first_step");
  if (is_first_step) {
    auto gcolp = dm_host.get<int const,1>("gcolp").createDeviceCopy();
    modules::perturb_temperature( coupler , gcolp );
  }

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  // Run the CRM
  real etime_crm = 0;
  while (etime_crm < gcm_dt) {
    if (crm_dt == 0.) { crm_dt = dycore.compute_time_step(coupler); }
    if (etime_crm + crm_dt > gcm_dt) { crm_dt = gcm_dt - etime_crm; }

    // run a PAM time step
    coupler.run_module( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.run_module( "sponge_layer"                 , modules::sponge_layer );
    coupler.run_module( "radiation"                    , [&] (pam::PamCoupler &coupler) {rad   .timeStep(coupler);} );
    coupler.run_module( "dycore"                       , [&] (pam::PamCoupler &coupler) {dycore.timeStep(coupler);} );
    coupler.run_module( "sgs"                          , [&] (pam::PamCoupler &coupler) {sgs   .timeStep(coupler);} );
    coupler.run_module( "micro"                        , [&] (pam::PamCoupler &coupler) {micro .timeStep(coupler);} );

    pam_radiation_timestep_aggregation(coupler);
    pam_statistics_timestep_aggregation(coupler);

    etime_crm += crm_dt;
  }
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------

  // Compute primary feedback tendencies and copy to GCM
  pam_feedback_compute_crm_feedback_tendencies( coupler, gcm_dt );
  
  // Compute horizontal means of CRM state variables that are not forced
  pam_feedback_compute_crm_mean_state(coupler);

  // Copy the output CRM state from the PAM coupler to the GCM
  pam_feedback_copy_output_state_to_gcm(coupler);

  // Copy radiation column quantities to host
  pam_radiation_copy_output_to_gcm(coupler);

  // copy aggregated statistical quantities to host
  pam_statistics_copy_to_host( coupler, gcm_dt );

  //------------------------------------------------------------------------------------------------
  // Finalize and clean up
  micro .finalize(coupler);
  sgs   .finalize(coupler);
  dycore.finalize(coupler);
  rad   .finalize(coupler);
  pam_interface::finalize();
  //------------------------------------------------------------------------------------------------
}
