#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "pam_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"
#include "broadcast_initial_gcm_column.h"
#include "output.h"

extern "C" void pam_driver() {
  //------------------------------------------------------------------------------------------------
  using yakl::intrinsics::abs;
  using yakl::intrinsics::maxval;
  //------------------------------------------------------------------------------------------------
  auto &coupler = pam_interface::get_coupler();
  //------------------------------------------------------------------------------------------------
  // retreive coupler options
  auto ncrms       = coupler.get_option<int>('ncrms');
  auto gcm_nlev    = coupler.get_option<int>('gcm_nlev');
  auto crm_nz      = coupler.get_option<int>('crm_nz');
  auto crm_nx      = coupler.get_option<int>('crm_nx');
  auto crm_ny      = coupler.get_option<int>('crm_ny');
  auto crm_dx      = coupler.get_option<double>('crm_dx');
  auto crm_dy      = coupler.get_option<double>('crm_dy');
  auto gcm_dt      = coupler.get_option<double>('gcm_dt');
  auto crm_dt      = coupler.get_option<double>('crm_dt');
  //------------------------------------------------------------------------------------------------
  // TODO: add check to ensure crm_dt evenly divides gcm_dt
  //------------------------------------------------------------------------------------------------
  // Allocate the coupler state and retrieve host/device data managers
  coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );
  auto &dm_device = coupler.get_data_manager_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  // wrap host data in YAKL array
  auto input_bflxls        = dm_host.get<real const,2>("input_bflxls").createDeviceCopy();
  auto input_wndls         = dm_host.get<real const,2>("input_wndls").createDeviceCopy();
  auto input_zmid          = dm_host.get<real const,2>("input_zmid").createDeviceCopy();
  auto input_zint          = dm_host.get<real const,2>("input_zint").createDeviceCopy();
  auto input_pmid          = dm_host.get<real const,2>("input_pmid").createDeviceCopy();
  auto input_pint          = dm_host.get<real const,2>("input_pint").createDeviceCopy();
  auto input_pdel          = dm_host.get<real const,2>("input_pdel").createDeviceCopy();
  auto input_ul            = dm_host.get<real const,2>("input_ul").createDeviceCopy();
  auto input_vl            = dm_host.get<real const,2>("input_vl").createDeviceCopy();
  auto input_tl            = dm_host.get<real const,2>("input_tl").createDeviceCopy();
  auto input_qccl          = dm_host.get<real const,2>("input_qccl").createDeviceCopy();
  auto input_qiil          = dm_host.get<real const,2>("input_qiil").createDeviceCopy();
  auto input_ql            = dm_host.get<real const,2>("input_ql").createDeviceCopy();
  auto input_tau00         = dm_host.get<real const,2>("input_tau00").createDeviceCopy();
  auto input_ul_esmt       = dm_host.get<real const,2>("input_ul_esmt").createDeviceCopy();
  auto input_vl_esmt       = dm_host.get<real const,2>("input_vl_esmt").createDeviceCopy();
  auto input_t_vt          = dm_host.get<real const,2>("input_t_vt").createDeviceCopy();
  auto input_q_vt          = dm_host.get<real const,2>("input_q_vt").createDeviceCopy();
  auto input_u_vt          = dm_host.get<real const,2>("input_u_vt").createDeviceCopy();

  auto state_u_wind        = dm_host.get<real const,2>("state_u_wind").createDeviceCopy();
  auto state_v_wind        = dm_host.get<real const,2>("state_v_wind").createDeviceCopy();
  auto state_w_wind        = dm_host.get<real const,2>("state_w_wind").createDeviceCopy();
  auto state_temperature   = dm_host.get<real const,2>("state_temperature").createDeviceCopy();
  auto state_qt            = dm_host.get<real const,2>("state_qt").createDeviceCopy();
  auto state_qp            = dm_host.get<real const,2>("state_qp").createDeviceCopy();
  auto state_qn            = dm_host.get<real const,2>("state_qn").createDeviceCopy();

  //------------------------------------------------------------------------------------------------
  // define GCM state used for forcing
  auto gcm_rho_d = dm.get<real,2>("gcm_density_dry");
  auto gcm_rho_v = dm.get<real,2>("gcm_water_vapor");
  auto gcm_uvel  = dm.get<real,2>("gcm_uvel"       );
  auto gcm_vvel  = dm.get<real,2>("gcm_vvel"       );
  auto gcm_wvel  = dm.get<real,2>("gcm_wvel"       );
  auto gcm_temp  = dm.get<real,2>("gcm_temp"       );

  parallel_for( Bounds<2>(crm_nz,nens) , YAKL_LAMBDA (int k, int iens) {
    gcm_rho_d(k,iens) = rho_d_col(k);
    gcm_rho_v(k,iens) = rho_v_col(k);
    gcm_uvel (k,iens) = uvel_col (k);
    gcm_vvel (k,iens) = vvel_col (k);
    gcm_wvel (k,iens) = wvel_col (k);
    gcm_temp (k,iens) = temp_col (k);
  });
  //------------------------------------------------------------------------------------------------
  // NORMALLY THIS WOULD BE DONE INSIDE THE CRM, BUT WE'RE USING CONSTANTS DEFINED BY THE CRM MICRO SCHEME
  // Create the dycore and the microphysics
  Dycore       dycore;
  Microphysics micro;
  SGS          sgs;

  // Set physical constants for coupler at thread 0 using microphysics data
  coupler.set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

  // Set the vertical grid in the coupler
  coupler.set_grid( crm_dx , crm_dy , input_zint );

  micro .init( coupler );
  sgs   .init( coupler );
  dycore.init( coupler );

  // Initialize the CRM internal state from the initial GCM column and random temperature perturbations
  modules::broadcast_initial_gcm_column( coupler );

  // Now that we have an initial state, define hydrostasis for each ensemble member
  coupler.update_hydrostasis();

  modules::perturb_temperature( coupler , 0 );

  coupler.add_pam_function( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
  coupler.add_pam_function( "dycore" , [&] (pam::PamCoupler &coupler, real dt) { dycore.timeStep(coupler,dt); } );
  coupler.add_pam_function( "sgs"    , [&] (pam::PamCoupler &coupler, real dt) { sgs   .timeStep(coupler,dt); } );
  coupler.add_pam_function( "micro"  , [&] (pam::PamCoupler &coupler, real dt) { micro .timeStep(coupler,dt); } );
  coupler.add_pam_function( "sponge_layer"                 , modules::sponge_layer                 );
  // coupler.add_dycore_function( "saturation_adjustment" , saturation_adjustment );

  yakl::timer_start("main_loop");
  while (etime_gcm < simTime) {
    if (etime_gcm + gcm_dt > simTime) { gcm_dt = simTime - etime_gcm; }

    modules::compute_gcm_forcing_tendencies( coupler , gcm_dt );

    real etime_crm = 0;
    real simTime_crm = gcm_dt;
    while (etime_crm < simTime_crm) {
      if (crm_dt == 0.) { crm_dt = dycore.compute_time_step(coupler); }
      if (etime_crm + crm_dt > simTime_crm) { crm_dt = simTime_crm - etime_crm; }

      coupler.run_pam_function( "apply_gcm_forcing_tendencies" , crm_dt );
      coupler.run_pam_function( "dycore"                       , crm_dt );
      coupler.run_pam_function( "sponge_layer"                 , crm_dt );
      coupler.run_pam_function( "sgs"                          , crm_dt );
      coupler.run_pam_function( "micro"                        , crm_dt );

      etime_crm += crm_dt;
      etime_gcm += crm_dt;
      // if (out_freq >= 0. && etime_gcm / out_freq >= num_out+1) {
      //   output( coupler , out_prefix , etime_gcm );
      //   real maxw = maxval(abs(dm.get_collapsed<real const>("wvel")));
      //   if (mainproc) {
      //     std::cout << "Etime , dtphys, maxw: " << etime_gcm << " , " 
      //                                           << crm_dt    << " , "
      //                                           << std::setw(10) << maxw << std::endl;
      //   }
      //   num_out++;
      // }
    }
  }

  if (mainproc) {
    std::cout << "Elapsed Time: " << etime_gcm << "\n";
  }

  dycore.finalize( coupler );

  pam_interface::finalize();

}