#pragma once

#include "pam_coupler.h"

// These routines were helpful for debugging the coupling 
// between PAM and E3SM, so we kept them here for future use


// main routine to encapsulate checking for bad values across many variables
void pam_debug_check_state( pam::PamCoupler &coupler, int id, int nstep ) {

  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nx   = coupler.get_option<int>("crm_nx");
  auto ny   = coupler.get_option<int>("crm_ny");
  auto nz   = coupler.get_option<int>("crm_nz");
  auto nens = coupler.get_option<int>("ncrms");
  auto temp  = dm_device.get<real,4>("temp");
  auto rho_d = dm_device.get<real,4>("density_dry");
  auto rho_v = dm_device.get<real,4>("water_vapor");
  auto rho_c = dm_device.get<real,4>("cloud_water");
  auto rho_i = dm_device.get<real,4>("ice");

  // Check for NaNs
  parallel_for("", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {

    const auto is_nan_t_atm = isnan( temp(k,j,i,iens) );
    const auto is_nan_r_atm = isnan( rho_d(k,j,i,iens) );
    const auto is_nan_q_atm = isnan( rho_v(k,j,i,iens) );
    if ( is_nan_t_atm || is_nan_r_atm || is_nan_q_atm ) {
      printf("PAM-DEBUG nan-found - st:%d  id:%d  k:%d  i:%d  n:%d  t: %g  rv: %g  rd: %g \n",
        nstep,id,k,i,nens,
        temp (k,j,i,iens),
        rho_v(k,j,i,iens),
        rho_d(k,j,i,iens)
      );
    }

    const auto is_neg_t_atm = temp(k,j,i,iens)  < 0;
    const auto is_neg_r_atm = rho_d(k,j,i,iens) < 0;
    const auto is_neg_q_atm = rho_v(k,j,i,iens) < 0;
    if ( is_neg_t_atm || is_neg_r_atm || is_neg_q_atm ) {
      printf("PAM-DEBUG negative-found - st:%d  id:%d  k:%d  i:%d  n:%d  t: %g  rv: %g  rd: %g \n",
        nstep,id,k,i,nens,
        temp (k,j,i,iens),
        rho_v(k,j,i,iens),
        rho_d(k,j,i,iens)
      );
    }

  });

}

// print the min and max of PAM state variables to help look for problems
// void print_state_min_max( pam::PamCoupler &coupler, std::string id ) {
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
  //   std::cout<<"PAM-DEBUG print_state_min_max - "<<id
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
// }


// print certain state variables
void pam_debug_print_state( pam::PamCoupler &coupler, int id ) {
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nz     = coupler.get_option<int>("crm_nz");
  auto nx     = coupler.get_option<int>("crm_nx");
  auto ny     = coupler.get_option<int>("crm_ny");
  auto nens   = coupler.get_option<int>("ncrms");
  // auto pmid   = coupler.compute_pressure_array();
  // auto zmid   = dm_device.get<real,2>("vertical_midpoint_height" );
  auto temp   = dm_device.get<real,4>("temp");
  auto rho_v  = dm_device.get<real,4>("water_vapor");
  auto rho_c  = dm_device.get<real,4>("cloud_water");
  auto rho_d  = dm_device.get<real,4>("density_dry");
  auto rho_i  = dm_device.get<real,4>("ice");
  // parallel_for("", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
  parallel_for("pam_debug_print_state", SimpleBounds<2>(nz,nx), YAKL_LAMBDA (int k, int i) {
    int j = 0;
    int n = 0;
    printf("PAM-DEBUG %d - k:%d  i:%d  temp : %g  rv: %g  rc: %g  ri: %g \n",
      id,k,i,
      temp(k,j,i,n),
      rho_v(k,j,i,n),
      rho_c(k,j,i,n),
      rho_i(k,j,i,n)
    );
  });
}