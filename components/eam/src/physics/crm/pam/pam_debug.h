#pragma once

#include "pam_coupler.h"

#if defined(__SYCL_DEVICE_ONLY__)
#define PRINTF(format, ...)                                       \
  do {                                                            \
    const __attribute__((opencl_constant)) char fmt[] = (format); \
    sycl::ext::oneapi::experimental::printf(fmt, ##__VA_ARGS__);  \
  } while (0)
#else
  #define PRINTF(format, ...)                                     \
    printf(format, ##__VA_ARGS__);
#endif

// These routines were helpful for debugging the coupling 
// between PAM and E3SM, so we kept them here for future use

// register and initialize various quantities for statistical calculations
inline void pam_debug_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nx   = coupler.get_option<int>("crm_nx");
  auto ny   = coupler.get_option<int>("crm_ny");
  auto nz   = coupler.get_option<int>("crm_nz");
  auto nens = coupler.get_option<int>("ncrms");
  //------------------------------------------------------------------------------------------------
  dm_device.register_and_allocate<real>("debug_save_temp", "saved temp for debug", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("debug_save_rhod", "saved rhod for debug", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("debug_save_rhov", "saved rhov for debug", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("debug_save_rhoc", "saved rhoc for debug", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("debug_save_rhoi", "saved rhoi for debug", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("debug_save_uvel", "saved uvel for debug", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("debug_save_wvel", "saved wvel for debug", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  auto debug_save_temp = dm_device.get<real,4>("debug_save_temp");
  auto debug_save_rhod = dm_device.get<real,4>("debug_save_rhod");
  auto debug_save_rhov = dm_device.get<real,4>("debug_save_rhov");
  auto debug_save_rhoc = dm_device.get<real,4>("debug_save_rhoc");
  auto debug_save_rhoi = dm_device.get<real,4>("debug_save_rhoi");
  auto debug_save_uvel = dm_device.get<real,4>("debug_save_uvel");
  auto debug_save_wvel = dm_device.get<real,4>("debug_save_wvel");
  //------------------------------------------------------------------------------------------------
  auto temp = dm_device.get<real const,4>("temp");
  auto rhod = dm_device.get<real const,4>("density_dry");
  auto rhov = dm_device.get<real const,4>("water_vapor");
  auto rhoc = dm_device.get<real const,4>("cloud_water");
  auto rhoi = dm_device.get<real const,4>("ice");
  auto uvel = dm_device.get<real const,4>("uvel");
  auto wvel = dm_device.get<real const,4>("wvel");
  //------------------------------------------------------------------------------------------------
  parallel_for("copy data to saved debug variables", SimpleBounds<4>(nz,ny,nx,nens), 
    YAKL_LAMBDA (int k, int j, int i, int iens) {
    debug_save_temp(k,j,i,iens) = temp(k,j,i,iens);
    debug_save_rhod(k,j,i,iens) = rhod(k,j,i,iens);
    debug_save_rhov(k,j,i,iens) = rhov(k,j,i,iens);
    debug_save_rhoc(k,j,i,iens) = rhoc(k,j,i,iens);
    debug_save_rhoi(k,j,i,iens) = rhoi(k,j,i,iens);
    debug_save_uvel(k,j,i,iens) = uvel(k,j,i,iens);
    debug_save_wvel(k,j,i,iens) = wvel(k,j,i,iens);
  });
  //------------------------------------------------------------------------------------------------
}

// main routine to encapsulate checking for bad values across many variables
void pam_debug_check_state( pam::PamCoupler &coupler, int id, int nstep ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readonly();
  auto nx   = coupler.get_option<int>("crm_nx");
  auto ny   = coupler.get_option<int>("crm_ny");
  auto nz   = coupler.get_option<int>("crm_nz");
  auto nens = coupler.get_option<int>("ncrms");
  auto temp = dm_device.get<real const,4>("temp");
  auto rhod = dm_device.get<real const,4>("density_dry");
  auto rhov = dm_device.get<real const,4>("water_vapor");
  auto rhoc = dm_device.get<real const,4>("cloud_water");
  auto rhoi = dm_device.get<real const,4>("ice");
  auto uvel = dm_device.get<real const,4>("uvel");
  auto wvel = dm_device.get<real const,4>("wvel");
  auto debug_save_temp = dm_device.get<real,4>("debug_save_temp");
  auto debug_save_rhod = dm_device.get<real,4>("debug_save_rhod");
  auto debug_save_rhov = dm_device.get<real,4>("debug_save_rhov");
  auto debug_save_rhoc = dm_device.get<real,4>("debug_save_rhoc");
  auto debug_save_rhoi = dm_device.get<real,4>("debug_save_rhoi");
  auto debug_save_uvel = dm_device.get<real,4>("debug_save_uvel");
  auto debug_save_wvel = dm_device.get<real,4>("debug_save_wvel");
  auto lat        = dm_host.get<real const,1>("latitude"  ).createDeviceCopy();
  auto lon        = dm_host.get<real const,1>("longitude" ).createDeviceCopy();
  auto input_phis = dm_host.get<real const,1>("input_phis").createDeviceCopy();
  real grav = 9.80616;
  //------------------------------------------------------------------------------------------------
  // Check for invalid values
  parallel_for("", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    auto phis = input_phis(iens)/grav;
    // Check for NaNs
    const auto is_nan_t_atm = isnan( temp(k,j,i,iens) );
    const auto is_nan_d_atm = isnan( rhod(k,j,i,iens) );
    const auto is_nan_q_atm = isnan( rhov(k,j,i,iens) );
    if ( is_nan_t_atm || is_nan_q_atm || is_nan_d_atm ) {
      auto phis = input_phis(iens)/grav;
      printf("PAM-DEBUG nan-found - st:%3.3d id:%2.2d k:%3.3d i:%3.3d n:%3.3d y:%5.1f x:%5.1f ph:%6.1f -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g \n",
        nstep,id,k,i,iens,lat(iens),lon(iens),phis,
        temp(k,j,i,iens),
        rhod(k,j,i,iens),
        rhov(k,j,i,iens),
        rhoc(k,j,i,iens),
        rhoi(k,j,i,iens),
        uvel(k,j,i,iens),
        wvel(k,j,i,iens),
        debug_save_temp(k,j,i,iens),
        debug_save_rhod(k,j,i,iens),
        debug_save_rhov(k,j,i,iens),
        debug_save_rhoc(k,j,i,iens),
        debug_save_rhoi(k,j,i,iens),
        debug_save_uvel(k,j,i,iens),
        debug_save_wvel(k,j,i,iens)
      );
    }
    // Check for negative values
    const auto is_neg_t_atm = temp(k,j,i,iens)<0;
    const auto is_neg_d_atm = rhod(k,j,i,iens)<0;
    const auto is_neg_q_atm = rhov(k,j,i,iens)<0;
    if ( is_neg_t_atm || is_neg_q_atm || is_neg_d_atm ) {
      auto phis = input_phis(iens)/grav;
      printf("PAM-DEBUG neg-found - st:%3.3d id:%2.2d k:%3.3d i:%3.3d n:%3.3d y:%5.1f x:%5.1f ph:%6.1f -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g \n",
        nstep,id,k,i,iens,lat(iens),lon(iens),phis,
        temp(k,j,i,iens),
        rhod(k,j,i,iens),
        rhov(k,j,i,iens),
        rhoc(k,j,i,iens),
        rhoi(k,j,i,iens),
        uvel(k,j,i,iens),
        wvel(k,j,i,iens),
        debug_save_temp(k,j,i,iens),
        debug_save_rhod(k,j,i,iens),
        debug_save_rhov(k,j,i,iens),
        debug_save_rhoc(k,j,i,iens),
        debug_save_rhoi(k,j,i,iens),
        debug_save_uvel(k,j,i,iens),
        debug_save_wvel(k,j,i,iens)
      );
    }
    // Check for low temperature
    const auto is_low_t = temp(k,j,i,iens)<100;
    if ( is_low_t ) {
      printf("PAM-DEBUG low-T - st:%3.3d id:%2.2d k:%3.3d i:%3.3d n:%3.3d y:%5.1f x:%5.1f ph:%6.1f -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g \n",
        nstep,id,k,i,iens,lat(iens),lon(iens),phis,
        temp(k,j,i,iens),
        rhod(k,j,i,iens),
        rhov(k,j,i,iens),
        rhoc(k,j,i,iens),
        rhoi(k,j,i,iens),
        uvel(k,j,i,iens),
        wvel(k,j,i,iens),
        debug_save_temp(k,j,i,iens),
        debug_save_rhod(k,j,i,iens),
        debug_save_rhov(k,j,i,iens),
        debug_save_rhoc(k,j,i,iens),
        debug_save_rhoi(k,j,i,iens),
        debug_save_uvel(k,j,i,iens),
        debug_save_wvel(k,j,i,iens)
      );
    }
    // Check for large vertical velocity
    const auto is_large_pos_w = wvel(k,j,i,iens)> 40;
    const auto is_large_neg_w = wvel(k,j,i,iens)<-40;
    if ( is_large_pos_w || is_large_neg_w ) {
      printf("PAM-DEBUG large-W - st:%3.3d id:%2.2d k:%3.3d i:%3.3d n:%3.3d y:%5.1f x:%5.1f ph:%6.1f -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g -- t:%8.2g rd:%8.2g rv:%8.2g rc:%8.2g ri:%8.2g u:%8.2g w:%8.2g \n",
        nstep,id,k,i,iens,lat(iens),lon(iens),phis,
        temp(k,j,i,iens),
        rhod(k,j,i,iens),
        rhov(k,j,i,iens),
        rhoc(k,j,i,iens),
        rhoi(k,j,i,iens),
        uvel(k,j,i,iens),
        wvel(k,j,i,iens),
        debug_save_temp(k,j,i,iens),
        debug_save_rhod(k,j,i,iens),
        debug_save_rhov(k,j,i,iens),
        debug_save_rhoc(k,j,i,iens),
        debug_save_rhoi(k,j,i,iens),
        debug_save_uvel(k,j,i,iens),
        debug_save_wvel(k,j,i,iens)
      );
    }
    // update saved previous values
    debug_save_temp(k,j,i,iens) = temp(k,j,i,iens);
    debug_save_rhod(k,j,i,iens) = rhod(k,j,i,iens);
    debug_save_rhov(k,j,i,iens) = rhov(k,j,i,iens);
    debug_save_rhoc(k,j,i,iens) = rhoc(k,j,i,iens);
    debug_save_rhoi(k,j,i,iens) = rhoi(k,j,i,iens);
    debug_save_uvel(k,j,i,iens) = uvel(k,j,i,iens);
    debug_save_wvel(k,j,i,iens) = wvel(k,j,i,iens);
  });
  //------------------------------------------------------------------------------------------------
}

// print the min and max of PAM state variables to help look for problems
// void print_state_min_max( pam::PamCoupler &coupler, std::string id ) {
  // auto &dm_device = coupler.get_data_manager_device_readonly();
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
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readonly();
  auto nz     = coupler.get_option<int>("crm_nz");
  auto nx     = coupler.get_option<int>("crm_nx");
  auto ny     = coupler.get_option<int>("crm_ny");
  auto nens   = coupler.get_option<int>("ncrms");
  // auto pmid   = coupler.compute_pressure_array();
  // auto zmid   = dm_device.get<real,2>("vertical_midpoint_height" );
  auto temp   = dm_device.get<real const,4>("temp");
  auto rho_v  = dm_device.get<real const,4>("water_vapor");
  auto rho_c  = dm_device.get<real const,4>("cloud_water");
  auto rho_d  = dm_device.get<real const,4>("density_dry");
  auto rho_i  = dm_device.get<real const,4>("ice");
  // parallel_for("", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
  parallel_for("pam_debug_print_state", SimpleBounds<2>(nz,nx), YAKL_LAMBDA (int k, int i) {
    int j = 0;
    int n = 0;
    PRINTF("PAM-DEBUG %d - k:%d  i:%d  temp : %g  rv: %g  rc: %g  ri: %g \n",
      id,k,i,
      temp(k,j,i,n),
      rho_v(k,j,i,n),
      rho_c(k,j,i,n),
      rho_i(k,j,i,n)
    );
  });
}
