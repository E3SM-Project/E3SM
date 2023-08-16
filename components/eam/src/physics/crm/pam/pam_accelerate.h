#pragma once

#include "pam_coupler.h"

void pam_accelerate_nstop( pam::PamCoupler &coupler, int &nstop) {
  auto crm_accel_factor = coupler.get_option<real>("crm_accel_factor");
  if(nstop%static_cast<int>((1+crm_accel_factor)) != 0) {
    printf("pam_accelerate_nstop: Error: (1+crm_accel_factor) does not divide equally into nstop: %4.4d  crm_accel_factor: %6.1f \n",nstop, crm_accel_factor);
    exit(-1);
  } else {
    nstop = nstop / (1 + crm_accel_factor);
  }
}


inline void pam_accelerate_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto nx           = coupler.get_option<int>("crm_nx");
  auto crm_accel_uv = coupler.get_option<bool>("crm_accel_uv");
  //------------------------------------------------------------------------------------------------
  dm_device.register_and_allocate<real>("accel_save_t", "saved temperature for MSA", {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("accel_save_r", "saved dry density for MSA", {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("accel_save_q", "saved total water for MSA", {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("accel_save_u", "saved uvel for MSA",        {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("accel_save_v", "saved vvel for MSA",        {nz,nens}, {"z","nens"} );
  //------------------------------------------------------------------------------------------------
}


inline void pam_accelerate_diagnose( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto ny         = coupler.get_option<int>("crm_ny");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto crm_accel_uv = coupler.get_option<bool>("crm_accel_uv");
  //------------------------------------------------------------------------------------------------
  auto temp = dm_device.get<real,4>("temp"       );
  auto rhod = dm_device.get<real,4>("density_dry");
  auto rhov = dm_device.get<real,4>("water_vapor");
  auto rhol = dm_device.get<real,4>("cloud_water");
  auto rhoi = dm_device.get<real,4>("ice"        );
  auto uvel = dm_device.get<real,4>("uvel"       );
  auto vvel = dm_device.get<real,4>("vvel"       );
  auto accel_save_t = dm_device.get<real,2>("accel_save_t");
  auto accel_save_r = dm_device.get<real,2>("accel_save_r");
  auto accel_save_q = dm_device.get<real,2>("accel_save_q");
  auto accel_save_u = dm_device.get<real,2>("accel_save_u");
  auto accel_save_v = dm_device.get<real,2>("accel_save_v");
  //------------------------------------------------------------------------------------------------
  // compute horizontal means needed later for mean-state acceleration
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int n) {
    accel_save_t(k,n) = 0.0;
    accel_save_r(k,n) = 0.0;
    accel_save_q(k,n) = 0.0;
    if (crm_accel_uv) {
      accel_save_u(k,n) = 0.0;
      accel_save_v(k,n) = 0.0;
    }
  });
  real r_nx_ny  = 1._fp/(nx*ny);  // precompute reciprocal to avoid costly divisions
  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
    yakl::atomicAdd( accel_save_t(k,n), temp(k,j,i,n) * r_nx_ny );
    yakl::atomicAdd( accel_save_r(k,n), rhod(k,j,i,n) * r_nx_ny );
    yakl::atomicAdd( accel_save_q(k,n), ( rhov(k,j,i,n) + rhol(k,j,i,n) + rhoi(k,j,i,n) ) * r_nx_ny );
    if (crm_accel_uv) {
      yakl::atomicAdd( accel_save_u(k,n), uvel(k,j,i,n) * r_nx_ny );
      yakl::atomicAdd( accel_save_v(k,n), vvel(k,j,i,n) * r_nx_ny );
    }
  });
  //------------------------------------------------------------------------------------------------
}


inline void pam_accelerate( pam::PamCoupler &coupler, int nstep, int &nstop ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  using yakl::ScalarLiveOut;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto ny         = coupler.get_option<int>("crm_ny");
  auto nx         = coupler.get_option<int>("crm_nx");
  //------------------------------------------------------------------------------------------------
  auto temp = dm_device.get<real,4>("temp"       );
  auto rhod = dm_device.get<real,4>("density_dry");
  auto rhov = dm_device.get<real,4>("water_vapor");
  auto rhol = dm_device.get<real,4>("cloud_water");
  auto rhoi = dm_device.get<real,4>("ice"        );
  auto uvel = dm_device.get<real,4>("uvel"       );
  auto vvel = dm_device.get<real,4>("vvel"       );
  auto accel_save_t = dm_device.get<real,2>("accel_save_t");
  auto accel_save_r = dm_device.get<real,2>("accel_save_r");
  auto accel_save_q = dm_device.get<real,2>("accel_save_q");
  auto accel_save_u = dm_device.get<real,2>("accel_save_u");
  auto accel_save_v = dm_device.get<real,2>("accel_save_v");
  //------------------------------------------------------------------------------------------------
  bool crm_accel_uv     = coupler.get_option<bool>("crm_accel_uv");
  real crm_accel_factor = coupler.get_option<real>("crm_accel_factor");
  //------------------------------------------------------------------------------------------------
  real2d hmean_t  ("hmean_t",   nz,nens);
  real2d hmean_r  ("hmean_r",   nz,nens);
  real2d hmean_q  ("hmean_q",   nz,nens);
  real2d hmean_u  ("hmean_u",   nz,nens);
  real2d hmean_v  ("hmean_v",   nz,nens);
  real2d ttend_acc("ttend_acc", nz,nens);
  real2d rtend_acc("rtend_acc", nz,nens);
  real2d qtend_acc("qtend_acc", nz,nens);
  real2d utend_acc("utend_acc", nz,nens);
  real2d vtend_acc("vtend_acc", nz,nens);
  real2d qpoz     ("qpoz",      nz,nens);
  real2d qneg     ("qneg",      nz,nens);
  //------------------------------------------------------------------------------------------------
  real constexpr dtemp_max = 5; // temperature tendency max threshold =>  5 K following UP-CAM
  real constexpr temp_min = 50; // temperature minimum minthreshold   => 50 K following UP-CAM
  //------------------------------------------------------------------------------------------------
  // Compute the horizontal mean for each variable
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int n) {
    hmean_t(k,n) = 0.0;
    hmean_r(k,n) = 0.0;
    hmean_q(k,n) = 0.0;
    if (crm_accel_uv) {
      hmean_u(k,n) = 0.0;
      hmean_v(k,n) = 0.0;
    }
  });
  real r_nx_ny  = 1._fp/(nx*ny);  // precompute reciprocal to avoid costly divisions
  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
    yakl::atomicAdd( hmean_t(k,n), temp(k,j,i,n) * r_nx_ny );
    yakl::atomicAdd( hmean_r(k,n), rhod(k,j,i,n) * r_nx_ny );
    yakl::atomicAdd( hmean_q(k,n), ( rhov(k,j,i,n) + rhol(k,j,i,n) + rhoi(k,j,i,n) ) * r_nx_ny );
    if (crm_accel_uv) {
      yakl::atomicAdd( hmean_u(k,n), uvel(k,j,i,n) * r_nx_ny );
      yakl::atomicAdd( hmean_v(k,n), vvel(k,j,i,n) * r_nx_ny );
    }
  });
  //------------------------------------------------------------------------------------------------
  // Compute the accelerated tendencies
  // NOTE - these are tendencies multiplied by the time step (i.e. d/dt * crm_dt)
  ScalarLiveOut<bool> ceaseflag_liveout(false);
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int n) {
    ttend_acc(k,n) = hmean_t(k,n) - accel_save_t(k,n);
    rtend_acc(k,n) = hmean_r(k,n) - accel_save_r(k,n);
    qtend_acc(k,n) = hmean_q(k,n) - accel_save_q(k,n);
    if (crm_accel_uv) {
      utend_acc(k,n) = hmean_u(k,n) - accel_save_u(k,n);
      vtend_acc(k,n) = hmean_v(k,n) - accel_save_v(k,n);
    }
    if (abs(ttend_acc(k,n)) > dtemp_max) {
      ceaseflag_liveout = true;
    }
  });
  bool ceaseflag = ceaseflag_liveout.hostRead();
  //------------------------------------------------------------------------------------------------
  // If acceleration tendencies are insane then just abort the acceleration
  if (ceaseflag) {
    // When temperature tendency threshold is triggered the acceleration will
    // not be applied for the remainder of the CRM integration.
    // The number of CRM steps for this integration (nstop) must be  updated
    // to ensure the CRM integration duration is unchanged.
    // 
    // The effective MSA timestep is dt_a = crm_dt * (1 + crm_accel_factor). When
    // ceaseflag is triggered at nstep, we've taken (nstep - 1) previous steps of
    // size crm_dt * (1 + crm_accel_factor). The current step, and all future
    // steps, will revert to size crm_dt. Therefore, the total crm integration
    // time remaining after this step is
    //     time_remaining = crm_run_time - (nstep - 1)* dt_a + crm_dt
    //     nsteps_remaining = time_remaining / crm_dt
    //     updated nstop = nstep + nsteps_remaining
    // Because we set nstop = crm_run_time / dt_a in crm_accel_nstop, subbing
    // crm_run_time = nstop * dt_a and working through algebra yields 
    //     updated nstop = nstop + (nstop - nstep + 1) * crm_accel_factor.
    // This only can happen once!
    coupler.set_option<bool>("crm_acceleration_ceaseflag",true);
    int nstop_old = nstop;
    int nstop_new = nstop_old + (nstop_old - nstep + 1)*crm_accel_factor;
    printf("pam_accelerate: MSA tendencies too large - nstop increased from %d to %d \n",nstop_old, nstop_new);
    nstop = nstop_new; 
    return;
  }
  //------------------------------------------------------------------------------------------------
  // Apply the accelerated tendencies
  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
    
    temp(k,j,i,n) = temp(k,j,i,n) + crm_accel_factor * ttend_acc(k,n);
    rhod(k,j,i,n) = rhod(k,j,i,n) + crm_accel_factor * rtend_acc(k,n);
    rhov(k,j,i,n) = rhov(k,j,i,n) + crm_accel_factor * qtend_acc(k,n);
    if (crm_accel_uv) {
      uvel(k,j,i,n) = uvel(k,j,i,n) + crm_accel_factor * utend_acc(k,n); 
      vvel(k,j,i,n) = vvel(k,j,i,n) + crm_accel_factor * vtend_acc(k,n); 
    }
    // apply limiter on temperature
    temp(k,j,i,n) = std::max( temp_min, temp(k,j,i,n) );
  });
  //------------------------------------------------------------------------------------------------
  // Evaluate and address instances of negative water vapor
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int n) {
    qpoz(k,n) = 0.0;
    qneg(k,n) = 0.0;
  });
  // accumulate positive and negative water vapor density values in each layer
  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
    if (rhov(k,j,i,n) < 0.0) {
      yakl::atomicAdd( qneg(k,n) , rhov(k,j,i,n) ); 
    } else {
      yakl::atomicAdd( qpoz(k,n) , rhov(k,j,i,n) );
    }
  });
  // clip or adjust points with negative water vapor density
  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
    real factor;
    // check if all moisture is depleted in the layer
    if (qpoz(k,n) + qneg(k,n) <= 0.0) {
      rhov(k,j,i,n) = 0.0;
    } else {
      // Clip vapor density values at 0 and remove the negative excess in each layer
      // proportionally from the positive qt fields in the layer
      factor = 1.0 + qneg(k,n) / qpoz(k,n);
      rhov(k,j,i,n) = std::max(0.0, rhov(k,j,i,n) * factor);
    } 
  });
  //------------------------------------------------------------------------------------------------
}

