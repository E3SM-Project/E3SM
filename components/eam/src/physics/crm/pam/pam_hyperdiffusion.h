#pragma once
#include "pam_coupler.h"

inline void pam_hyperdiffusion( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  using yakl::ScalarLiveOut;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto ny         = coupler.get_option<int>("crm_ny");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto dx         = coupler.get_option<real>("crm_dx");
  auto dy         = coupler.get_option<real>("crm_dy");
  auto dt         = coupler.get_option<real>("crm_dt");
  //------------------------------------------------------------------------------------------------
  auto temp = dm_device.get<real,4>("temp"           );
  auto rhod = dm_device.get<real,4>("density_dry"    );
  auto rhov = dm_device.get<real,4>("water_vapor"    );
  auto rhol = dm_device.get<real,4>("cloud_water"    );
  auto rhoi = dm_device.get<real,4>("ice"            );
  auto nliq = dm_device.get<real,4>("cloud_water_num");
  auto nice = dm_device.get<real,4>("ice_num"        );
  auto uvel = dm_device.get<real,4>("uvel"           );
  auto vvel = dm_device.get<real,4>("vvel"           );
  //------------------------------------------------------------------------------------------------
  #ifdef MMF_PAM_HDT
  real constexpr hd_timescale = MMF_PAM_HDT;  // damping time scale [sec]
  #else
  real constexpr hd_timescale = 10.0;  // damping time scale [sec]
  #endif
  //------------------------------------------------------------------------------------------------
  real4d hd_temp("hd_temp",nz,ny,nx,nens);
  real4d hd_rhod("hd_rhod",nz,ny,nx,nens);
  real4d hd_rhov("hd_rhov",nz,ny,nx,nens);
  real4d hd_rhol("hd_rhol",nz,ny,nx,nens);
  real4d hd_rhoi("hd_rhoi",nz,ny,nx,nens);
  real4d hd_nliq("hd_nliq",nz,ny,nx,nens);
  real4d hd_nice("hd_nice",nz,ny,nx,nens);
  real4d hd_uvel("hd_uvel",nz,ny,nx,nens);
  real4d hd_vvel("hd_vvel",nz,ny,nx,nens);
  //------------------------------------------------------------------------------------------------
  if (ny==1) {
    // Calculate hyperdiffusion to all fields in 2D domain configuration
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
      int im2 = i-2;
      int im1 = i-1;
      int ip1 = i+1;
      int ip2 = i+2;
      if (im2<0   ) { im2 += nx; } 
      if (im1<0   ) { im1 += nx; } 
      if (ip1>nx-1) { ip1 -= nx; } 
      if (ip2>nx-1) { ip2 -= nx; } 
      hd_temp(k,j,i,n) = temp(k,j,im2,n) - 4*temp(k,j,im1,n) + 6*temp(k,j,i,n) - 4*temp(k,j,ip1,n) + temp(k,j,ip2,n);
      hd_rhod(k,j,i,n) = rhod(k,j,im2,n) - 4*rhod(k,j,im1,n) + 6*rhod(k,j,i,n) - 4*rhod(k,j,ip1,n) + rhod(k,j,ip2,n);
      hd_rhov(k,j,i,n) = rhov(k,j,im2,n) - 4*rhov(k,j,im1,n) + 6*rhov(k,j,i,n) - 4*rhov(k,j,ip1,n) + rhov(k,j,ip2,n);
      hd_rhol(k,j,i,n) = rhol(k,j,im2,n) - 4*rhol(k,j,im1,n) + 6*rhol(k,j,i,n) - 4*rhol(k,j,ip1,n) + rhol(k,j,ip2,n);
      hd_rhoi(k,j,i,n) = rhoi(k,j,im2,n) - 4*rhoi(k,j,im1,n) + 6*rhoi(k,j,i,n) - 4*rhoi(k,j,ip1,n) + rhoi(k,j,ip2,n);
      hd_nliq(k,j,i,n) = nliq(k,j,im2,n) - 4*nliq(k,j,im1,n) + 6*nliq(k,j,i,n) - 4*nliq(k,j,ip1,n) + nliq(k,j,ip2,n);
      hd_nice(k,j,i,n) = nice(k,j,im2,n) - 4*nice(k,j,im1,n) + 6*nice(k,j,i,n) - 4*nice(k,j,ip1,n) + nice(k,j,ip2,n);
      hd_uvel(k,j,i,n) = uvel(k,j,im2,n) - 4*uvel(k,j,im1,n) + 6*uvel(k,j,i,n) - 4*uvel(k,j,ip1,n) + uvel(k,j,ip2,n);
      hd_vvel(k,j,i,n) = vvel(k,j,im2,n) - 4*vvel(k,j,im1,n) + 6*vvel(k,j,i,n) - 4*vvel(k,j,ip1,n) + vvel(k,j,ip2,n);
    });
  // } else {
    // Implement support for 3D domain here
  }
  //------------------------------------------------------------------------------------------------
  // Apply hyperdiffusion
  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
    temp(k,j,i,n) += -1 * dt * hd_temp(k,j,i,n) / ( 16 * hd_timescale );
    rhod(k,j,i,n) += -1 * dt * hd_rhod(k,j,i,n) / ( 16 * hd_timescale );
    rhov(k,j,i,n) += -1 * dt * hd_rhov(k,j,i,n) / ( 16 * hd_timescale );
    rhol(k,j,i,n) += -1 * dt * hd_rhol(k,j,i,n) / ( 16 * hd_timescale );
    rhoi(k,j,i,n) += -1 * dt * hd_rhoi(k,j,i,n) / ( 16 * hd_timescale );
    nliq(k,j,i,n) += -1 * dt * hd_nliq(k,j,i,n) / ( 16 * hd_timescale );
    nice(k,j,i,n) += -1 * dt * hd_nice(k,j,i,n) / ( 16 * hd_timescale );
    uvel(k,j,i,n) += -1 * dt * hd_uvel(k,j,i,n) / ( 16 * hd_timescale );
    vvel(k,j,i,n) += -1 * dt * hd_vvel(k,j,i,n) / ( 16 * hd_timescale );
  });
  //------------------------------------------------------------------------------------------------
}
