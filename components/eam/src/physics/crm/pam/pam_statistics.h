#pragma once

#include "pam_coupler.h"
#include "saturation_adjustment.h"

// These routines are used to encapsulate the aggregation
// of various quantities, such as precipitation


// register and initialize various quantities for statistical calculations
inline void pam_statistics_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  //------------------------------------------------------------------------------------------------
  // aggregted quantities
  dm_device.register_and_allocate<real>("stat_aggregation_cnt", "number of aggregated samples",  {nens},{"nens"});
  dm_device.register_and_allocate<real>("precip_liq_aggregated","aggregated sfc liq precip rate",{nens},{"nens"});
  dm_device.register_and_allocate<real>("precip_ice_aggregated","aggregated sfc ice precip rate",{nens},{"nens"});
  dm_device.register_and_allocate<real>("cldfrac_aggregated",   "aggregated cloud fraction",     {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("clear_rh"       ,      "clear air rel humidity",        {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("clear_rh_cnt"   ,      "clear air count",               {nz,nens},{"z","nens"});
  //------------------------------------------------------------------------------------------------
  auto stat_aggregation_cnt  = dm_device.get<real,1>("stat_aggregation_cnt");
  auto precip_liq_aggregated = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice_aggregated = dm_device.get<real,1>("precip_ice_aggregated");
  auto cldfrac_aggregated    = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh              = dm_device.get<real,2>("clear_rh");
  auto clear_rh_cnt          = dm_device.get<real,2>("clear_rh_cnt");
  parallel_for("Initialize aggregated precipitation", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    stat_aggregation_cnt(iens)  = 0;
    precip_liq_aggregated(iens) = 0;
    precip_ice_aggregated(iens) = 0;
  });
  parallel_for("Initialize aggregated precipitation", SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k, int iens) {
    cldfrac_aggregated(k,iens) = 0;
    clear_rh          (k,iens) = 0;
    clear_rh_cnt      (k,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
}


// aggregate quanties for each PAM time step
inline void pam_statistics_timestep_aggregation( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  real R_v        = coupler.get_option<real>("R_v");
  //------------------------------------------------------------------------------------------------
  // get CRM variables to be aggregated
  auto precip_liq = dm_device.get<real,3>("precip_liq_surf_out");
  auto precip_ice = dm_device.get<real,3>("precip_ice_surf_out");
  auto temp       = dm_device.get<real,4>("temp"       );
  auto rho_v      = dm_device.get<real,4>("water_vapor");
  auto rho_l      = dm_device.get<real,4>("cloud_water");
  auto rho_i      = dm_device.get<real,4>("ice"        );
  //------------------------------------------------------------------------------------------------
  // get aggregation variables
  auto stat_aggregation_cnt  = dm_device.get<real,1>("stat_aggregation_cnt");
  auto precip_liq_aggregated = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice_aggregated = dm_device.get<real,1>("precip_ice_aggregated");
  auto cldfrac_aggregated    = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh              = dm_device.get<real,2>("clear_rh");
  auto clear_rh_cnt          = dm_device.get<real,2>("clear_rh_cnt");
  //------------------------------------------------------------------------------------------------
  // perform aggregation
  parallel_for("update aggregation count", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    stat_aggregation_cnt(iens) = stat_aggregation_cnt(iens) + 1;
  });
  real r_nx_ny  = 1._fp / (nx*ny);
  parallel_for("aggregate 0D statistics", SimpleBounds<3>(ny,nx,nens), YAKL_LAMBDA (int j, int i, int iens) {
    atomicAdd( precip_liq_aggregated(iens), precip_liq(j,i,iens)*r_nx_ny );
    atomicAdd( precip_ice_aggregated(iens), precip_ice(j,i,iens)*r_nx_ny );
  });
  parallel_for("aggregate 1D statistics", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    if ( rho_l(k,j,i,iens) + rho_i(k,j,i,iens ) > 0) {
      atomicAdd( cldfrac_aggregated(k,iens), 1.*r_nx_ny );
    } else {
      // calculate RH from vapor pressure and saturation vapor pressure
      real pv = rho_v(k,j,i,iens) * R_v * temp(k,j,i,iens);
      real svp = modules::saturation_vapor_pressure( temp(k,j,i,iens) );
      atomicAdd( clear_rh    (k,iens), pv/svp );
      atomicAdd( clear_rh_cnt(k,iens), 1. );
    }
  });
  //------------------------------------------------------------------------------------------------
}


// copy aggregated statistical quantities to host
inline void pam_statistics_compute_means( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto crm_nz     = coupler.get_option<int>("crm_nz");
  auto gcm_nlev   = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  // convert aggregated values to time means
  auto aggregation_cnt = dm_device.get<real,1>("stat_aggregation_cnt");
  auto precip_liq      = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice      = dm_device.get<real,1>("precip_ice_aggregated");
  auto cldfrac         = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh        = dm_device.get<real,2>("clear_rh");
  auto clear_rh_cnt    = dm_device.get<real,2>("clear_rh_cnt");
  parallel_for("finalize aggregated variables", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    precip_liq(iens) = precip_liq(iens) / aggregation_cnt(iens);
    precip_ice(iens) = precip_ice(iens) / aggregation_cnt(iens);
  });
  parallel_for("Initialize aggregated precipitation", SimpleBounds<2>(crm_nz,nens), YAKL_LAMBDA (int k, int iens) {
    cldfrac(k,iens)  = cldfrac(k,iens)  / aggregation_cnt(iens);
    if (clear_rh_cnt(k,iens)>0) {
      clear_rh(k,iens) = clear_rh(k,iens) / clear_rh_cnt(k,iens);
    }
  });
  //------------------------------------------------------------------------------------------------
}


// copy aggregated statistical quantities to host
inline void pam_statistics_copy_to_host( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto crm_nz     = coupler.get_option<int>("crm_nz");
  auto gcm_nlev   = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  // convert aggregated values to time means
  auto precip_liq      = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice      = dm_device.get<real,1>("precip_ice_aggregated");
  auto cldfrac         = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh        = dm_device.get<real,2>("clear_rh");
  //------------------------------------------------------------------------------------------------
  // calculate total precip
  real1d precip_tot("precip_tot",nens);
  parallel_for("finalize aggregated variables", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    precip_tot(iens) = precip_liq(iens) + precip_ice(iens);
  });
  //------------------------------------------------------------------------------------------------
  // convert variables to GCM vertical grid
  real2d cldfrac_gcm("cldfrac_gcm",gcm_nlev,nens);
  parallel_for("Initialize aggregated precipitation", SimpleBounds<2>(gcm_nlev,nens), YAKL_LAMBDA (int k_gcm, int iens) {
    int k_crm = gcm_nlev-1-k_gcm;
    if (k_crm<crm_nz) {
      cldfrac_gcm(k_gcm,iens) = cldfrac(k_crm,iens);
    } else {
      cldfrac_gcm(k_gcm,iens) = 0.;
    }
  });
  //------------------------------------------------------------------------------------------------
  // copy data to host
  auto precip_tot_host = dm_host.get<real,1>("output_precc");
  auto precip_ice_host = dm_host.get<real,1>("output_precsc");
  auto cldfrac_host    = dm_host.get<real,2>("output_cld");
  auto clear_rh_host   = dm_host.get<real,2>("output_clear_rh");
  precip_tot    .deep_copy_to(precip_tot_host);
  precip_ice    .deep_copy_to(precip_ice_host);
  cldfrac_gcm   .deep_copy_to(cldfrac_host);
  clear_rh      .deep_copy_to(clear_rh_host);
  //------------------------------------------------------------------------------------------------
}

