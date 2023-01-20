#pragma once

#include "pam_coupler.h"

// These routines are used to encapsulate the agggregation 
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
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  //------------------------------------------------------------------------------------------------
  // aggregted quantities
  dm_device.register_and_allocate<real>("precip_liq_aggregated","aggregated sfc liq precip rate",{nens},{"nens"});
  dm_device.register_and_allocate<real>("precip_ice_aggregated","aggregated sfc ice precip rate",{nens},{"nens"});
  dm_device.register_and_allocate<real>("stat_aggregation_cnt","number of aggregated samples",{nens},{"nens"});
  //------------------------------------------------------------------------------------------------
  auto precip_liq_aggregated = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice_aggregated = dm_device.get<real,1>("precip_ice_aggregated");
  auto stat_aggregation_cnt  = dm_device.get<real,1>("stat_aggregation_cnt");
  parallel_for("Initialize aggregated precipitation", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    precip_liq_aggregated(iens) = 0;
    precip_ice_aggregated(iens) = 0;
    stat_aggregation_cnt(iens) = 0;
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
  //------------------------------------------------------------------------------------------------
  auto precip_liq = dm_device.get<real,3>("precip_liq_surf_out");
  auto precip_ice = dm_device.get<real,3>("precip_ice_surf_out");
  auto precip_liq_aggregated = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice_aggregated = dm_device.get<real,1>("precip_ice_aggregated");
  auto stat_aggregation_cnt  = dm_device.get<real,1>("stat_aggregation_cnt");
  // aggregate surface precipitation
  real r_nx_ny  = 1._fp / (nx*ny);
  parallel_for("aggregate statistics", SimpleBounds<3>(ny,nx,nens), YAKL_LAMBDA (int j, int i, int iens) {
    atomicAdd( precip_liq_aggregated(iens), precip_liq(j,i,iens)*r_nx_ny );
    atomicAdd( precip_ice_aggregated(iens), precip_ice(j,i,iens)*r_nx_ny );
  });
  parallel_for("update statistics aggregation count", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    stat_aggregation_cnt(iens) = stat_aggregation_cnt(iens) + 1;
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
  auto nz         = coupler.get_option<int>("crm_nz");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  //------------------------------------------------------------------------------------------------
  // convert aggregated values to time means
  auto precip_liq      = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice      = dm_device.get<real,1>("precip_ice_aggregated");
  auto aggregation_cnt = dm_device.get<real,1>("stat_aggregation_cnt");
  real1d precip_tot("precip_tot",nens);
  parallel_for("calculate total precipitation", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    precip_liq(iens) = precip_liq(iens) / aggregation_cnt(iens);
    precip_ice(iens) = precip_ice(iens) / aggregation_cnt(iens);
    precip_tot(iens) = precip_liq(iens) + precip_ice(iens);
  });
  //------------------------------------------------------------------------------------------------
  // copy data to host
  auto precip_tot_host = dm_host.get<real,1>("output_precc");
  auto precip_ice_host = dm_host.get<real,1>("output_precsc");
  precip_tot.deep_copy_to(precip_tot_host);
  precip_ice.deep_copy_to(precip_ice_host);
  //------------------------------------------------------------------------------------------------
}

