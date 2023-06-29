#pragma once

#include "pam_coupler.h"

// Copy the CRM radiation tendencies into the PAM coupler
inline void pam_radiation_copy_input_to_coupler( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto crm_ny     = coupler.get_option<int>("crm_ny");
  auto crm_nx     = coupler.get_option<int>("crm_nx");
  auto rad_ny     = coupler.get_option<int>("rad_ny");
  auto rad_nx     = coupler.get_option<int>("rad_nx");
  auto cp_d       = coupler.get_option<real>("cp_d");
  //------------------------------------------------------------------------------------------------
  // get the coupler rad tendency variable
  auto rad_enthalpy_tend = dm_device.get<real,4>("rad_enthalpy_tend");
  //------------------------------------------------------------------------------------------------
  // wrap the host CRM state data in YAKL arrays
  auto gcm_qrad = dm_host.get<real const,4>("rad_qrad").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // Copy the host CRM data to the coupler
  parallel_for("copy rad tendencies to CRM", SimpleBounds<4>(nz,rad_ny,rad_nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    rad_enthalpy_tend(k,j,i,iens) = gcm_qrad(k,j,i,iens)*cp_d;
  });
  //------------------------------------------------------------------------------------------------
}


// register and initialize various quantities for radiation
inline void pam_radiation_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm = coupler.get_data_manager_device_readwrite();
  auto nens   = coupler.get_option<int>("ncrms");
  auto nz     = coupler.get_option<int>("crm_nz");
  auto crm_ny = coupler.get_option<int>("crm_ny");
  auto crm_nx = coupler.get_option<int>("crm_nx");
  auto rad_ny = coupler.get_option<int>("rad_ny");
  auto rad_nx = coupler.get_option<int>("rad_nx");
  //------------------------------------------------------------------------------------------------
  real rad_nx_fac;
  real rad_ny_fac;
  if ( crm_nx%rad_nx==0 || crm_ny%rad_ny==0  ) {
    rad_nx_fac = static_cast<real>(rad_nx)/static_cast<real>(crm_nx);
    rad_ny_fac = static_cast<real>(rad_ny)/static_cast<real>(crm_ny);
  } else {
    std::cout << "crm_nx_rad and crm_ny_rad need to be divisible by nx and ny";
    exit(-1);
  }
  coupler.set_option<real>("rad_nx_fac",rad_nx_fac);
  coupler.set_option<real>("rad_ny_fac",rad_ny_fac);
  //------------------------------------------------------------------------------------------------
  // register aggregted quantities
  dm.register_and_allocate<real>("rad_aggregation_cnt","number of aggregated samples",{nens},{"nens"});
  dm.register_and_allocate<real>("rad_temperature","rad column mean temperature",      {nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  dm.register_and_allocate<real>("rad_qv"         ,"rad column mean water vapor",      {nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  dm.register_and_allocate<real>("rad_qc"         ,"rad column mean cloud liq amount", {nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  dm.register_and_allocate<real>("rad_qi"         ,"rad column mean cloud ice amount", {nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  dm.register_and_allocate<real>("rad_nc"         ,"rad column mean cloud liq number", {nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  dm.register_and_allocate<real>("rad_ni"         ,"rad column mean cloud ice number", {nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  dm.register_and_allocate<real>("rad_cld"        ,"rad column mean cloud fraction",   {nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  //------------------------------------------------------------------------------------------------
  // initialize aggregted quantities
  auto rad_aggregation_cnt = dm.get<real,1>("rad_aggregation_cnt");
  auto rad_temperature     = dm.get<real,4>("rad_temperature");
  auto rad_qv              = dm.get<real,4>("rad_qv");
  auto rad_qc              = dm.get<real,4>("rad_qc");
  auto rad_qi              = dm.get<real,4>("rad_qi");
  auto rad_nc              = dm.get<real,4>("rad_nc");
  auto rad_ni              = dm.get<real,4>("rad_ni");
  auto rad_cld             = dm.get<real,4>("rad_cld");
  parallel_for("Initialize output for radiation", SimpleBounds<4>(nz,rad_ny,rad_nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    rad_temperature(k,j,i,iens) = 0.;
    rad_qv         (k,j,i,iens) = 0.;
    rad_qc         (k,j,i,iens) = 0.;
    rad_qi         (k,j,i,iens) = 0.;
    rad_nc         (k,j,i,iens) = 0.;
    rad_ni         (k,j,i,iens) = 0.;
    rad_cld        (k,j,i,iens) = 0.;
  });
  parallel_for("update radiation aggregation count", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    rad_aggregation_cnt(iens) = 0.;
  });
  //------------------------------------------------------------------------------------------------
}


// aggregate rad quanties for each PAM time step
inline void pam_radiation_timestep_aggregation( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  using yakl::intrinsics::maxval;
  auto &dm = coupler.get_data_manager_device_readwrite();
  auto nens   = coupler.get_option<int>("ncrms");
  auto nz     = coupler.get_option<int>("crm_nz");
  auto crm_ny = coupler.get_option<int>("crm_ny");
  auto crm_nx = coupler.get_option<int>("crm_nx");
  auto rad_ny = coupler.get_option<int>("rad_ny");
  auto rad_nx = coupler.get_option<int>("rad_nx");
  //------------------------------------------------------------------------------------------------
  // Get current CRM state
  auto temp    = dm.get<real const,4>("temp"           );
  auto rho_d   = dm.get<real const,4>("density_dry"    );
  auto rho_v   = dm.get<real const,4>("water_vapor"    );
  auto rho_l   = dm.get<real const,4>("cloud_water"    );
  auto rho_i   = dm.get<real const,4>("ice"            );
  auto num_l   = dm.get<real const,4>("cloud_water_num");
  auto num_i   = dm.get<real const,4>("ice_num"        );
  //------------------------------------------------------------------------------------------------
  // Get aggregated rad variables
  auto rad_temperature     = dm.get<real,4>("rad_temperature");
  auto rad_qv              = dm.get<real,4>("rad_qv");
  auto rad_qc              = dm.get<real,4>("rad_qc");
  auto rad_qi              = dm.get<real,4>("rad_qi");
  auto rad_nc              = dm.get<real,4>("rad_nc");
  auto rad_ni              = dm.get<real,4>("rad_ni");
  auto rad_cld             = dm.get<real,4>("rad_cld");
  auto rad_aggregation_cnt = dm.get<real,1>("rad_aggregation_cnt");
  //------------------------------------------------------------------------------------------------
  // aggregate radiation column data
  auto rad_nx_fac = coupler.get_option<real>("rad_nx_fac");
  auto rad_ny_fac = coupler.get_option<real>("rad_ny_fac");
  real r_nx_ny  = rad_nx_fac * rad_ny_fac;
  parallel_for("aggregate rad state", SimpleBounds<4>(nz,crm_ny,crm_nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    int i_rad = i / (crm_nx/rad_nx);
    int j_rad = j / (crm_nx/rad_ny);
    real rho_total = rho_d(k,j,i,iens) + rho_v(k,j,i,iens);
    atomicAdd( rad_temperature(k,j_rad,i_rad,iens), temp(k,j,i,iens)                 * r_nx_ny );
    atomicAdd( rad_qv         (k,j_rad,i_rad,iens), std::max(0.0,rho_v(k,j,i,iens)/rho_total) * r_nx_ny );
    atomicAdd( rad_qc         (k,j_rad,i_rad,iens), std::max(0.0,rho_l(k,j,i,iens)/rho_total) * r_nx_ny );
    atomicAdd( rad_qi         (k,j_rad,i_rad,iens), std::max(0.0,rho_i(k,j,i,iens)/rho_total) * r_nx_ny );
    atomicAdd( rad_nc         (k,j_rad,i_rad,iens), num_l(k,j,i,iens)                * r_nx_ny );
    atomicAdd( rad_ni         (k,j_rad,i_rad,iens), num_i(k,j,i,iens)                * r_nx_ny );
    if ( rho_l(k,j,i,iens) + rho_i(k,j,i,iens ) > 0) {
      atomicAdd( rad_cld(k,j_rad,i_rad,iens), 1.* r_nx_ny );
    }
  });
  parallel_for("update radiation aggregation count", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    rad_aggregation_cnt(iens) += 1;
  });
  //------------------------------------------------------------------------------------------------
}


// convert aggregated quantites to means
inline void pam_radiation_compute_means( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto crm_ny     = coupler.get_option<int>("crm_ny");
  auto crm_nx     = coupler.get_option<int>("crm_nx");
  auto rad_ny     = coupler.get_option<int>("rad_ny");
  auto rad_nx     = coupler.get_option<int>("rad_nx");
  //------------------------------------------------------------------------------------------------
  // get the coupler rad tendency variable
  auto rad_temperature     = dm_device.get<real,4>("rad_temperature");
  auto rad_qv              = dm_device.get<real,4>("rad_qv");
  auto rad_qc              = dm_device.get<real,4>("rad_qc");
  auto rad_qi              = dm_device.get<real,4>("rad_qi");
  auto rad_nc              = dm_device.get<real,4>("rad_nc");
  auto rad_ni              = dm_device.get<real,4>("rad_ni");
  auto rad_cld             = dm_device.get<real,4>("rad_cld");
  auto rad_aggregation_cnt = dm_device.get<real,1>("rad_aggregation_cnt");
  //------------------------------------------------------------------------------------------------
  // Convert sum to time mean and copy the CRM data to the GCM
  parallel_for("copy rad state to GCM", SimpleBounds<4>(nz,rad_ny,rad_nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    rad_temperature(k,j,i,iens) = rad_temperature(k,j,i,iens) / rad_aggregation_cnt(iens);
    rad_qv         (k,j,i,iens) = rad_qv         (k,j,i,iens) / rad_aggregation_cnt(iens);
    rad_qc         (k,j,i,iens) = rad_qc         (k,j,i,iens) / rad_aggregation_cnt(iens);
    rad_qi         (k,j,i,iens) = rad_qi         (k,j,i,iens) / rad_aggregation_cnt(iens);
    rad_nc         (k,j,i,iens) = rad_nc         (k,j,i,iens) / rad_aggregation_cnt(iens);
    rad_ni         (k,j,i,iens) = rad_ni         (k,j,i,iens) / rad_aggregation_cnt(iens);
    rad_cld        (k,j,i,iens) = rad_cld        (k,j,i,iens) / rad_aggregation_cnt(iens);
  });
  //------------------------------------------------------------------------------------------------
}

// Copy the CRM radiation state into the PAM coupler
inline void pam_radiation_copy_to_host( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  // get the coupler rad tendency variable
  auto rad_temperature     = dm_device.get<real,4>("rad_temperature");
  auto rad_qv              = dm_device.get<real,4>("rad_qv");
  auto rad_qc              = dm_device.get<real,4>("rad_qc");
  auto rad_qi              = dm_device.get<real,4>("rad_qi");
  auto rad_nc              = dm_device.get<real,4>("rad_nc");
  auto rad_ni              = dm_device.get<real,4>("rad_ni");
  auto rad_cld             = dm_device.get<real,4>("rad_cld");
  auto rad_aggregation_cnt = dm_device.get<real,1>("rad_aggregation_cnt");
  //------------------------------------------------------------------------------------------------
  // copy rad column data to host
  auto gcm_rad_temperature = dm_host.get<real,4>("rad_temperature");
  auto gcm_rad_qv          = dm_host.get<real,4>("rad_qv");
  auto gcm_rad_qc          = dm_host.get<real,4>("rad_qc");
  auto gcm_rad_qi          = dm_host.get<real,4>("rad_qi");
  auto gcm_rad_nc          = dm_host.get<real,4>("rad_nc");
  auto gcm_rad_ni          = dm_host.get<real,4>("rad_ni");
  auto gcm_rad_cld         = dm_host.get<real,4>("rad_cld");
  rad_temperature.deep_copy_to(gcm_rad_temperature);
  rad_qv         .deep_copy_to(gcm_rad_qv         );
  rad_qc         .deep_copy_to(gcm_rad_qc         );
  rad_qi         .deep_copy_to(gcm_rad_qi         );
  rad_nc         .deep_copy_to(gcm_rad_nc         );
  rad_ni         .deep_copy_to(gcm_rad_ni         );
  rad_cld        .deep_copy_to(gcm_rad_cld        );
  yakl::fence();
  //------------------------------------------------------------------------------------------------
}