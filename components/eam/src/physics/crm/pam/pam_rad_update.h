#pragma once

#include "pam_coupler.h"

// Copy the CRM radiation tendencies into the PAM coupler
inline void copy_input_rad_to_coupler( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int nens = dm_device.get_dimension_size("nens");
  int nz   = dm_device.get_dimension_size("z"   );
  int ny   = dm_device.get_dimension_size("y"   );
  int nx   = dm_device.get_dimension_size("x"   );
  auto ny_rad = coupler.get_option<int>("ny_rad");
  auto nx_rad = coupler.get_option<int>("nx_rad");
  //------------------------------------------------------------------------------------------------
  // get the coupler rad tendency variable
  auto rad_enthalpy_tend = dm_device.get<real,4>("crm_rad_tend");
  //------------------------------------------------------------------------------------------------
  // wrap the host CRM state data in YAKL arrays
  auto gcm_qrad = dm_host.get<real const,4>("rad_qrad").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // Copy the host CRM data to the coupler
  parallel_for("copy rad tendencies to CRM", SimpleBounds<4>(nz,ny_rad,nx_rad,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    rad_enthalpy_tend(k,j,i,iens) = gcm_qrad(k,j,i,iens);
  });
  //------------------------------------------------------------------------------------------------
}


// Copy the CRM radiation state into the PAM coupler
inline void copy_output_rad_to_gcm( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int nens = dm_device.get_dimension_size("nens");
  int nz   = dm_device.get_dimension_size("z"   );
  int ny   = dm_device.get_dimension_size("y"   );
  int nx   = dm_device.get_dimension_size("x"   );
  auto ny_rad = coupler.get_option<int>("ny_rad");
  auto nx_rad = coupler.get_option<int>("nx_rad");
  //------------------------------------------------------------------------------------------------
  // get the coupler rad tendency variable
  auto crm_rad_temperature = dm_device.get<real const,4>("crm_rad_temperature");
  auto crm_rad_qv          = dm_device.get<real const,4>("crm_rad_qv");
  auto crm_rad_qc          = dm_device.get<real const,4>("crm_rad_qc");
  auto crm_rad_qi          = dm_device.get<real const,4>("crm_rad_qi");
  auto crm_rad_cld         = dm_device.get<real const,4>("crm_rad_cld");
  //------------------------------------------------------------------------------------------------
  // wrap the host CRM state data in YAKL arrays
  auto gcm_rad_temperature = dm_host.get<real,4>("rad_temperature").createDeviceCopy();
  auto gcm_rad_qv          = dm_host.get<real,4>("rad_qv").createDeviceCopy();
  auto gcm_rad_qc          = dm_host.get<real,4>("rad_qc").createDeviceCopy();
  auto gcm_rad_qi          = dm_host.get<real,4>("rad_qi").createDeviceCopy();
  auto gcm_rad_cld         = dm_host.get<real,4>("rad_cld").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // Copy the CRM data to the GCM
  parallel_for("copy rad state to GCM", SimpleBounds<4>(nz,ny_rad,nx_rad,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    int i_rad = i / (nx/nx_rad);
    int j_rad = j / (nx/ny_rad);
    gcm_rad_temperature(k,j,i,iens) = crm_rad_temperature(k,j,i,iens);
    gcm_rad_qv         (k,j,i,iens) = crm_rad_qv         (k,j,i,iens);
    gcm_rad_qc         (k,j,i,iens) = crm_rad_qc         (k,j,i,iens);
    gcm_rad_qi         (k,j,i,iens) = crm_rad_qi         (k,j,i,iens);
    gcm_rad_cld        (k,j,i,iens) = crm_rad_cld        (k,j,i,iens);
  });
  //------------------------------------------------------------------------------------------------
}