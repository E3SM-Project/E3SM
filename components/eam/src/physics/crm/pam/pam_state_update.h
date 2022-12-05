#pragma once

#include "pam_coupler.h"

// update the coupler GCM state variables using the input GCM state
inline void update_gcm_state( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int nz   = dm_device.get_dimension_size("z"   );
  int ny   = dm_device.get_dimension_size("y"   );
  int nx   = dm_device.get_dimension_size("x"   );
  int nens = dm_device.get_dimension_size("nens");
  //------------------------------------------------------------------------------------------------
  // get the coupler GCM state arrays used to force the CRM
  real Lv = coupler.get_option<real>("latvap") ;
  real Lf = coupler.get_option<real>("latice") ;
  auto gcm_rho_d = dm_device.get<real,2>("gcm_density_dry");
  auto gcm_uvel  = dm_device.get<real,2>("gcm_uvel"       );
  auto gcm_vvel  = dm_device.get<real,2>("gcm_vvel"       );
  auto gcm_temp  = dm_device.get<real,2>("gcm_temp"       );
  auto gcm_rho_v = dm_device.get<real,2>("gcm_water_vapor");
  //------------------------------------------------------------------------------------------------
  // wrap the host GCM state data in YAKL arrays
  auto input_ul   = dm_host.get<real const,2>("input_ul").createDeviceCopy();
  auto input_vl   = dm_host.get<real const,2>("input_vl").createDeviceCopy();
  auto input_tl   = dm_host.get<real const,2>("input_tl").createDeviceCopy();
  auto input_qccl = dm_host.get<real const,2>("input_qccl").createDeviceCopy();
  auto input_qiil = dm_host.get<real const,2>("input_qiil").createDeviceCopy();
  auto input_ql   = dm_host.get<real const,2>("input_ql").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // Define GCM state for forcing - adjusted to avoid directly forcing cloud liquid and ice fields
  parallel_for( Bounds<2>(crm_nz,nens) , YAKL_LAMBDA (int k, int iens) {
    gcm_rho_d(k,iens) = input_pmid(k,iens) / ( micro.R_d * input_tl(k,iens) );
    gcm_uvel (k,iens) = input_ul(k,iens);
    gcm_vvel (k,iens) = input_vl(k,iens);
    // convert total water mixing ratio to water vapor density
    real input_qt = input_ql(k,iens) + input_qccl(k,iens) + input_qiil(k,iens);
    gcm_rho_v(k,iens) = input_qt * gcm_rho_d(k,iens) / ( 1 - input_qt );
    // adjust temperature to account for liq/ice to vapor conversion
    real input_t_adj = input_tl(k,iens) - ( input_qccl(k,iens)*Lv + input_qiil(k,iens)*Lf ) / micro.cp_d ;
    gcm_temp(k,iens) = input_t_adj;
  });
  //------------------------------------------------------------------------------------------------
}


// Copy the CRM state saved by the GCM into the PAM coupler
inline void copy_input_state_to_coupler( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int nz   = dm_device.get_dimension_size("z"   );
  int ny   = dm_device.get_dimension_size("y"   );
  int nx   = dm_device.get_dimension_size("x"   );
  int nens = dm_device.get_dimension_size("nens");
  //------------------------------------------------------------------------------------------------
  // get the coupler state variables
  auto crm_rho_d  = dm_device.get<real,4>("density_dry");
  auto crm_uvel   = dm_device.get<real,4>("uvel");
  auto crm_vvel   = dm_device.get<real,4>("vvel");
  auto crm_wvel   = dm_device.get<real,4>("wvel");
  auto crm_temp   = dm_device.get<real,4>("temp");
  auto crm_qv     = dm_device.get<real,4>("water_vapor");
  auto crm_qc     = dm_device.get<real,4>("cloud_water");
  auto crm_nc     = dm_device.get<real,4>("cloud_water_num");
  auto crm_qr     = dm_device.get<real,4>("rain");
  auto crm_nr     = dm_device.get<real,4>("rain_num");
  auto crm_qi     = dm_device.get<real,4>("ice");
  auto crm_ni     = dm_device.get<real,4>("ice_num");
  auto crm_qm     = dm_device.get<real,4>("ice_rime");
  auto crm_bm     = dm_device.get<real,4>("ice_rime_vol");
  auto crm_t_prev = dm_device.get<real,4>("qv_prev");
  auto crm_q_prev = dm_device.get<real,4>("t_prev");
  //------------------------------------------------------------------------------------------------
  // wrap the host CRM state data in YAKL arrays
  auto state_u_wind       = dm_host.get<real const,2>("state_u_wind").createDeviceCopy();
  auto state_v_wind       = dm_host.get<real const,2>("state_v_wind").createDeviceCopy();
  auto state_w_wind       = dm_host.get<real const,2>("state_w_wind").createDeviceCopy();
  auto state_temperature  = dm_host.get<real const,2>("state_temperature").createDeviceCopy();
  auto state_qv           = dm_host.get<real const,2>("state_qv").createDeviceCopy();
  auto state_qc           = dm_host.get<real const,2>("state_qc").createDeviceCopy();
  auto state_nc           = dm_host.get<real const,2>("state_nc").createDeviceCopy();
  auto state_qr           = dm_host.get<real const,2>("state_qr").createDeviceCopy();
  auto state_nr           = dm_host.get<real const,2>("state_nr").createDeviceCopy();
  auto state_qi           = dm_host.get<real const,2>("state_qi").createDeviceCopy();
  auto state_ni           = dm_host.get<real const,2>("state_ni").createDeviceCopy();
  auto state_qm           = dm_host.get<real const,2>("state_qm").createDeviceCopy();
  auto state_bm           = dm_host.get<real const,2>("state_bm").createDeviceCopy();
  auto state_t_prev       = dm_host.get<real const,2>("state_t_prev").createDeviceCopy();
  auto state_q_prev       = dm_host.get<real const,2>("state_q_prev").createDeviceCopy();  
  //------------------------------------------------------------------------------------------------
  // Copy the host CRM data to the coupler
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    crm_uvel  (k,j,i,iens) = state_u_wind      (iens,i,j,k);
    crm_vvel  (k,j,i,iens) = state_v_wind      (iens,i,j,k);
    crm_wvel  (k,j,i,iens) = state_w_wind      (iens,i,j,k);
    crm_temp  (k,j,i,iens) = state_temperature (iens,i,j,k);
    crm_qv    (k,j,i,iens) = state_qv          (iens,i,j,k);
    crm_qc    (k,j,i,iens) = state_qc          (iens,i,j,k);
    crm_nc    (k,j,i,iens) = state_nc          (iens,i,j,k);
    crm_qr    (k,j,i,iens) = state_qr          (iens,i,j,k);
    crm_nr    (k,j,i,iens) = state_nr          (iens,i,j,k);
    crm_qi    (k,j,i,iens) = state_qi          (iens,i,j,k);
    crm_ni    (k,j,i,iens) = state_ni          (iens,i,j,k);
    crm_qm    (k,j,i,iens) = state_qm          (iens,i,j,k);
    crm_bm    (k,j,i,iens) = state_bm          (iens,i,j,k);
    crm_t_prev(k,j,i,iens) = state_t_prev      (iens,i,j,k);
    crm_q_prev(k,j,i,iens) = state_q_prev      (iens,i,j,k);
  });
  //------------------------------------------------------------------------------------------------
}
