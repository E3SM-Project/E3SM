#pragma once

#include "pam_coupler.h"

// update the coupler GCM state variables using the input GCM state
inline void pam_state_update_gcm_state( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int nz   = dm_device.get_dimension_size("z"   );
  int ny   = dm_device.get_dimension_size("y"   );
  int nx   = dm_device.get_dimension_size("x"   );
  int nens = dm_device.get_dimension_size("nens");
  real R_d = coupler.get_option<real>("R_d");
  real cp_d = coupler.get_option<real>("cp_d");
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
  auto input_pmid = dm_host.get<real const,2>("input_pmid").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // Define GCM state for forcing - adjusted to avoid directly forcing cloud liquid and ice fields
  parallel_for( Bounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
    gcm_rho_d(k,iens) = ( input_pmid(k,iens) / ( R_d * input_tl(k,iens) ) ) * (1.-input_ql(k,iens)) ;
    gcm_uvel (k,iens) = input_ul(k,iens);
    gcm_vvel (k,iens) = input_vl(k,iens);
    // convert total water mixing ratio to water vapor density
    real input_qt = input_ql(k,iens) + input_qccl(k,iens) + input_qiil(k,iens);
    gcm_rho_v(k,iens) = input_qt * gcm_rho_d(k,iens) / ( 1 - input_qt );
    // adjust temperature to account for liq/ice to vapor conversion
    real input_t_adj = input_tl(k,iens) - ( input_qccl(k,iens)*Lv + input_qiil(k,iens)*Lf ) / cp_d ;
    gcm_temp(k,iens) = input_t_adj;
  });
  //------------------------------------------------------------------------------------------------
}


// Copy the CRM state saved by the GCM into the PAM coupler
inline void pam_state_copy_input_to_coupler( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int nz   = dm_device.get_dimension_size("z"   );
  int ny   = dm_device.get_dimension_size("y"   );
  int nx   = dm_device.get_dimension_size("x"   );
  int nens = dm_device.get_dimension_size("nens");
  //------------------------------------------------------------------------------------------------
  // get the coupler state variables
  auto crm_rho_d         = dm_device.get<real,4>("density_dry");
  auto crm_uvel          = dm_device.get<real,4>("uvel");
  auto crm_vvel          = dm_device.get<real,4>("vvel");
  auto crm_wvel          = dm_device.get<real,4>("wvel");
  auto crm_temp          = dm_device.get<real,4>("temp");
  auto crm_qv            = dm_device.get<real,4>("water_vapor");
  auto crm_qc            = dm_device.get<real,4>("cloud_water");
  auto crm_nc            = dm_device.get<real,4>("cloud_water_num");
  auto crm_qr            = dm_device.get<real,4>("rain");
  auto crm_nr            = dm_device.get<real,4>("rain_num");
  auto crm_qi            = dm_device.get<real,4>("ice");
  auto crm_ni            = dm_device.get<real,4>("ice_num");
  auto crm_qm            = dm_device.get<real,4>("ice_rime");
  auto crm_bm            = dm_device.get<real,4>("ice_rime_vol");
  auto crm_t_prev        = dm_device.get<real,4>("qv_prev");
  auto crm_q_prev        = dm_device.get<real,4>("t_prev");
  auto crm_shoc_wthv_sec = dm_device.get<real,4>("wthv_sec");
  auto crm_shoc_tk       = dm_device.get<real,4>("tk");
  auto crm_shoc_tkh      = dm_device.get<real,4>("tkh");
  auto crm_shoc_cldfrac  = dm_device.get<real,4>("cldfrac");
  auto crm_shoc_relvar   = dm_device.get<real,4>("relvar");
  //------------------------------------------------------------------------------------------------
  // wrap the host CRM state data in YAKL arrays
  auto state_u_wind        = dm_host.get<real const,4>("state_u_wind").createDeviceCopy();
  auto state_v_wind        = dm_host.get<real const,4>("state_v_wind").createDeviceCopy();
  auto state_w_wind        = dm_host.get<real const,4>("state_w_wind").createDeviceCopy();
  auto state_temperature   = dm_host.get<real const,4>("state_temperature").createDeviceCopy();
  auto state_qv            = dm_host.get<real const,4>("state_qv").createDeviceCopy();
  auto state_qc            = dm_host.get<real const,4>("state_qc").createDeviceCopy();
  auto state_nc            = dm_host.get<real const,4>("state_nc").createDeviceCopy();
  auto state_qr            = dm_host.get<real const,4>("state_qr").createDeviceCopy();
  auto state_nr            = dm_host.get<real const,4>("state_nr").createDeviceCopy();
  auto state_qi            = dm_host.get<real const,4>("state_qi").createDeviceCopy();
  auto state_ni            = dm_host.get<real const,4>("state_ni").createDeviceCopy();
  auto state_qm            = dm_host.get<real const,4>("state_qm").createDeviceCopy();
  auto state_bm            = dm_host.get<real const,4>("state_bm").createDeviceCopy();
  auto state_t_prev        = dm_host.get<real const,4>("state_t_prev").createDeviceCopy();
  auto state_q_prev        = dm_host.get<real const,4>("state_q_prev").createDeviceCopy();  
  auto state_shoc_wthv_sec = dm_host.get<real const,4>("state_shoc_wthv_sec").createDeviceCopy();  
  auto state_shoc_tk       = dm_host.get<real const,4>("state_shoc_tk").createDeviceCopy();  
  auto state_shoc_tkh      = dm_host.get<real const,4>("state_shoc_tkh").createDeviceCopy();  
  auto state_shoc_cldfrac  = dm_host.get<real const,4>("state_shoc_cldfrac").createDeviceCopy();  
  auto state_shoc_relvar   = dm_host.get<real const,4>("state_shoc_relvar").createDeviceCopy();  
  //------------------------------------------------------------------------------------------------
  // Copy the host CRM data to the coupler
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    // convert specific mixing ratios to density
    real rho_v = state_qv(k,j,i,iens) * crm_rho_d(k,j,i,iens) / ( 1 - state_qv(k,j,i,iens) ) ;
    real rho_c = state_qc(k,j,i,iens) * ( crm_rho_d(k,j,i,iens) + rho_v ) ;
    real rho_r = state_qr(k,j,i,iens) * ( crm_rho_d(k,j,i,iens) + rho_v ) ;
    real rho_i = state_qi(k,j,i,iens) * ( crm_rho_d(k,j,i,iens) + rho_v ) ;

    crm_uvel         (k,j,i,iens) = state_u_wind       (k,j,i,iens);
    crm_vvel         (k,j,i,iens) = state_v_wind       (k,j,i,iens);
    crm_wvel         (k,j,i,iens) = state_w_wind       (k,j,i,iens);
    crm_temp         (k,j,i,iens) = state_temperature  (k,j,i,iens);
    // crm_qv           (k,j,i,iens) = state_qv           (k,j,i,iens);
    // crm_qc           (k,j,i,iens) = state_qc           (k,j,i,iens);
    // crm_qr           (k,j,i,iens) = state_qr           (k,j,i,iens);
    // crm_qi           (k,j,i,iens) = state_qi           (k,j,i,iens);
    crm_qv           (k,j,i,iens) = rho_v;
    crm_qc           (k,j,i,iens) = rho_c;
    crm_qr           (k,j,i,iens) = rho_r;
    crm_qi           (k,j,i,iens) = rho_i;
    crm_nc           (k,j,i,iens) = state_nc           (k,j,i,iens);
    crm_nr           (k,j,i,iens) = state_nr           (k,j,i,iens);
    crm_ni           (k,j,i,iens) = state_ni           (k,j,i,iens);
    crm_qm           (k,j,i,iens) = state_qm           (k,j,i,iens);
    crm_bm           (k,j,i,iens) = state_bm           (k,j,i,iens);
    crm_t_prev       (k,j,i,iens) = state_t_prev       (k,j,i,iens);
    crm_q_prev       (k,j,i,iens) = state_q_prev       (k,j,i,iens);
    crm_shoc_wthv_sec(k,j,i,iens) = state_shoc_wthv_sec(k,j,i,iens);
    crm_shoc_tk      (k,j,i,iens) = state_shoc_tk      (k,j,i,iens);
    crm_shoc_tkh     (k,j,i,iens) = state_shoc_tkh     (k,j,i,iens);
    crm_shoc_cldfrac (k,j,i,iens) = state_shoc_cldfrac (k,j,i,iens);
    crm_shoc_relvar  (k,j,i,iens) = state_shoc_relvar  (k,j,i,iens);
  });
  //------------------------------------------------------------------------------------------------
}

// 
inline void pam_state_copy_output_to_gcm( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int nz   = dm_device.get_dimension_size("z"   );
  int ny   = dm_device.get_dimension_size("y"   );
  int nx   = dm_device.get_dimension_size("x"   );
  int nens = dm_device.get_dimension_size("nens");
  //------------------------------------------------------------------------------------------------
  auto crm_rho_d                = dm_device.get<real,4>("density_dry");
  auto crm_uvel                 = dm_device.get<real,4>("uvel");
  auto crm_vvel                 = dm_device.get<real,4>("vvel");
  auto crm_wvel                 = dm_device.get<real,4>("wvel");
  auto crm_temp                 = dm_device.get<real,4>("temp");
  auto crm_rho_v                = dm_device.get<real,4>("water_vapor");
  auto crm_rho_c                = dm_device.get<real,4>("cloud_water");
  auto crm_rho_r                = dm_device.get<real,4>("rain");
  auto crm_rho_i                = dm_device.get<real,4>("ice");
  auto crm_num_c                = dm_device.get<real,4>("cloud_water_num");
  auto crm_num_r                = dm_device.get<real,4>("rain_num");
  auto crm_num_i                = dm_device.get<real,4>("ice_num");
  auto crm_qm                   = dm_device.get<real,4>("ice_rime");
  auto crm_bm                   = dm_device.get<real,4>("ice_rime_vol");
  auto crm_t_prev               = dm_device.get<real,4>("qv_prev");
  auto crm_q_prev               = dm_device.get<real,4>("t_prev");
  auto crm_shoc_wthv_sec        = dm_device.get<real,4>("wthv_sec");
  auto crm_shoc_tk              = dm_device.get<real,4>("tk");
  auto crm_shoc_tkh             = dm_device.get<real,4>("tkh");
  auto crm_shoc_cldfrac         = dm_device.get<real,4>("cldfrac");
  auto crm_shoc_relvar          = dm_device.get<real,4>("relvar");
  //------------------------------------------------------------------------------------------------
  // wrap the host CRM state data in YAKL arrays
  auto host_state_u_wind        = dm_host.get<real,4>("state_u_wind");
  auto host_state_v_wind        = dm_host.get<real,4>("state_v_wind");
  auto host_state_w_wind        = dm_host.get<real,4>("state_w_wind");
  auto host_state_temperature   = dm_host.get<real,4>("state_temperature");
  auto host_state_qv            = dm_host.get<real,4>("state_qv");
  auto host_state_qc            = dm_host.get<real,4>("state_qc");
  auto host_state_qr            = dm_host.get<real,4>("state_qr");
  auto host_state_qi            = dm_host.get<real,4>("state_qi");
  auto host_state_nc            = dm_host.get<real,4>("state_nc");
  auto host_state_nr            = dm_host.get<real,4>("state_nr");
  auto host_state_ni            = dm_host.get<real,4>("state_ni");
  auto host_state_qm            = dm_host.get<real,4>("state_qm");
  auto host_state_bm            = dm_host.get<real,4>("state_bm");
  auto host_state_t_prev        = dm_host.get<real,4>("state_t_prev");
  auto host_state_q_prev        = dm_host.get<real,4>("state_q_prev");
  auto host_state_shoc_wthv_sec = dm_host.get<real,4>("state_shoc_wthv_sec");
  auto host_state_shoc_tk       = dm_host.get<real,4>("state_shoc_tk");
  auto host_state_shoc_tkh      = dm_host.get<real,4>("state_shoc_tkh");
  auto host_state_shoc_cldfrac  = dm_host.get<real,4>("state_shoc_cldfrac");
  auto host_state_shoc_relvar   = dm_host.get<real,4>("state_shoc_relvar");
  //------------------------------------------------------------------------------------------------
  // convert densities back to specific mixing ratios
  real4d qv("qv",nz,ny,nx,nens);
  real4d qc("qc",nz,ny,nx,nens);
  real4d qr("qr",nz,ny,nx,nens);
  real4d qi("qi",nz,ny,nx,nens);
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    // convert specific mixing ratios to density
    qv(k,j,i,iens) = crm_rho_v(k,j,i,iens) / (crm_rho_d(k,j,i,iens)+crm_rho_v(k,j,i,iens));
    qc(k,j,i,iens) = crm_rho_c(k,j,i,iens) / (crm_rho_d(k,j,i,iens)+crm_rho_c(k,j,i,iens));
    qr(k,j,i,iens) = crm_rho_r(k,j,i,iens) / (crm_rho_d(k,j,i,iens)+crm_rho_r(k,j,i,iens));
    qi(k,j,i,iens) = crm_rho_i(k,j,i,iens) / (crm_rho_d(k,j,i,iens)+crm_rho_i(k,j,i,iens));
  });
  //------------------------------------------------------------------------------------------------
  // Copy the CRM state to host arrays
  crm_uvel          .deep_copy_to( host_state_u_wind        );
  crm_vvel          .deep_copy_to( host_state_v_wind        );
  crm_wvel          .deep_copy_to( host_state_w_wind        );
  crm_temp          .deep_copy_to( host_state_temperature   );
  // crm_rho_v         .deep_copy_to( host_state_qv            );
  // crm_rho_c         .deep_copy_to( host_state_qc            );
  // crm_rho_i         .deep_copy_to( host_state_qi            );
  // crm_rho_r         .deep_copy_to( host_state_qr            );
  crm_num_c         .deep_copy_to( host_state_nc            );
  crm_num_r         .deep_copy_to( host_state_nr            );
  crm_num_i         .deep_copy_to( host_state_ni            );
  crm_qm            .deep_copy_to( host_state_qm            );
  crm_bm            .deep_copy_to( host_state_bm            );
  crm_t_prev        .deep_copy_to( host_state_t_prev        );
  crm_q_prev        .deep_copy_to( host_state_q_prev        );
  crm_shoc_wthv_sec .deep_copy_to( host_state_shoc_wthv_sec );
  crm_shoc_tk       .deep_copy_to( host_state_shoc_tk       );
  crm_shoc_tkh      .deep_copy_to( host_state_shoc_tkh      );
  crm_shoc_cldfrac  .deep_copy_to( host_state_shoc_cldfrac  );
  crm_shoc_relvar   .deep_copy_to( host_state_shoc_relvar   );
  //------------------------------------------------------------------------------------------------
}
