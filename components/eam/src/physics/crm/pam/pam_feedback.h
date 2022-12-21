#pragma once

#include "pam_coupler.h"

// These routines are only called once at the end of the CRM call
// to provide the tendencies and fields to couple the CRM and GCM

// Terminology: 
//     - GCM state at the beginning of the CRM call:       state_gcm
//     - instantaneous horizontally-averaged CRM state:    state_crm
// The GCM forces the CRM in the MMF by computing 
//     (state_gcm - state_crm) / gcm_dt
// The forcing is applied to the CRM state at each CRM physics time step. 
// Similarly, the CRM provides a feedback tendency to the GCM by computing
//     (state_crm - state_gcm) / gcm_dt
// using the final state of the CRM at the end of the integration


// Compute feedback tendencies for the GCM
inline void compute_crm_feedback_tendencies( pam::PamCoupler &coupler , real gcm_dt ) {
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
  // Get current CRM state
  auto rho_d = dm_device.get<real,4>("density_dry");
  auto uvel  = dm_device.get<real,4>("uvel"       );
  auto vvel  = dm_device.get<real,4>("vvel"       );
  auto temp  = dm_device.get<real,4>("temp"       );
  auto rho_v = dm_device.get<real,4>("water_vapor");
  auto rho_l = dm_device.get<real,4>("cloud_water");
  auto rho_i = dm_device.get<real,4>("ice"        );
  //------------------------------------------------------------------------------------------------
  // Get input GCM state
  auto gcm_ul   = dm_host.get<real const,2>("input_ul").createDeviceCopy();
  auto gcm_vl   = dm_host.get<real const,2>("input_vl").createDeviceCopy();
  auto gcm_tl   = dm_host.get<real const,2>("input_tl").createDeviceCopy();
  auto gcm_ql   = dm_host.get<real const,2>("input_ql").createDeviceCopy();
  auto gcm_qccl = dm_host.get<real const,2>("input_qccl").createDeviceCopy();
  auto gcm_qiil = dm_host.get<real const,2>("input_qiil").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // Create arrays to hold the current column average of the CRM internal columns
  real2d crm_hmean_rho_d("crm_hmean_rho_d",nz,nens);
  real2d crm_hmean_uvel ("crm_hmean_uvel" ,nz,nens);
  real2d crm_hmean_vvel ("crm_hmean_vvel" ,nz,nens);
  real2d crm_hmean_temp ("crm_hmean_temp" ,nz,nens);
  real2d crm_hmean_rho_v("crm_hmean_rho_v",nz,nens);
  real2d crm_hmean_rho_l("crm_hmean_rho_l",nz,nens);
  real2d crm_hmean_rho_i("crm_hmean_rho_i",nz,nens);
  //------------------------------------------------------------------------------------------------
  // We will be essentially reducing a summation to these variables, so initialize them to zero
  parallel_for("Initialize horzontal means", SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k, int iens) {
    crm_hmean_rho_d(k,iens) = 0;
    crm_hmean_uvel (k,iens) = 0;
    crm_hmean_vvel (k,iens) = 0;
    crm_hmean_temp (k,iens) = 0;
    crm_hmean_rho_v(k,iens) = 0;
    crm_hmean_rho_l(k,iens) = 0;
    crm_hmean_rho_i(k,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
  // Compute horizontal means
  real r_nx_ny  = 1._fp / (nx*ny);  // precompute reciprocal to avoid costly divisions
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    // yakl::atomicAdd ensures only one thread performs an update at a time to avoid data races and wrong answers
    atomicAdd( crm_hmean_rho_d(k,iens), rho_d(k,j,i,iens) * r_nx_ny );
    atomicAdd( crm_hmean_uvel (k,iens), uvel (k,j,i,iens) * r_nx_ny );
    atomicAdd( crm_hmean_vvel (k,iens), vvel (k,j,i,iens) * r_nx_ny );
    atomicAdd( crm_hmean_temp (k,iens), temp (k,j,i,iens) * r_nx_ny );
    atomicAdd( crm_hmean_rho_v(k,iens), rho_v(k,j,i,iens) * r_nx_ny );
    atomicAdd( crm_hmean_rho_l(k,iens), rho_l(k,j,i,iens) * r_nx_ny );
    atomicAdd( crm_hmean_rho_i(k,iens), rho_i(k,j,i,iens) * r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
  // Compute feedback tendencies
  real2d crm_feedback_tend_uvel("",nz,nens);
  real2d crm_feedback_tend_vvel("",nz,nens);
  real2d crm_feedback_tend_dse ("",nz,nens);
  real2d crm_feedback_tend_qv  ("",nz,nens);
  real2d crm_feedback_tend_ql  ("",nz,nens);
  real2d crm_feedback_tend_qi  ("",nz,nens);
  real cp_d = coupler.get_option<real>("cp_d");
  real r_gcm_dt = 1._fp / gcm_dt;  // precompute reciprocal to avoid costly divisions
  parallel_for( "Compute CRM feedback tendencies", SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k, int iens) {
    real crm_dse = crm_hmean_temp (k,iens) * cp_d;
    real crm_qv  = crm_hmean_rho_v(k,iens) / (crm_hmean_rho_d(k,iens) + crm_hmean_rho_v(k,iens));
    real crm_ql  = crm_hmean_rho_l(k,iens) / (crm_hmean_rho_d(k,iens) + crm_hmean_rho_l(k,iens));
    real crm_qi  = crm_hmean_rho_i(k,iens) / (crm_hmean_rho_d(k,iens) + crm_hmean_rho_i(k,iens));
    crm_feedback_tend_uvel(k,iens) = ( crm_hmean_uvel (k,iens) - gcm_ul  (k,iens) ) * r_gcm_dt;
    crm_feedback_tend_vvel(k,iens) = ( crm_hmean_vvel (k,iens) - gcm_vl  (k,iens) ) * r_gcm_dt;
    crm_feedback_tend_dse (k,iens) = ( crm_dse                 - gcm_tl  (k,iens) ) * r_gcm_dt;
    crm_feedback_tend_qv  (k,iens) = ( crm_qv                  - gcm_ql  (k,iens) ) * r_gcm_dt;
    crm_feedback_tend_ql  (k,iens) = ( crm_ql                  - gcm_qccl(k,iens) ) * r_gcm_dt;
    crm_feedback_tend_qi  (k,iens) = ( crm_qi                  - gcm_qiil(k,iens) ) * r_gcm_dt;
  });
  //------------------------------------------------------------------------------------------------
  // Copy the CRM feedback tendencies to host arrays
  auto output_sltend_host  = dm_host.get<real,1>("output_sltend");
  auto output_qltend_host  = dm_host.get<real,1>("output_qltend");
  auto output_qcltend_host = dm_host.get<real,1>("output_qcltend");
  auto output_qiltend_host = dm_host.get<real,1>("output_qiltend");
  auto output_ultend_host  = dm_host.get<real,1>("output_ultend");
  auto output_vltend_host  = dm_host.get<real,1>("output_vltend");
  crm_feedback_tend_uvel.deep_copy_to(output_ultend_host);
  crm_feedback_tend_vvel.deep_copy_to(output_vltend_host);
  crm_feedback_tend_dse .deep_copy_to(output_sltend_host);
  crm_feedback_tend_qv  .deep_copy_to(output_qltend_host);
  crm_feedback_tend_ql  .deep_copy_to(output_qcltend_host);
  crm_feedback_tend_qi  .deep_copy_to(output_qiltend_host);
  //------------------------------------------------------------------------------------------------
}


// Compute horizontal means for feedback tendencies of variables that are not forced
inline void compute_crm_mean_state( pam::PamCoupler &coupler ) {
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
  // Get current CRM state
  auto nc = dm_device.get<real,4>("cloud_water_num");
  auto ni = dm_device.get<real,4>("ice_num");
  auto qr = dm_device.get<real,4>("rain");
  auto nr = dm_device.get<real,4>("rain_num");
  auto qm = dm_device.get<real,4>("ice_rime");
  auto bm = dm_device.get<real,4>("ice_rime_vol");
  //------------------------------------------------------------------------------------------------
  // Create arrays to hold the current column average of the CRM internal columns
  real2d nc_mean("nc_mean",nz,nens);
  real2d ni_mean("ni_mean",nz,nens);
  real2d qr_mean("qr_mean",nz,nens);
  real2d nr_mean("nr_mean",nz,nens);
  real2d qm_mean("qm_mean",nz,nens);
  real2d bm_mean("bm_mean",nz,nens);
  //------------------------------------------------------------------------------------------------
  // We will be essentially reducing a summation to these variables, so initialize them to zero
  parallel_for("Initialize horzontal means", SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k, int iens) {
    nc_mean(k,iens) = 0;
    ni_mean(k,iens) = 0;
    qr_mean(k,iens) = 0;
    nr_mean(k,iens) = 0;
    qm_mean(k,iens) = 0;
    bm_mean(k,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
  // Compute horizontal means
  real r_nx_ny  = 1._fp / (nx*ny);  // precompute reciprocal to avoid costly divisions
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    // yakl::atomicAdd ensures only one thread performs an update at a time to avoid data races and wrong answers
    atomicAdd( nc_mean(k,iens), nc(k,j,i,iens) * r_nx_ny );
    atomicAdd( ni_mean(k,iens), ni(k,j,i,iens) * r_nx_ny );
    atomicAdd( qr_mean(k,iens), qr(k,j,i,iens) * r_nx_ny );
    atomicAdd( nr_mean(k,iens), nr(k,j,i,iens) * r_nx_ny );
    atomicAdd( qm_mean(k,iens), qm(k,j,i,iens) * r_nx_ny );
    atomicAdd( bm_mean(k,iens), bm(k,j,i,iens) * r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
  // Copy the CRM mean data to host arrays
  auto output_nc_mean = dm_host.get<real,1>("output_nc_mean");
  auto output_ni_mean = dm_host.get<real,1>("output_ni_mean");
  auto output_qr_mean = dm_host.get<real,1>("output_qr_mean");
  auto output_nr_mean = dm_host.get<real,1>("output_nr_mean");
  auto output_qm_mean = dm_host.get<real,1>("output_qm_mean");
  auto output_bm_mean = dm_host.get<real,1>("output_bm_mean");
  nc_mean.deep_copy_to(output_nc_mean);
  ni_mean.deep_copy_to(output_ni_mean);
  qr_mean.deep_copy_to(output_qr_mean);
  nr_mean.deep_copy_to(output_nr_mean);
  qm_mean.deep_copy_to(output_qm_mean);
  bm_mean.deep_copy_to(output_bm_mean);
}
