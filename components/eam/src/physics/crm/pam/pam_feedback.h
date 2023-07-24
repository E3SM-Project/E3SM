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
inline void pam_feedback_compute_tendencies( pam::PamCoupler &coupler , real gcm_dt ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  int crm_nz      = dm_device.get_dimension_size("z"   );
  int crm_ny      = dm_device.get_dimension_size("y"   );
  int crm_nx      = dm_device.get_dimension_size("x"   );
  int nens        = dm_device.get_dimension_size("nens");
  int gcm_nlev    = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  // Get current CRM state
  auto rho_d = dm_device.get<real,4>("density_dry");
  auto uvel  = dm_device.get<real,4>("uvel"       );
  auto vvel  = dm_device.get<real,4>("vvel"       );
  auto temp  = dm_device.get<real,4>("temp"       );
  auto rho_v = dm_device.get<real,4>("water_vapor");
  auto rho_c = dm_device.get<real,4>("cloud_water");
  auto rho_i = dm_device.get<real,4>("ice"        );
  //------------------------------------------------------------------------------------------------
  // Get input GCM state
  auto gcm_ul = dm_host.get<real const,2>("input_ul").createDeviceCopy();
  auto gcm_vl = dm_host.get<real const,2>("input_vl").createDeviceCopy();
  auto gcm_tl = dm_host.get<real const,2>("input_tl").createDeviceCopy();
  auto gcm_qv = dm_host.get<real const,2>("input_ql").createDeviceCopy();
  auto gcm_qc = dm_host.get<real const,2>("input_qccl").createDeviceCopy();
  auto gcm_qi = dm_host.get<real const,2>("input_qiil").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // Create arrays to hold the current column average of the CRM internal columns
  real2d crm_hmean_uvel ("crm_hmean_uvel" ,crm_nz,nens);
  real2d crm_hmean_vvel ("crm_hmean_vvel" ,crm_nz,nens);
  real2d crm_hmean_temp ("crm_hmean_temp" ,crm_nz,nens);
  real2d crm_hmean_qv   ("crm_hmean_qv"   ,crm_nz,nens);
  real2d crm_hmean_qc   ("crm_hmean_qc"   ,crm_nz,nens);
  real2d crm_hmean_qi   ("crm_hmean_qi"   ,crm_nz,nens);
  //------------------------------------------------------------------------------------------------
  // We will be essentially reducing a summation to these variables, so initialize them to zero
  parallel_for("Initialize horzontal means", SimpleBounds<2>(crm_nz,nens), YAKL_LAMBDA (int k_crm, int iens) {
    crm_hmean_uvel (k_crm,iens) = 0;
    crm_hmean_vvel (k_crm,iens) = 0;
    crm_hmean_temp (k_crm,iens) = 0;
    crm_hmean_qv   (k_crm,iens) = 0;
    crm_hmean_qc   (k_crm,iens) = 0;
    crm_hmean_qi   (k_crm,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
  // Compute horizontal means
  real r_nx_ny  = 1._fp/(crm_nx*crm_ny);  // precompute reciprocal to avoid costly divisions
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(crm_nz,crm_ny,crm_nx,nens), YAKL_LAMBDA (int k_crm, int j, int i, int iens) {
    real rho_total = rho_d(k_crm,j,i,iens) + rho_v(k_crm,j,i,iens);
    atomicAdd( crm_hmean_uvel(k_crm,iens), uvel (k_crm,j,i,iens)            * r_nx_ny );
    atomicAdd( crm_hmean_vvel(k_crm,iens), vvel (k_crm,j,i,iens)            * r_nx_ny );
    atomicAdd( crm_hmean_temp(k_crm,iens), temp (k_crm,j,i,iens)            * r_nx_ny );
    atomicAdd( crm_hmean_qv  (k_crm,iens),(rho_v(k_crm,j,i,iens)/rho_total) * r_nx_ny );
    atomicAdd( crm_hmean_qc  (k_crm,iens),(rho_c(k_crm,j,i,iens)/rho_total) * r_nx_ny );
    atomicAdd( crm_hmean_qi  (k_crm,iens),(rho_i(k_crm,j,i,iens)/rho_total) * r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
  // Create arrays to hold the feedback tendencies
  dm_device.register_and_allocate<real>("crm_feedback_tend_uvel", "feedback tendency of uvel", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("crm_feedback_tend_vvel", "feedback tendency of vvel", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("crm_feedback_tend_dse" , "feedback tendency of dse",  {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("crm_feedback_tend_qv"  , "feedback tendency of qv",   {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("crm_feedback_tend_qc"  , "feedback tendency of qc",   {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("crm_feedback_tend_qi"  , "feedback tendency of qi",   {gcm_nlev,nens},{"gcm_lev","nens"});
  auto crm_feedback_tend_uvel = dm_device.get<real,2>("crm_feedback_tend_uvel");
  auto crm_feedback_tend_vvel = dm_device.get<real,2>("crm_feedback_tend_vvel");
  auto crm_feedback_tend_dse  = dm_device.get<real,2>("crm_feedback_tend_dse");
  auto crm_feedback_tend_qv   = dm_device.get<real,2>("crm_feedback_tend_qv");
  auto crm_feedback_tend_qc   = dm_device.get<real,2>("crm_feedback_tend_qc");
  auto crm_feedback_tend_qi   = dm_device.get<real,2>("crm_feedback_tend_qi");
  //------------------------------------------------------------------------------------------------
  // Compute feedback tendencies
  real cp_d = coupler.get_option<real>("cp_d");
  real r_gcm_dt = 1._fp / gcm_dt;  // precompute reciprocal to avoid costly divisions
  parallel_for( "Compute CRM feedback tendencies", SimpleBounds<2>(gcm_nlev,nens), YAKL_LAMBDA (int k_gcm, int iens) {
    int k_crm = gcm_nlev-1-k_gcm;
    // if (k_crm<crm_nz-2) { // avoid coupling top 2 layers (things get weird up there)
    if (k_crm<crm_nz) {
      crm_feedback_tend_uvel(k_gcm,iens) = ( crm_hmean_uvel(k_crm,iens) - gcm_ul(k_gcm,iens) )*r_gcm_dt;
      crm_feedback_tend_vvel(k_gcm,iens) = ( crm_hmean_vvel(k_crm,iens) - gcm_vl(k_gcm,iens) )*r_gcm_dt;
      crm_feedback_tend_dse (k_gcm,iens) = ( crm_hmean_temp(k_crm,iens) - gcm_tl(k_gcm,iens) )*r_gcm_dt * cp_d;
      crm_feedback_tend_qv  (k_gcm,iens) = ( crm_hmean_qv  (k_crm,iens) - gcm_qv(k_gcm,iens) )*r_gcm_dt;
      crm_feedback_tend_qc  (k_gcm,iens) = ( crm_hmean_qc  (k_crm,iens) - gcm_qc(k_gcm,iens) )*r_gcm_dt;
      crm_feedback_tend_qi  (k_gcm,iens) = ( crm_hmean_qi  (k_crm,iens) - gcm_qi(k_gcm,iens) )*r_gcm_dt;
    } else {
      crm_feedback_tend_uvel(k_gcm,iens) = 0.;
      crm_feedback_tend_vvel(k_gcm,iens) = 0.;
      crm_feedback_tend_dse (k_gcm,iens) = 0.;
      crm_feedback_tend_qv  (k_gcm,iens) = 0.;
      crm_feedback_tend_qc  (k_gcm,iens) = 0.;
      crm_feedback_tend_qi  (k_gcm,iens) = 0.;
    }
  });
  //------------------------------------------------------------------------------------------------
}


// Compute feedback tendencies for the GCM
inline void pam_feedback_copy_to_host( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  auto crm_feedback_tend_uvel = dm_device.get<real,2>("crm_feedback_tend_uvel");
  auto crm_feedback_tend_vvel = dm_device.get<real,2>("crm_feedback_tend_vvel");
  auto crm_feedback_tend_dse  = dm_device.get<real,2>("crm_feedback_tend_dse");
  auto crm_feedback_tend_qv   = dm_device.get<real,2>("crm_feedback_tend_qv");
  auto crm_feedback_tend_qc   = dm_device.get<real,2>("crm_feedback_tend_qc");
  auto crm_feedback_tend_qi   = dm_device.get<real,2>("crm_feedback_tend_qi");
  //------------------------------------------------------------------------------------------------
  auto output_ultend_host  = dm_host.get<real,2>("output_ultend");
  auto output_vltend_host  = dm_host.get<real,2>("output_vltend");
  auto output_sltend_host  = dm_host.get<real,2>("output_sltend");
  auto output_qvltend_host = dm_host.get<real,2>("output_qltend");
  auto output_qcltend_host = dm_host.get<real,2>("output_qcltend");
  auto output_qiltend_host = dm_host.get<real,2>("output_qiltend");
  //------------------------------------------------------------------------------------------------
  // Copy the data to host
  crm_feedback_tend_uvel.deep_copy_to(output_ultend_host);
  crm_feedback_tend_vvel.deep_copy_to(output_vltend_host);
  crm_feedback_tend_dse .deep_copy_to(output_sltend_host);
  crm_feedback_tend_qv  .deep_copy_to(output_qvltend_host);
  crm_feedback_tend_qc  .deep_copy_to(output_qcltend_host);
  crm_feedback_tend_qi  .deep_copy_to(output_qiltend_host);
  yakl::fence();
  //------------------------------------------------------------------------------------------------
}


