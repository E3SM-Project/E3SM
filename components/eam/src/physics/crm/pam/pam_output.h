#pragma once

#include "pam_coupler.h"

// Compute horizontal means for feedback tendencies of variables that are not forced
inline void pam_output_compute_means( pam::PamCoupler &coupler ) {
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
  auto crm_rho_d = dm_device.get<real,4>("density_dry");
  auto crm_rho_v = dm_device.get<real,4>("water_vapor");
  auto crm_rho_c = dm_device.get<real,4>("cloud_water");
  auto crm_rho_r = dm_device.get<real,4>("rain");
  auto crm_rho_i = dm_device.get<real,4>("ice");
  auto crm_nc    = dm_device.get<real,4>("cloud_water_num");
  auto crm_ni    = dm_device.get<real,4>("ice_num");
  auto crm_nr    = dm_device.get<real,4>("rain_num");
  auto crm_qm    = dm_device.get<real,4>("ice_rime");
  auto crm_bm    = dm_device.get<real,4>("ice_rime_vol");
  //------------------------------------------------------------------------------------------------
  // Create arrays to hold the current column average of the CRM internal columns
  dm_device.register_and_allocate<real>("qv_mean", "domain mean qv", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("qc_mean", "domain mean qc", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("qi_mean", "domain mean qi", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("qr_mean", "domain mean qr", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("nc_mean", "domain mean nc", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("ni_mean", "domain mean ni", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("nr_mean", "domain mean nr", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("qm_mean", "domain mean qm", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("bm_mean", "domain mean bm", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("rho_d_mean", "domain mean rho_d", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("rho_v_mean", "domain mean rho_v", {gcm_nlev,nens},{"gcm_lev","nens"});
  auto qv_mean = dm_device.get<real,2>("qv_mean");
  auto qc_mean = dm_device.get<real,2>("qc_mean");
  auto qi_mean = dm_device.get<real,2>("qi_mean");
  auto qr_mean = dm_device.get<real,2>("qr_mean");
  auto nc_mean = dm_device.get<real,2>("nc_mean");
  auto ni_mean = dm_device.get<real,2>("ni_mean");
  auto nr_mean = dm_device.get<real,2>("nr_mean");
  auto qm_mean = dm_device.get<real,2>("qm_mean");
  auto bm_mean = dm_device.get<real,2>("bm_mean");
  auto rho_d_mean = dm_device.get<real,2>("rho_d_mean");
  auto rho_v_mean = dm_device.get<real,2>("rho_v_mean");
  //------------------------------------------------------------------------------------------------
  // We will be essentially reducing a summation to these variables, so initialize them to zero
  parallel_for("Initialize horzontal means", SimpleBounds<2>(gcm_nlev,nens), YAKL_LAMBDA (int k_gcm, int iens) {
    qv_mean(k_gcm,iens) = 0;
    qc_mean(k_gcm,iens) = 0;
    qi_mean(k_gcm,iens) = 0;
    qr_mean(k_gcm,iens) = 0;
    nc_mean(k_gcm,iens) = 0;
    ni_mean(k_gcm,iens) = 0;
    nr_mean(k_gcm,iens) = 0;
    qm_mean(k_gcm,iens) = 0;
    bm_mean(k_gcm,iens) = 0;
    rho_d_mean(k_gcm,iens) = 0;
    rho_v_mean(k_gcm,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
  // Compute horizontal means
  real r_nx_ny  = 1._fp/(crm_nx*crm_ny);  // precompute reciprocal to avoid costly divisions
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(crm_nz,crm_ny,crm_nx,nens), YAKL_LAMBDA (int k_crm, int j, int i, int iens) {
    int k_gcm = gcm_nlev-1-k_crm;
    real rho_total = crm_rho_d(k_crm,j,i,iens) + crm_rho_v(k_crm,j,i,iens);
    atomicAdd( qv_mean(k_gcm,iens), (crm_rho_v(k_crm,j,i,iens) / rho_total) * r_nx_ny );
    atomicAdd( qc_mean(k_gcm,iens), (crm_rho_c(k_crm,j,i,iens) / rho_total) * r_nx_ny );
    atomicAdd( qi_mean(k_gcm,iens), (crm_rho_i(k_crm,j,i,iens) / rho_total) * r_nx_ny );
    atomicAdd( qr_mean(k_gcm,iens), (crm_rho_r(k_crm,j,i,iens) / rho_total) * r_nx_ny );
    atomicAdd( nc_mean(k_gcm,iens),  crm_nc   (k_crm,j,i,iens)              * r_nx_ny );
    atomicAdd( ni_mean(k_gcm,iens),  crm_ni   (k_crm,j,i,iens)              * r_nx_ny );
    atomicAdd( nr_mean(k_gcm,iens),  crm_nr   (k_crm,j,i,iens)              * r_nx_ny );
    atomicAdd( qm_mean(k_gcm,iens),  crm_qm   (k_crm,j,i,iens)              * r_nx_ny );
    atomicAdd( bm_mean(k_gcm,iens),  crm_bm   (k_crm,j,i,iens)              * r_nx_ny );
    atomicAdd( rho_d_mean(k_gcm,iens), crm_rho_d(k_crm,j,i,iens)              * r_nx_ny );
    atomicAdd( rho_v_mean(k_gcm,iens), crm_rho_v(k_crm,j,i,iens)              * r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
}


inline void pam_output_copy_to_host( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto crm_nz     = coupler.get_option<int>("crm_nz");
  auto crm_nx     = coupler.get_option<int>("crm_nx");
  auto crm_ny     = coupler.get_option<int>("crm_ny");
  auto gcm_nlev   = coupler.get_option<int>("gcm_nlev");
  //------------------------------------------------------------------------------------------------
  auto qv_mean                = dm_device.get<real const,2>("qv_mean");
  auto qc_mean                = dm_device.get<real const,2>("qc_mean");
  auto qi_mean                = dm_device.get<real const,2>("qi_mean");
  auto qr_mean                = dm_device.get<real const,2>("qr_mean");
  auto nc_mean                = dm_device.get<real const,2>("nc_mean");
  auto ni_mean                = dm_device.get<real const,2>("ni_mean");
  auto nr_mean                = dm_device.get<real const,2>("nr_mean");
  auto qm_mean                = dm_device.get<real const,2>("qm_mean");
  auto bm_mean                = dm_device.get<real const,2>("bm_mean");
  auto rho_d_mean = dm_device.get<real const,2>("rho_d_mean");
  auto rho_v_mean = dm_device.get<real const,2>("rho_v_mean");
  auto gcm_forcing_tend_temp  = dm_device.get<real const,2>("gcm_forcing_tend_temp" );
  auto gcm_forcing_tend_rho_d = dm_device.get<real const,2>("gcm_forcing_tend_rho_d");
  auto gcm_forcing_tend_qtot  = dm_device.get<real const,2>("gcm_forcing_tend_qtot" );
  //------------------------------------------------------------------------------------------------
  // calculate quantites needed for forcing of total water mixing ratio
  auto gcm_qv = dm_host.get<real const,2>("input_ql").createDeviceCopy();
  auto gcm_qc = dm_host.get<real const,2>("input_qccl").createDeviceCopy();
  auto gcm_qi = dm_host.get<real const,2>("input_qiil").createDeviceCopy();
  auto dt_gcm = coupler.get_option<real>("gcm_physics_dt");
  real r_dt_gcm = 1._fp / dt_gcm;
  auto rho_d               = dm_device.get<real,4>("density_dry");
  auto state_rho_d         = dm_host.get<real const,4>("state_rho_dry").createDeviceCopy();
  auto state_qv            = dm_host.get<real const,4>("state_qv").createDeviceCopy();
  auto state_qc            = dm_host.get<real const,4>("state_qc").createDeviceCopy();
  auto state_qi            = dm_host.get<real const,4>("state_qi").createDeviceCopy();
  real2d crm_hmean_qt   ("crm_hmean_qt"   ,crm_nz,nens);
  parallel_for("Initialize horzontal means", SimpleBounds<2>(crm_nz,nens), YAKL_LAMBDA (int k_crm, int iens) {
    crm_hmean_qt(k_crm,iens) = 0;
  });
  real r_nx_ny  = 1._fp/(crm_nx*crm_ny);  // precompute reciprocal to avoid costly divisions
  parallel_for("Horz mean of CRM state", SimpleBounds<4>(crm_nz,crm_ny,crm_nx,nens), YAKL_LAMBDA (int k_crm, int j, int i, int iens) {
    real tmp_qt = state_qv(k_crm,j,i,iens) + state_qc(k_crm,j,i,iens) + state_qi(k_crm,j,i,iens);
    atomicAdd( crm_hmean_qt(k_crm,iens), tmp_qt * r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
  // convert variables to GCM vertical grid
  real2d forcing_tend_out_temp ("forcing_tend_out_temp ",gcm_nlev,nens);
  real2d forcing_tend_out_rho_d("forcing_tend_out_rho_d",gcm_nlev,nens);
  real2d forcing_tend_out_qt("forcing_tend_out_qt ",gcm_nlev,nens);
  parallel_for("set output forcing tendencies", SimpleBounds<2>(gcm_nlev,nens), YAKL_LAMBDA (int k_gcm, int iens) {
    int k_crm = gcm_nlev-1-k_gcm;
    if (k_crm<crm_nz) {
      forcing_tend_out_temp (k_gcm,iens) = gcm_forcing_tend_temp (k_crm,iens);
      forcing_tend_out_rho_d(k_gcm,iens) = gcm_forcing_tend_rho_d(k_crm,iens);
      forcing_tend_out_qt   (k_gcm,iens) = gcm_forcing_tend_qtot (k_crm,iens);
    } else {
      forcing_tend_out_temp (k_gcm,iens) = 0.;
      forcing_tend_out_rho_d(k_gcm,iens) = 0.;
      forcing_tend_out_qt   (k_gcm,iens) = 0.;
    }
  });
  //------------------------------------------------------------------------------------------------
  auto output_qv_mean   = dm_host.get<real,2>("output_qv_mean");
  auto output_qc_mean   = dm_host.get<real,2>("output_qc_mean");
  auto output_qi_mean   = dm_host.get<real,2>("output_qi_mean");
  auto output_qr_mean   = dm_host.get<real,2>("output_qr_mean");
  auto output_nc_mean   = dm_host.get<real,2>("output_nc_mean");
  auto output_ni_mean   = dm_host.get<real,2>("output_ni_mean");
  auto output_nr_mean   = dm_host.get<real,2>("output_nr_mean");
  auto output_qm_mean   = dm_host.get<real,2>("output_qm_mean");
  auto output_bm_mean   = dm_host.get<real,2>("output_bm_mean");
  auto output_rho_d_mean= dm_host.get<real,2>("output_rho_d_mean");
  auto output_rho_v_mean= dm_host.get<real,2>("output_rho_v_mean");
  auto output_t_ls      = dm_host.get<real,2>("output_t_ls");
  auto output_rho_d_ls  = dm_host.get<real,2>("output_rho_d_ls");
  auto output_qt_ls     = dm_host.get<real,2>("output_qt_ls");
  //------------------------------------------------------------------------------------------------
  // Copy the data to host
  qv_mean                 .deep_copy_to(output_qv_mean);
  qc_mean                 .deep_copy_to(output_qc_mean);
  qi_mean                 .deep_copy_to(output_qi_mean);
  qr_mean                 .deep_copy_to(output_qr_mean);
  nc_mean                 .deep_copy_to(output_nc_mean);
  ni_mean                 .deep_copy_to(output_ni_mean);
  nr_mean                 .deep_copy_to(output_nr_mean);
  qm_mean                 .deep_copy_to(output_qm_mean);
  bm_mean                 .deep_copy_to(output_bm_mean);
  rho_d_mean              .deep_copy_to(output_rho_d_mean);
  rho_v_mean              .deep_copy_to(output_rho_v_mean);
  forcing_tend_out_temp   .deep_copy_to(output_t_ls);
  forcing_tend_out_rho_d  .deep_copy_to(output_rho_d_ls);
  forcing_tend_out_qt     .deep_copy_to(output_qt_ls);
  yakl::fence();
  //------------------------------------------------------------------------------------------------
}

