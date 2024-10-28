#pragma once

#include "pam_coupler.h"


inline void pam_variance_transport_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto nx           = coupler.get_option<int>("crm_nx");
  //------------------------------------------------------------------------------------------------
  dm_device.register_and_allocate<real>("vt_temp",      "temperature variance", {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("vt_rhov",      "water vapor variance", {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("vt_uvel",      "u momentum variance",  {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("vt_temp_pert", "temperature perturbation from horz mean", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("vt_rhov_pert", "water vapor perturbation from horz mean", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("vt_uvel_pert", "u momentum perturbation from horz mean",  {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("vt_temp_forcing_tend", "temperature variance forcing tendency", {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("vt_rhov_forcing_tend", "water vapor variance forcing tendency", {nz,nens}, {"z","nens"} );
  dm_device.register_and_allocate<real>("vt_uvel_forcing_tend", "u momentum variance forcing tendency",  {nz,nens}, {"z","nens"} );
  //------------------------------------------------------------------------------------------------
}


inline void pam_variance_transport_diagnose( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto nx           = coupler.get_option<int>("crm_nx");
  //------------------------------------------------------------------------------------------------
  auto temp         = dm_device.get<real const,4>("temp"       );
  auto rhov         = dm_device.get<real const,4>("water_vapor");
  auto rhoc         = dm_device.get<real const,4>("cloud_water");
  auto rhoi         = dm_device.get<real const,4>("ice"        );
  auto uvel         = dm_device.get<real const,4>("uvel"       );
  auto vt_temp      = dm_device.get<real,2>("vt_temp"       );
  auto vt_rhov      = dm_device.get<real,2>("vt_rhov"       );
  auto vt_uvel      = dm_device.get<real,2>("vt_uvel"       );
  auto vt_temp_pert = dm_device.get<real,4>("vt_temp_pert"  );
  auto vt_rhov_pert = dm_device.get<real,4>("vt_rhov_pert"  );
  auto vt_uvel_pert = dm_device.get<real,4>("vt_uvel_pert"  );
  //------------------------------------------------------------------------------------------------
  // local variables
  real2d temp_mean("temp_mean", nz, nens);
  real2d rhov_mean("rhov_mean", nz, nens);
  real2d uvel_mean("uvel_mean", nz, nens);
  //------------------------------------------------------------------------------------------------
  // initialize variance and horizontal mean
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int n) {
    temp_mean(k,n) = 0.0;
    rhov_mean(k,n) = 0.0;
    uvel_mean(k,n) = 0.0;
    vt_temp(k,n)   = 0.0;
    vt_rhov(k,n)   = 0.0;
    vt_uvel(k,n)   = 0.0;
  });
  //------------------------------------------------------------------------------------------------
  // calculate horizontal mean
  real r_nx_ny  = 1._fp/(nx*ny);  // precompute reciprocal to avoid costly divisions
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int n) {
    real rhot = rhov(k,j,i,n) + rhoc(k,j,i,n) + rhoi(k,j,i,n);
    yakl::atomicAdd( temp_mean(k,n), temp(k,j,i,n)*r_nx_ny );
    yakl::atomicAdd( rhov_mean(k,n), rhot         *r_nx_ny );
    yakl::atomicAdd( uvel_mean(k,n), uvel(k,j,i,n)*r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
  // calculate fluctuations from horz mean
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int n) {
    real rhot = rhov(k,j,i,n) + rhoc(k,j,i,n) + rhoi(k,j,i,n);
    vt_temp_pert(k,j,i,n) = temp(k,j,i,n) - temp_mean(k,n);
    vt_rhov_pert(k,j,i,n) = rhot          - rhov_mean(k,n);
    vt_uvel_pert(k,j,i,n) = uvel(k,j,i,n) - uvel_mean(k,n);
  });
  //------------------------------------------------------------------------------------------------
  // calculate variance
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int n) {
    yakl::atomicAdd( vt_temp(k,n), ( vt_temp_pert(k,j,i,n) * vt_temp_pert(k,j,i,n) ) * r_nx_ny );
    yakl::atomicAdd( vt_rhov(k,n), ( vt_rhov_pert(k,j,i,n) * vt_rhov_pert(k,j,i,n) ) * r_nx_ny );
    yakl::atomicAdd( vt_uvel(k,n), ( vt_uvel_pert(k,j,i,n) * vt_uvel_pert(k,j,i,n) ) * r_nx_ny );
  });
  //------------------------------------------------------------------------------------------------
}


inline void pam_variance_transport_compute_forcing( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto &dm_host     = coupler.get_data_manager_host_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto nx           = coupler.get_option<int>("crm_nx");
  auto gcm_nlev     = coupler.get_option<int>("gcm_nlev");
  auto gcm_dt       = coupler.get_option<real>("gcm_dt");
  //------------------------------------------------------------------------------------------------
  // update CRM variance values
  pam_variance_transport_diagnose(coupler);
  //------------------------------------------------------------------------------------------------
  auto vt_temp              = dm_device.get<real,2>("vt_temp"             );
  auto vt_rhov              = dm_device.get<real,2>("vt_rhov"             );
  auto vt_uvel              = dm_device.get<real,2>("vt_uvel"             );
  auto vt_temp_forcing_tend = dm_device.get<real,2>("vt_temp_forcing_tend");
  auto vt_rhov_forcing_tend = dm_device.get<real,2>("vt_rhov_forcing_tend");
  auto vt_uvel_forcing_tend = dm_device.get<real,2>("vt_uvel_forcing_tend");
  auto gcm_vt_temp          = dm_host.get<real const,2>("input_vt_t").createDeviceCopy();
  auto gcm_vt_rhov          = dm_host.get<real const,2>("input_vt_q").createDeviceCopy();
  auto gcm_vt_uvel          = dm_host.get<real const,2>("input_vt_u").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  // calculate variance transport forcing
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k_crm, int n) {
    int k_gcm = gcm_nlev-1-k_crm;
    vt_temp_forcing_tend(k_crm,n) = ( gcm_vt_temp(k_gcm,n) - vt_temp(k_crm,n) ) / gcm_dt ;
    vt_rhov_forcing_tend(k_crm,n) = ( gcm_vt_rhov(k_gcm,n) - vt_rhov(k_crm,n) ) / gcm_dt ;
    vt_uvel_forcing_tend(k_crm,n) = ( gcm_vt_uvel(k_gcm,n) - vt_uvel(k_crm,n) ) / gcm_dt ;
  });
  //------------------------------------------------------------------------------------------------
}

inline void pam_variance_transport_apply_forcing( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto nx           = coupler.get_option<int>("crm_nx");
  auto crm_dt       = coupler.get_option<real>("crm_dt");
  //------------------------------------------------------------------------------------------------
  // update CRM variance values
  pam_variance_transport_diagnose(coupler);
  //------------------------------------------------------------------------------------------------
  // min and max perturbation scaling values are used to limit the 
  // large-scale forcing from variance transport. This is meant to 
  // protect against creating unstable situations, although 
  // problematic scenarios were extremely rare in testing.
  // A scaling limit of +/- 10% was found to be adequate.
  real constexpr pert_scale_min = 1.0 - 0.05;
  real constexpr pert_scale_max = 1.0 + 0.05;
  //------------------------------------------------------------------------------------------------
  auto temp                 = dm_device.get<real,4>("temp"                );
  auto rhov                 = dm_device.get<real,4>("water_vapor"         );
  auto uvel                 = dm_device.get<real,4>("uvel"                );
  auto vt_temp              = dm_device.get<real,2>("vt_temp"             );
  auto vt_rhov              = dm_device.get<real,2>("vt_rhov"             );
  auto vt_uvel              = dm_device.get<real,2>("vt_uvel"             );
  auto vt_temp_pert         = dm_device.get<real,4>("vt_temp_pert"        );
  auto vt_rhov_pert         = dm_device.get<real,4>("vt_rhov_pert"        );
  auto vt_uvel_pert         = dm_device.get<real,4>("vt_uvel_pert"        );
  auto vt_temp_forcing_tend = dm_device.get<real,2>("vt_temp_forcing_tend");
  auto vt_rhov_forcing_tend = dm_device.get<real,2>("vt_rhov_forcing_tend");
  auto vt_uvel_forcing_tend = dm_device.get<real,2>("vt_uvel_forcing_tend");
  //------------------------------------------------------------------------------------------------
  // local variables
  real2d temp_pert_scale("temp_pert_scale", nz, nens);
  real2d rhov_pert_scale("rhov_pert_scale", nz, nens);
  real2d uvel_pert_scale("uvel_pert_scale", nz, nens);
  //------------------------------------------------------------------------------------------------
  // calculate scaling factor for local perturbations
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int n) {
    // initialize scaling factors to 1.0
    temp_pert_scale(k,n) = 1.0;
    rhov_pert_scale(k,n) = 1.0;
    uvel_pert_scale(k,n) = 1.0;
    // calculate variance scaling factor
    real tmp_t_scale = -1.0;
    real tmp_q_scale = -1.0;
    real tmp_u_scale = -1.0;
    if (vt_temp(k,n)>0.0) { tmp_t_scale = 1. + crm_dt * vt_temp_forcing_tend(k,n) / vt_temp(k,n); }
    if (vt_rhov(k,n)>0.0) { tmp_q_scale = 1. + crm_dt * vt_rhov_forcing_tend(k,n) / vt_rhov(k,n); }
    if (vt_uvel(k,n)>0.0) { tmp_u_scale = 1. + crm_dt * vt_uvel_forcing_tend(k,n) / vt_uvel(k,n); }
    if (tmp_t_scale>0.0){ temp_pert_scale(k,n) = sqrt(tmp_t_scale); }
    if (tmp_q_scale>0.0){ rhov_pert_scale(k,n) = sqrt(tmp_q_scale); }
    if (tmp_u_scale>0.0){ uvel_pert_scale(k,n) = sqrt(tmp_u_scale); }
    // enforce minimum scaling
    temp_pert_scale(k,n) = std::max( temp_pert_scale(k,n), pert_scale_min );
    rhov_pert_scale(k,n) = std::max( rhov_pert_scale(k,n), pert_scale_min );
    uvel_pert_scale(k,n) = std::max( uvel_pert_scale(k,n), pert_scale_min );
    // enforce maximum scaling
    temp_pert_scale(k,n) = std::min( temp_pert_scale(k,n), pert_scale_max );
    rhov_pert_scale(k,n) = std::min( rhov_pert_scale(k,n), pert_scale_max );
    uvel_pert_scale(k,n) = std::min( uvel_pert_scale(k,n), pert_scale_max );
  });
  //------------------------------------------------------------------------------------------------
  // apply variance forcing tendency
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int n) {
    real temp_tend_loc = ( temp_pert_scale(k,n) * vt_temp_pert(k,j,i,n) - vt_temp_pert(k,j,i,n) );
    real rhov_tend_loc = ( rhov_pert_scale(k,n) * vt_rhov_pert(k,j,i,n) - vt_rhov_pert(k,j,i,n) );
    real uvel_tend_loc = ( uvel_pert_scale(k,n) * vt_uvel_pert(k,j,i,n) - vt_uvel_pert(k,j,i,n) );
    temp(k,j,i,n) = temp(k,j,i,n) + temp_tend_loc;
    rhov(k,j,i,n) = rhov(k,j,i,n) + rhov_tend_loc;
    uvel(k,j,i,n) = uvel(k,j,i,n) + uvel_tend_loc;
  });
  //------------------------------------------------------------------------------------------------
}

inline void pam_variance_transport_compute_feedback( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device   = coupler.get_data_manager_device_readwrite();
  auto &dm_host     = coupler.get_data_manager_host_readwrite();
  auto nens         = coupler.get_option<int>("ncrms");
  auto nz           = coupler.get_option<int>("crm_nz");
  auto ny           = coupler.get_option<int>("crm_ny");
  auto nx           = coupler.get_option<int>("crm_nx");
  auto gcm_nlev     = coupler.get_option<int>("gcm_nlev");
  auto gcm_dt       = coupler.get_option<real>("gcm_dt");
  //------------------------------------------------------------------------------------------------
  // update CRM variance values
  pam_variance_transport_diagnose(coupler);
  //------------------------------------------------------------------------------------------------
  auto vt_temp      = dm_device.get<real,2>("vt_temp"    );
  auto vt_rhov      = dm_device.get<real,2>("vt_rhov"    );
  auto vt_uvel      = dm_device.get<real,2>("vt_uvel"    );
  auto gcm_vt_temp  = dm_host.get<real const,2>("input_vt_t").createDeviceCopy();
  auto gcm_vt_rhov  = dm_host.get<real const,2>("input_vt_q").createDeviceCopy();
  auto gcm_vt_uvel  = dm_host.get<real const,2>("input_vt_u").createDeviceCopy();
  //------------------------------------------------------------------------------------------------
  dm_device.register_and_allocate<real>("vt_temp_feedback_tend", "feedback tend of temp variance", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("vt_rhov_feedback_tend", "feedback tend of rhov variance", {gcm_nlev,nens},{"gcm_lev","nens"});
  dm_device.register_and_allocate<real>("vt_uvel_feedback_tend", "feedback tend of uvel variance", {gcm_nlev,nens},{"gcm_lev","nens"});
  auto vt_temp_feedback_tend = dm_device.get<real,2>("vt_temp_feedback_tend"  );
  auto vt_rhov_feedback_tend = dm_device.get<real,2>("vt_rhov_feedback_tend"  );
  auto vt_uvel_feedback_tend = dm_device.get<real,2>("vt_uvel_feedback_tend"  );
  //------------------------------------------------------------------------------------------------
  // calculate variance transport forcing
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k_crm, int n) {
    int k_gcm = gcm_nlev-1-k_crm;
    if (k_crm<nz) {
      vt_temp_feedback_tend(k_gcm,n) = ( vt_temp(k_crm,n) - gcm_vt_temp(k_gcm,n) ) / gcm_dt;
      vt_rhov_feedback_tend(k_gcm,n) = ( vt_rhov(k_crm,n) - gcm_vt_rhov(k_gcm,n) ) / gcm_dt;
      vt_uvel_feedback_tend(k_gcm,n) = ( vt_uvel(k_crm,n) - gcm_vt_uvel(k_gcm,n) ) / gcm_dt;
    } else {
      vt_temp_feedback_tend(k_gcm,n) = 0.0;
      vt_rhov_feedback_tend(k_gcm,n) = 0.0;
      vt_uvel_feedback_tend(k_gcm,n) = 0.0;
    }
  });
  //------------------------------------------------------------------------------------------------
}

// Copy feedback tendencies to GCM host
inline void pam_variance_transport_copy_to_host( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  //------------------------------------------------------------------------------------------------
  auto vt_temp_feedback_tend = dm_device.get<real,2>("vt_temp_feedback_tend"  );
  auto vt_rhov_feedback_tend = dm_device.get<real,2>("vt_rhov_feedback_tend"  );
  auto vt_uvel_feedback_tend = dm_device.get<real,2>("vt_uvel_feedback_tend"  );
  //------------------------------------------------------------------------------------------------
  auto output_vt_temp_tend_host = dm_host.get<real,2>("output_t_vt_tend"  );
  auto output_vt_rhov_tend_host = dm_host.get<real,2>("output_q_vt_tend"  );
  auto output_vt_uvel_tend_host = dm_host.get<real,2>("output_u_vt_tend"  );
  //------------------------------------------------------------------------------------------------
  // Copy the data to host
  vt_temp_feedback_tend.deep_copy_to(output_vt_temp_tend_host);
  vt_rhov_feedback_tend.deep_copy_to(output_vt_rhov_tend_host);
  vt_uvel_feedback_tend.deep_copy_to(output_vt_uvel_tend_host);
  //------------------------------------------------------------------------------------------------
}

