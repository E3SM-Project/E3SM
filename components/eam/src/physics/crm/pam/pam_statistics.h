#pragma once

#include "pam_coupler.h"
#include "saturation_adjustment.h"

// These routines are used to encapsulate the aggregation
// of various quantities, such as precipitation


real constexpr cld_threshold = .001; // condensate thresdhold for diagnostic cloud fraction


// register and initialize various quantities for statistical calculations
inline void pam_statistics_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto ny         = coupler.get_option<int>("crm_ny");
  auto nx         = coupler.get_option<int>("crm_nx");
  //------------------------------------------------------------------------------------------------
  // aggregated quantities
  dm_device.register_and_allocate<real>("stat_aggregation_cnt",       "number of aggregated samples",  {nens},{"nens"});
  dm_device.register_and_allocate<real>("precip_liq_aggregated",      "aggregated sfc liq precip rate",{nens},{"nens"});
  dm_device.register_and_allocate<real>("precip_ice_aggregated",      "aggregated sfc ice precip rate",{nens},{"nens"});
  dm_device.register_and_allocate<real>("liqwp_aggregated",           "aggregated liquid water path",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("icewp_aggregated",           "aggregated ice water path",     {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("liq_ice_exchange_aggregated","aggregated liq_ice_exchange",   {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("vap_liq_exchange_aggregated","aggregated vap_liq_exchange",   {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("vap_ice_exchange_aggregated","aggregated vap_ice_exchange",   {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("rho_v_forcing_aggregated",   "aggregated rho_v_forcing",      {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("rho_l_forcing_aggregated",   "aggregated rho_l_forcing",      {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("rho_i_forcing_aggregated",   "aggregated rho_i_forcing",      {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("cldfrac_aggregated",         "aggregated cloud fraction",     {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("clear_rh"       ,            "clear air rel humidity",        {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("clear_rh_cnt"   ,            "clear air count",               {nz,nens},{"z","nens"});
  //------------------------------------------------------------------------------------------------
  // aggregated physics tendencies
  // temporary state variables
  dm_device.register_and_allocate<real>("phys_tend_save_temp",  "saved state for tendency", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("phys_tend_save_qv",    "saved state for tendency", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("phys_tend_save_qc",    "saved state for tendency", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("phys_tend_save_qi",    "saved state for tendency", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  dm_device.register_and_allocate<real>("phys_tend_save_qr",    "saved state for tendency", {nz,ny,nx,nens}, {"z","y","x","nens"} );
  // SGS tendencies
  dm_device.register_and_allocate<real>("phys_tend_sgs_cnt",   "count for aggregated SGS tendency ",  {nens},{"nens"});
  dm_device.register_and_allocate<real>("phys_tend_sgs_temp",  "aggregated temperature tend from SGS",{nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sgs_qv",    "aggregated qv tend from SGS",         {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sgs_qc",    "aggregated qc tend from SGS",         {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sgs_qi",    "aggregated qi tend from SGS",         {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sgs_qr",    "aggregated qr tend from SGS",         {nz,nens},{"z","nens"});
  // micro tendencies
  dm_device.register_and_allocate<real>("phys_tend_micro_cnt", "count for aggregated micro tendency ",  {nens},{"nens"});
  dm_device.register_and_allocate<real>("phys_tend_micro_temp","aggregated temperature tend from micro",{nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_micro_qv",  "aggregated qv tend from microphysics",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_micro_qc",  "aggregated qc tend from microphysics",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_micro_qi",  "aggregated qi tend from microphysics",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_micro_qr",  "aggregated qr tend from microphysics",  {nz,nens},{"z","nens"});
  // dycor tendencies
  dm_device.register_and_allocate<real>("phys_tend_dycor_cnt", "count for aggregated dycor tendency ",  {nens},{"nens"});
  dm_device.register_and_allocate<real>("phys_tend_dycor_temp","aggregated temperature tend from dycor",{nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_dycor_qv",  "aggregated qv tend from dycor",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_dycor_qc",  "aggregated qc tend from dycor",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_dycor_qi",  "aggregated qi tend from dycor",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_dycor_qr",  "aggregated qi tend from dycor",  {nz,nens},{"z","nens"});
  // sponge layer tendencies
  dm_device.register_and_allocate<real>("phys_tend_sponge_cnt", "count for aggregated sponge tendency ",  {nens},{"nens"});
  dm_device.register_and_allocate<real>("phys_tend_sponge_temp","aggregated temperature tend from sponge",{nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sponge_qv",  "aggregated qv tend from sponge",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sponge_qc",  "aggregated qc tend from sponge",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sponge_qi",  "aggregated qi tend from sponge",  {nz,nens},{"z","nens"});
  dm_device.register_and_allocate<real>("phys_tend_sponge_qr",  "aggregated qi tend from sponge",  {nz,nens},{"z","nens"});
  //------------------------------------------------------------------------------------------------
  auto stat_aggregation_cnt        = dm_device.get<real,1>("stat_aggregation_cnt");
  auto precip_liq_aggregated       = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice_aggregated       = dm_device.get<real,1>("precip_ice_aggregated");
  auto liqwp_aggregated            = dm_device.get<real,2>("liqwp_aggregated");
  auto icewp_aggregated            = dm_device.get<real,2>("icewp_aggregated");
  auto liq_ice_exchange_aggregated = dm_device.get<real,2>("liq_ice_exchange_aggregated");
  auto vap_liq_exchange_aggregated = dm_device.get<real,2>("vap_liq_exchange_aggregated");
  auto vap_ice_exchange_aggregated = dm_device.get<real,2>("vap_ice_exchange_aggregated");
  auto rho_v_forcing_aggregated    = dm_device.get<real,2>("rho_v_forcing_aggregated");
  auto rho_l_forcing_aggregated    = dm_device.get<real,2>("rho_l_forcing_aggregated");
  auto rho_i_forcing_aggregated    = dm_device.get<real,2>("rho_i_forcing_aggregated");
  auto cldfrac_aggregated          = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh                    = dm_device.get<real,2>("clear_rh");
  auto clear_rh_cnt                = dm_device.get<real,2>("clear_rh_cnt");

  auto phys_tend_sgs_cnt     = dm_device.get<real,1>("phys_tend_sgs_cnt");
  auto phys_tend_sgs_temp    = dm_device.get<real,2>("phys_tend_sgs_temp");
  auto phys_tend_sgs_qv      = dm_device.get<real,2>("phys_tend_sgs_qv");
  auto phys_tend_sgs_qc      = dm_device.get<real,2>("phys_tend_sgs_qc");
  auto phys_tend_sgs_qi      = dm_device.get<real,2>("phys_tend_sgs_qi");
  auto phys_tend_sgs_qr      = dm_device.get<real,2>("phys_tend_sgs_qr");
  auto phys_tend_micro_cnt   = dm_device.get<real,1>("phys_tend_micro_cnt");
  auto phys_tend_micro_temp  = dm_device.get<real,2>("phys_tend_micro_temp");
  auto phys_tend_micro_qv    = dm_device.get<real,2>("phys_tend_micro_qv");
  auto phys_tend_micro_qc    = dm_device.get<real,2>("phys_tend_micro_qc");
  auto phys_tend_micro_qi    = dm_device.get<real,2>("phys_tend_micro_qi");
  auto phys_tend_micro_qr    = dm_device.get<real,2>("phys_tend_micro_qr");
  auto phys_tend_dycor_cnt   = dm_device.get<real,1>("phys_tend_dycor_cnt");
  auto phys_tend_dycor_temp  = dm_device.get<real,2>("phys_tend_dycor_temp");
  auto phys_tend_dycor_qv    = dm_device.get<real,2>("phys_tend_dycor_qv");
  auto phys_tend_dycor_qc    = dm_device.get<real,2>("phys_tend_dycor_qc");
  auto phys_tend_dycor_qi    = dm_device.get<real,2>("phys_tend_dycor_qi");
  auto phys_tend_dycor_qr    = dm_device.get<real,2>("phys_tend_dycor_qr");
  auto phys_tend_sponge_cnt  = dm_device.get<real,1>("phys_tend_sponge_cnt");
  auto phys_tend_sponge_temp = dm_device.get<real,2>("phys_tend_sponge_temp");
  auto phys_tend_sponge_qv   = dm_device.get<real,2>("phys_tend_sponge_qv");
  auto phys_tend_sponge_qc   = dm_device.get<real,2>("phys_tend_sponge_qc");
  auto phys_tend_sponge_qi   = dm_device.get<real,2>("phys_tend_sponge_qi");
  auto phys_tend_sponge_qr   = dm_device.get<real,2>("phys_tend_sponge_qr");
  // Initialize 0D aggregated quantities
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    stat_aggregation_cnt (iens) = 0;
    phys_tend_sgs_cnt    (iens) = 0;
    phys_tend_micro_cnt  (iens) = 0;
    phys_tend_dycor_cnt  (iens) = 0;
    phys_tend_sponge_cnt (iens) = 0;
    precip_liq_aggregated(iens) = 0;
    precip_ice_aggregated(iens) = 0;
  });
  // Initialize 1D aggregated quantities
  parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k, int iens) {
    liqwp_aggregated           (k,iens) = 0;
    icewp_aggregated           (k,iens) = 0;
    liq_ice_exchange_aggregated(k,iens) = 0;
    vap_liq_exchange_aggregated(k,iens) = 0;
    vap_ice_exchange_aggregated(k,iens) = 0;
    rho_v_forcing_aggregated   (k,iens) = 0;
    rho_l_forcing_aggregated   (k,iens) = 0;
    rho_i_forcing_aggregated   (k,iens) = 0;
    cldfrac_aggregated         (k,iens) = 0;
    clear_rh                   (k,iens) = 0;
    clear_rh_cnt               (k,iens) = 0;

    phys_tend_sgs_temp  (k,iens) = 0;
    phys_tend_sgs_qv    (k,iens) = 0;
    phys_tend_sgs_qc    (k,iens) = 0;
    phys_tend_sgs_qi    (k,iens) = 0;
    phys_tend_sgs_qr    (k,iens) = 0;

    phys_tend_micro_temp(k,iens) = 0;
    phys_tend_micro_qv  (k,iens) = 0;
    phys_tend_micro_qc  (k,iens) = 0;
    phys_tend_micro_qi  (k,iens) = 0;
    phys_tend_micro_qr  (k,iens) = 0;

    phys_tend_dycor_temp (k,iens) = 0;
    phys_tend_dycor_qv   (k,iens) = 0;
    phys_tend_dycor_qc   (k,iens) = 0;
    phys_tend_dycor_qi   (k,iens) = 0;
    phys_tend_dycor_qr   (k,iens) = 0;

    phys_tend_sponge_temp(k,iens) = 0;
    phys_tend_sponge_qv  (k,iens) = 0;
    phys_tend_sponge_qc  (k,iens) = 0;
    phys_tend_sponge_qi  (k,iens) = 0;
    phys_tend_sponge_qr  (k,iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
}

inline void pam_statistics_save_state( pam::PamCoupler &coupler ) {
using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  //------------------------------------------------------------------------------------------------
  // get CRM variables to be aggregated
  auto temp       = dm_device.get<real const,4>("temp"       );
  auto rho_d      = dm_device.get<real const,4>("density_dry");
  auto rho_v      = dm_device.get<real const,4>("water_vapor");
  auto rho_l      = dm_device.get<real const,4>("cloud_water");
  auto rho_i      = dm_device.get<real const,4>("ice"        );
  auto rho_r      = dm_device.get<real const,4>("rain"       );
  //------------------------------------------------------------------------------------------------
  // get temporary saved state variables
  auto phys_tend_save_temp    = dm_device.get<real,4>("phys_tend_save_temp");
  auto phys_tend_save_qv      = dm_device.get<real,4>("phys_tend_save_qv");
  auto phys_tend_save_qc      = dm_device.get<real,4>("phys_tend_save_qc");
  auto phys_tend_save_qi      = dm_device.get<real,4>("phys_tend_save_qi");
  auto phys_tend_save_qr      = dm_device.get<real,4>("phys_tend_save_qr");
  //------------------------------------------------------------------------------------------------
  // save temporary state for physics tendency calculation
  real r_nx_ny  = 1._fp / (nx*ny);
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    phys_tend_save_temp(k,j,i,iens) = temp(k,j,i,iens);
    real rho_total = rho_d(k,j,i,iens) + rho_v(k,j,i,iens);
    phys_tend_save_qv(k,j,i,iens) = rho_v(k,j,i,iens) / rho_total;
    phys_tend_save_qc(k,j,i,iens) = rho_l(k,j,i,iens) / rho_total;
    phys_tend_save_qi(k,j,i,iens) = rho_i(k,j,i,iens) / rho_total;
    phys_tend_save_qr(k,j,i,iens) = rho_r(k,j,i,iens) / rho_total;
  });
  //------------------------------------------------------------------------------------------------
}

inline void pam_statistics_aggregate_tendency( pam::PamCoupler &coupler, std::string scheme ) {
  using yakl::c::SimpleBounds;
  using yakl::atomicAdd;
  auto &dm_device = coupler.get_data_manager_device_readwrite();
  auto &dm_host   = coupler.get_data_manager_host_readwrite();
  auto nens       = coupler.get_option<int>("ncrms");
  auto nz         = coupler.get_option<int>("crm_nz");
  auto nx         = coupler.get_option<int>("crm_nx");
  auto ny         = coupler.get_option<int>("crm_ny");
  auto crm_dt     = coupler.get_option<double>("crm_dt");
  //------------------------------------------------------------------------------------------------
  // get CRM variables to be aggregated
  auto temp       = dm_device.get<real const,4>("temp"       );
  auto rho_d      = dm_device.get<real const,4>("density_dry");
  auto rho_v      = dm_device.get<real const,4>("water_vapor");
  auto rho_l      = dm_device.get<real const,4>("cloud_water");
  auto rho_i      = dm_device.get<real const,4>("ice"        );
  auto rho_r      = dm_device.get<real const,4>("rain"       );
  //------------------------------------------------------------------------------------------------
  // get temporary saved state variables
  auto phys_tend_save_temp    = dm_device.get<real,4>("phys_tend_save_temp");
  auto phys_tend_save_qv      = dm_device.get<real,4>("phys_tend_save_qv");
  auto phys_tend_save_qc      = dm_device.get<real,4>("phys_tend_save_qc");
  auto phys_tend_save_qi      = dm_device.get<real,4>("phys_tend_save_qi");
  auto phys_tend_save_qr      = dm_device.get<real,4>("phys_tend_save_qr");
  //------------------------------------------------------------------------------------------------
  real1d phys_tend_cnt ("phys_tend_cnt" ,nens);
  real2d phys_tend_temp("phys_tend_temp",nz,nens);
  real2d phys_tend_qv  ("phys_tend_qv"  ,nz,nens);
  real2d phys_tend_qc  ("phys_tend_qc"  ,nz,nens);
  real2d phys_tend_qi  ("phys_tend_qi"  ,nz,nens);
  real2d phys_tend_qr  ("phys_tend_qr"  ,nz,nens);
  if (scheme=="sgs"){
    phys_tend_cnt     = dm_device.get<real,1>("phys_tend_sgs_cnt");
    phys_tend_temp    = dm_device.get<real,2>("phys_tend_sgs_temp");
    phys_tend_qv      = dm_device.get<real,2>("phys_tend_sgs_qv");
    phys_tend_qc      = dm_device.get<real,2>("phys_tend_sgs_qc");
    phys_tend_qi      = dm_device.get<real,2>("phys_tend_sgs_qi");
    phys_tend_qr      = dm_device.get<real,2>("phys_tend_sgs_qr");
  }
  if (scheme=="micro"){
    phys_tend_cnt   = dm_device.get<real,1>("phys_tend_micro_cnt");
    phys_tend_temp  = dm_device.get<real,2>("phys_tend_micro_temp");
    phys_tend_qv    = dm_device.get<real,2>("phys_tend_micro_qv");
    phys_tend_qc    = dm_device.get<real,2>("phys_tend_micro_qc");
    phys_tend_qi    = dm_device.get<real,2>("phys_tend_micro_qi");
    phys_tend_qr    = dm_device.get<real,2>("phys_tend_micro_qr");
  }
  if (scheme=="dycor"){
    phys_tend_cnt   = dm_device.get<real,1>("phys_tend_dycor_cnt");
    phys_tend_temp  = dm_device.get<real,2>("phys_tend_dycor_temp");
    phys_tend_qv    = dm_device.get<real,2>("phys_tend_dycor_qv");
    phys_tend_qc    = dm_device.get<real,2>("phys_tend_dycor_qc");
    phys_tend_qi    = dm_device.get<real,2>("phys_tend_dycor_qi");
    phys_tend_qr    = dm_device.get<real,2>("phys_tend_dycor_qr");
  }
  if (scheme=="sponge"){
    phys_tend_cnt   = dm_device.get<real,1>("phys_tend_sponge_cnt");
    phys_tend_temp  = dm_device.get<real,2>("phys_tend_sponge_temp");
    phys_tend_qv    = dm_device.get<real,2>("phys_tend_sponge_qv");
    phys_tend_qc    = dm_device.get<real,2>("phys_tend_sponge_qc");
    phys_tend_qi    = dm_device.get<real,2>("phys_tend_sponge_qi");
    phys_tend_qr    = dm_device.get<real,2>("phys_tend_sponge_qr");
  }
  //------------------------------------------------------------------------------------------------
  // save temporary state for physics tendency calculation
  real r_crm_dt = 1._fp / crm_dt;  // precompute reciprocal to avoid costly divisions
  real r_nx_ny  = 1._fp / (nx*ny);  // precompute reciprocal to avoid costly divisions
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    real rho_total = rho_d(k,j,i,iens) + rho_v(k,j,i,iens);
    real qv_tmp = rho_v(k,j,i,iens) / rho_total;
    real qc_tmp = rho_l(k,j,i,iens) / rho_total;
    real qi_tmp = rho_i(k,j,i,iens) / rho_total;
    real qr_tmp = rho_r(k,j,i,iens) / rho_total;
    real tmp_tend_temp = ( temp(k,j,i,iens) - phys_tend_save_temp(k,j,i,iens) )*r_crm_dt;
    real tmp_tend_qv   = ( qv_tmp           - phys_tend_save_qv  (k,j,i,iens) )*r_crm_dt;
    real tmp_tend_qc   = ( qc_tmp           - phys_tend_save_qc  (k,j,i,iens) )*r_crm_dt;
    real tmp_tend_qi   = ( qi_tmp           - phys_tend_save_qi  (k,j,i,iens) )*r_crm_dt;
    real tmp_tend_qr   = ( qr_tmp           - phys_tend_save_qr  (k,j,i,iens) )*r_crm_dt;
    atomicAdd( phys_tend_temp(k,iens) ,  tmp_tend_temp*r_nx_ny );
    atomicAdd( phys_tend_qv  (k,iens) ,  tmp_tend_qv  *r_nx_ny );
    atomicAdd( phys_tend_qc  (k,iens) ,  tmp_tend_qc  *r_nx_ny );
    atomicAdd( phys_tend_qi  (k,iens) ,  tmp_tend_qi  *r_nx_ny );
    atomicAdd( phys_tend_qr  (k,iens) ,  tmp_tend_qr  *r_nx_ny );
  });
  // update aggregation count for phyics tendencies
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    phys_tend_cnt(iens) += 1;
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
  real grav       = coupler.get_option<real>("grav");
  auto crm_dt     = coupler.get_option<real>("crm_dt");
  auto gcm_nlev   = coupler.get_option<int>("gcm_nlev");
  auto zint       = dm_device.get<real const,2>("vertical_interface_height");
  //------------------------------------------------------------------------------------------------
  // get CRM variables to be aggregated
  auto input_pdel       = dm_host.get<real const,2>("input_pdel").createDeviceCopy();
  auto precip_liq       = dm_device.get<real const,3>("precip_liq_surf_out");
  auto precip_ice       = dm_device.get<real const,3>("precip_ice_surf_out");
  auto temp             = dm_device.get<real const,4>("temp"       );
  auto rho_d            = dm_device.get<real const,4>("density_dry");
  auto rho_v            = dm_device.get<real const,4>("water_vapor");
  auto rho_l            = dm_device.get<real const,4>("cloud_water");
  auto rho_i            = dm_device.get<real const,4>("ice"        );
  auto liq_ice_exchange = dm_device.get<real const,4>("liq_ice_exchange_out");
  auto vap_liq_exchange = dm_device.get<real const,4>("vap_liq_exchange_out");
  auto vap_ice_exchange = dm_device.get<real const,4>("vap_ice_exchange_out");
  auto gcm_forcing_tend_rho_v = dm_device.get<real const,2>("gcm_forcing_tend_rho_v");
  auto gcm_forcing_tend_rho_l = dm_device.get<real const,2>("gcm_forcing_tend_rho_l");
  auto gcm_forcing_tend_rho_i = dm_device.get<real const,2>("gcm_forcing_tend_rho_i");
  //------------------------------------------------------------------------------------------------
  // get aggregation variables
  auto stat_aggregation_cnt        = dm_device.get<real,1>("stat_aggregation_cnt");
  auto precip_liq_aggregated       = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice_aggregated       = dm_device.get<real,1>("precip_ice_aggregated");
  auto liqwp_aggregated            = dm_device.get<real,2>("liqwp_aggregated");
  auto icewp_aggregated            = dm_device.get<real,2>("icewp_aggregated");
  auto liq_ice_exchange_aggregated = dm_device.get<real,2>("liq_ice_exchange_aggregated");
  auto vap_liq_exchange_aggregated = dm_device.get<real,2>("vap_liq_exchange_aggregated");
  auto vap_ice_exchange_aggregated = dm_device.get<real,2>("vap_ice_exchange_aggregated");
  auto rho_v_forcing_aggregated    = dm_device.get<real,2>("rho_v_forcing_aggregated");
  auto rho_l_forcing_aggregated    = dm_device.get<real,2>("rho_l_forcing_aggregated");
  auto rho_i_forcing_aggregated    = dm_device.get<real,2>("rho_i_forcing_aggregated");
  auto cldfrac_aggregated          = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh                    = dm_device.get<real,2>("clear_rh");
  auto clear_rh_cnt                = dm_device.get<real,2>("clear_rh_cnt");
  
  //------------------------------------------------------------------------------------------------
  // update aggregation count
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    stat_aggregation_cnt(iens) = stat_aggregation_cnt(iens) + 1;
  });
  real r_nx_ny  = 1._fp/(nx*ny);
  // aggregate 0D statistics
  parallel_for(SimpleBounds<3>(ny,nx,nens), YAKL_LAMBDA (int j, int i, int iens) {
    // NOTE - precip is already in m/s
    atomicAdd( precip_liq_aggregated(iens), precip_liq(j,i,iens) * r_nx_ny );
    atomicAdd( precip_ice_aggregated(iens), precip_ice(j,i,iens) * r_nx_ny );
  });
  // aggregate 1D statistics
  parallel_for(SimpleBounds<4>(nz,ny,nx,nens), YAKL_LAMBDA (int k, int j, int i, int iens) {
    int k_gcm = gcm_nlev-1-k;
    real rho_total = rho_d(k,j,i,iens) + rho_v(k,j,i,iens);
    atomicAdd( liqwp_aggregated(k,iens), (rho_l(k,j,i,iens)/rho_total) * r_nx_ny * input_pdel(k_gcm,iens)*1000.0/grav );
    atomicAdd( icewp_aggregated(k,iens), (rho_i(k,j,i,iens)/rho_total) * r_nx_ny * input_pdel(k_gcm,iens)*1000.0/grav );
    atomicAdd( liq_ice_exchange_aggregated(k,iens), liq_ice_exchange(k,j,i,iens) * r_nx_ny );
    atomicAdd( vap_liq_exchange_aggregated(k,iens), vap_liq_exchange(k,j,i,iens) * r_nx_ny );
    atomicAdd( vap_ice_exchange_aggregated(k,iens), vap_ice_exchange(k,j,i,iens) * r_nx_ny );
    real rho_sum = rho_l(k,j,i,iens)+rho_i(k,j,i,iens);
    real dz = (zint(k+1,iens) - zint(k,iens));
    if ( ( rho_sum*dz ) > cld_threshold) {
      atomicAdd( cldfrac_aggregated(k,iens), 1.*r_nx_ny );
    } else {
      // calculate RH from vapor pressure and saturation vapor pressure
      real pv = rho_v(k,j,i,iens) * R_v * temp(k,j,i,iens);
      real svp = modules::saturation_vapor_pressure( temp(k,j,i,iens) );
      atomicAdd( clear_rh    (k,iens), pv/svp );
      atomicAdd( clear_rh_cnt(k,iens), 1. );
    }
  });
  parallel_for(SimpleBounds<2>(nz,nens), YAKL_LAMBDA (int k, int iens) {
    rho_v_forcing_aggregated(k,iens) += gcm_forcing_tend_rho_v(k,iens);
    rho_l_forcing_aggregated(k,iens) += gcm_forcing_tend_rho_l(k,iens);
    rho_i_forcing_aggregated(k,iens) += gcm_forcing_tend_rho_i(k,iens);
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
  auto aggregation_cnt       = dm_device.get<real,1>("stat_aggregation_cnt");
  auto precip_liq            = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice            = dm_device.get<real,1>("precip_ice_aggregated");
  auto liqwp                 = dm_device.get<real,2>("liqwp_aggregated");
  auto icewp                 = dm_device.get<real,2>("icewp_aggregated");
  auto liq_ice_exchange      = dm_device.get<real,2>("liq_ice_exchange_aggregated");
  auto vap_liq_exchange      = dm_device.get<real,2>("vap_liq_exchange_aggregated");
  auto vap_ice_exchange      = dm_device.get<real,2>("vap_ice_exchange_aggregated");
  auto rho_v_forcing         = dm_device.get<real,2>("rho_v_forcing_aggregated");
  auto rho_l_forcing         = dm_device.get<real,2>("rho_l_forcing_aggregated");
  auto rho_i_forcing         = dm_device.get<real,2>("rho_i_forcing_aggregated");
  auto cldfrac               = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh              = dm_device.get<real,2>("clear_rh");
  auto clear_rh_cnt          = dm_device.get<real,2>("clear_rh_cnt");

  auto phys_tend_sgs_cnt     = dm_device.get<real,1>("phys_tend_sgs_cnt");
  auto phys_tend_sgs_temp    = dm_device.get<real,2>("phys_tend_sgs_temp");
  auto phys_tend_sgs_qv      = dm_device.get<real,2>("phys_tend_sgs_qv");
  auto phys_tend_sgs_qc      = dm_device.get<real,2>("phys_tend_sgs_qc");
  auto phys_tend_sgs_qi      = dm_device.get<real,2>("phys_tend_sgs_qi");
  auto phys_tend_sgs_qr      = dm_device.get<real,2>("phys_tend_sgs_qr");

  auto phys_tend_micro_cnt   = dm_device.get<real,1>("phys_tend_micro_cnt");
  auto phys_tend_micro_temp  = dm_device.get<real,2>("phys_tend_micro_temp");
  auto phys_tend_micro_qv    = dm_device.get<real,2>("phys_tend_micro_qv");
  auto phys_tend_micro_qc    = dm_device.get<real,2>("phys_tend_micro_qc");
  auto phys_tend_micro_qi    = dm_device.get<real,2>("phys_tend_micro_qi");
  auto phys_tend_micro_qr    = dm_device.get<real,2>("phys_tend_micro_qr");

  auto phys_tend_dycor_cnt   = dm_device.get<real,1>("phys_tend_dycor_cnt");
  auto phys_tend_dycor_temp  = dm_device.get<real,2>("phys_tend_dycor_temp");
  auto phys_tend_dycor_qv    = dm_device.get<real,2>("phys_tend_dycor_qv");
  auto phys_tend_dycor_qc    = dm_device.get<real,2>("phys_tend_dycor_qc");
  auto phys_tend_dycor_qi    = dm_device.get<real,2>("phys_tend_dycor_qi");
  auto phys_tend_dycor_qr    = dm_device.get<real,2>("phys_tend_dycor_qr");

  auto phys_tend_sponge_cnt  = dm_device.get<real,1>("phys_tend_sponge_cnt");
  auto phys_tend_sponge_temp = dm_device.get<real,2>("phys_tend_sponge_temp");
  auto phys_tend_sponge_qv   = dm_device.get<real,2>("phys_tend_sponge_qv");
  auto phys_tend_sponge_qc   = dm_device.get<real,2>("phys_tend_sponge_qc");
  auto phys_tend_sponge_qi   = dm_device.get<real,2>("phys_tend_sponge_qi");
  auto phys_tend_sponge_qr   = dm_device.get<real,2>("phys_tend_sponge_qr");
  // finalize 1D aggregated variables
  parallel_for(SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    precip_liq(iens) = precip_liq(iens) / aggregation_cnt(iens);
    precip_ice(iens) = precip_ice(iens) / aggregation_cnt(iens);
  });
  // finalize 2D aggregated variables
  parallel_for(SimpleBounds<2>(crm_nz,nens), YAKL_LAMBDA (int k, int iens) {
    liqwp           (k,iens) = liqwp           (k,iens) / aggregation_cnt(iens);
    icewp           (k,iens) = icewp           (k,iens) / aggregation_cnt(iens);
    liq_ice_exchange(k,iens) = liq_ice_exchange(k,iens) / aggregation_cnt(iens);
    vap_liq_exchange(k,iens) = vap_liq_exchange(k,iens) / aggregation_cnt(iens);
    vap_ice_exchange(k,iens) = vap_ice_exchange(k,iens) / aggregation_cnt(iens);
    rho_v_forcing   (k,iens) = rho_v_forcing   (k,iens) / aggregation_cnt(iens);
    rho_l_forcing   (k,iens) = rho_l_forcing   (k,iens) / aggregation_cnt(iens);
    rho_i_forcing   (k,iens) = rho_i_forcing   (k,iens) / aggregation_cnt(iens);
    cldfrac         (k,iens) = cldfrac         (k,iens) / aggregation_cnt(iens);
    if (clear_rh_cnt(k,iens)>0) {
      clear_rh(k,iens) = clear_rh(k,iens) / clear_rh_cnt(k,iens);
    }
    phys_tend_sgs_temp  (k,iens) = phys_tend_sgs_temp  (k,iens) / phys_tend_sgs_cnt  (iens);
    phys_tend_sgs_qv    (k,iens) = phys_tend_sgs_qv    (k,iens) / phys_tend_sgs_cnt  (iens);
    phys_tend_sgs_qc    (k,iens) = phys_tend_sgs_qc    (k,iens) / phys_tend_sgs_cnt  (iens);
    phys_tend_sgs_qi    (k,iens) = phys_tend_sgs_qi    (k,iens) / phys_tend_sgs_cnt  (iens);
    phys_tend_sgs_qr    (k,iens) = phys_tend_sgs_qr    (k,iens) / phys_tend_sgs_cnt  (iens);

    phys_tend_micro_temp(k,iens) = phys_tend_micro_temp(k,iens) / phys_tend_micro_cnt(iens);
    phys_tend_micro_qv  (k,iens) = phys_tend_micro_qv  (k,iens) / phys_tend_micro_cnt(iens);
    phys_tend_micro_qc  (k,iens) = phys_tend_micro_qc  (k,iens) / phys_tend_micro_cnt(iens);
    phys_tend_micro_qi  (k,iens) = phys_tend_micro_qi  (k,iens) / phys_tend_micro_cnt(iens);
    phys_tend_micro_qr  (k,iens) = phys_tend_micro_qr  (k,iens) / phys_tend_micro_cnt(iens);

    phys_tend_dycor_temp (k,iens) = phys_tend_dycor_temp (k,iens) / phys_tend_dycor_cnt(iens);
    phys_tend_dycor_qv   (k,iens) = phys_tend_dycor_qv   (k,iens) / phys_tend_dycor_cnt(iens);
    phys_tend_dycor_qc   (k,iens) = phys_tend_dycor_qc   (k,iens) / phys_tend_dycor_cnt(iens);
    phys_tend_dycor_qi   (k,iens) = phys_tend_dycor_qi   (k,iens) / phys_tend_dycor_cnt(iens);
    phys_tend_dycor_qr   (k,iens) = phys_tend_dycor_qr   (k,iens) / phys_tend_dycor_cnt(iens);

    phys_tend_sponge_temp(k,iens) = phys_tend_sponge_temp(k,iens) / phys_tend_sponge_cnt(iens);
    phys_tend_sponge_qv  (k,iens) = phys_tend_sponge_qv  (k,iens) / phys_tend_sponge_cnt(iens);
    phys_tend_sponge_qc  (k,iens) = phys_tend_sponge_qc  (k,iens) / phys_tend_sponge_cnt(iens);
    phys_tend_sponge_qi  (k,iens) = phys_tend_sponge_qi  (k,iens) / phys_tend_sponge_cnt(iens);
    phys_tend_sponge_qr  (k,iens) = phys_tend_sponge_qr  (k,iens) / phys_tend_sponge_cnt(iens);
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
  auto precip_liq            = dm_device.get<real,1>("precip_liq_aggregated");
  auto precip_ice            = dm_device.get<real,1>("precip_ice_aggregated");
  auto liqwp                 = dm_device.get<real,2>("liqwp_aggregated");
  auto icewp                 = dm_device.get<real,2>("icewp_aggregated");
  auto liq_ice_exchange      = dm_device.get<real,2>("liq_ice_exchange_aggregated");
  auto vap_liq_exchange      = dm_device.get<real,2>("vap_liq_exchange_aggregated");
  auto vap_ice_exchange      = dm_device.get<real,2>("vap_ice_exchange_aggregated");
  auto rho_v_forcing         = dm_device.get<real,2>("rho_v_forcing_aggregated");
  auto rho_l_forcing         = dm_device.get<real,2>("rho_l_forcing_aggregated");
  auto rho_i_forcing         = dm_device.get<real,2>("rho_i_forcing_aggregated");
  auto cldfrac               = dm_device.get<real,2>("cldfrac_aggregated");
  auto clear_rh              = dm_device.get<real,2>("clear_rh");

  auto phys_tend_sgs_cnt     = dm_device.get<real,1>("phys_tend_sgs_cnt");
  auto phys_tend_sgs_temp    = dm_device.get<real,2>("phys_tend_sgs_temp");
  auto phys_tend_sgs_qv      = dm_device.get<real,2>("phys_tend_sgs_qv");
  auto phys_tend_sgs_qc      = dm_device.get<real,2>("phys_tend_sgs_qc");
  auto phys_tend_sgs_qi      = dm_device.get<real,2>("phys_tend_sgs_qi");
  auto phys_tend_sgs_qr      = dm_device.get<real,2>("phys_tend_sgs_qr");

  auto phys_tend_micro_cnt   = dm_device.get<real,1>("phys_tend_micro_cnt");
  auto phys_tend_micro_temp  = dm_device.get<real,2>("phys_tend_micro_temp");
  auto phys_tend_micro_qv    = dm_device.get<real,2>("phys_tend_micro_qv");
  auto phys_tend_micro_qc    = dm_device.get<real,2>("phys_tend_micro_qc");
  auto phys_tend_micro_qi    = dm_device.get<real,2>("phys_tend_micro_qi");
  auto phys_tend_micro_qr    = dm_device.get<real,2>("phys_tend_micro_qr");

  auto phys_tend_dycor_cnt   = dm_device.get<real,1>("phys_tend_dycor_cnt");
  auto phys_tend_dycor_temp  = dm_device.get<real,2>("phys_tend_dycor_temp");
  auto phys_tend_dycor_qv    = dm_device.get<real,2>("phys_tend_dycor_qv");
  auto phys_tend_dycor_qc    = dm_device.get<real,2>("phys_tend_dycor_qc");
  auto phys_tend_dycor_qi    = dm_device.get<real,2>("phys_tend_dycor_qi");
  auto phys_tend_dycor_qr    = dm_device.get<real,2>("phys_tend_dycor_qr");

  auto phys_tend_sponge_cnt  = dm_device.get<real,1>("phys_tend_sponge_cnt");
  auto phys_tend_sponge_temp = dm_device.get<real,2>("phys_tend_sponge_temp");
  auto phys_tend_sponge_qv   = dm_device.get<real,2>("phys_tend_sponge_qv");
  auto phys_tend_sponge_qc   = dm_device.get<real,2>("phys_tend_sponge_qc");
  auto phys_tend_sponge_qi   = dm_device.get<real,2>("phys_tend_sponge_qi");
  auto phys_tend_sponge_qr   = dm_device.get<real,2>("phys_tend_sponge_qr");
  //------------------------------------------------------------------------------------------------
  // calculate output precip variables
  // NOTE: the MMF doesn't need to distinguish between "convective" and "large-scale"
  // so just put all precip into the convective category for now
  real1d precip_tot_c("precip_liq_c",nens); // "convective" surface precipitation
  real1d precip_ice_c("precip_ice_c",nens); // "convective" surface precipitation of ice (snow)
  real1d precip_tot_l("precip_liq_l",nens); // "large-scale" surface precipitation
  real1d precip_ice_l("precip_ice_l",nens); // "large-scale" surface precipitation of ice (snow)
  parallel_for("finalize aggregated variables", SimpleBounds<1>(nens), YAKL_LAMBDA (int iens) {
    precip_tot_c(iens) = precip_liq(iens) + precip_ice(iens);
    precip_ice_c(iens) = precip_ice(iens);
    precip_tot_l(iens) = 0;
    precip_ice_l(iens) = 0;
  });
  //------------------------------------------------------------------------------------------------
  // convert variables to GCM vertical grid
  real2d liqwp_gcm           ("liqwp_gcm",            gcm_nlev,nens);
  real2d icewp_gcm           ("icewp_gcm",            gcm_nlev,nens);
  real2d liq_ice_exchange_gcm("liq_ice_exchange_gcm", gcm_nlev,nens);
  real2d vap_liq_exchange_gcm("vap_liq_exchange_gcm", gcm_nlev,nens);
  real2d vap_ice_exchange_gcm("vap_ice_exchange_gcm", gcm_nlev,nens);
  real2d rho_v_forcing_gcm   ("rho_v_forcing_gcm",    gcm_nlev,nens);
  real2d rho_l_forcing_gcm   ("rho_l_forcing_gcm",    gcm_nlev,nens);
  real2d rho_i_forcing_gcm   ("rho_i_forcing_gcm",    gcm_nlev,nens);
  real2d cldfrac_gcm         ("cldfrac_gcm",          gcm_nlev,nens);

  real2d phys_tend_sgs_temp_gcm  ("phys_tend_sgs_temp_gcm",  gcm_nlev,nens);
  real2d phys_tend_sgs_qv_gcm    ("phys_tend_sgs_qv_gcm",    gcm_nlev,nens);
  real2d phys_tend_sgs_qc_gcm    ("phys_tend_sgs_qc_gcm",    gcm_nlev,nens);
  real2d phys_tend_sgs_qi_gcm    ("phys_tend_sgs_qi_gcm",    gcm_nlev,nens);
  real2d phys_tend_sgs_qr_gcm    ("phys_tend_sgs_qr_gcm",    gcm_nlev,nens);

  real2d phys_tend_micro_temp_gcm("phys_tend_micro_temp_gcm",gcm_nlev,nens);
  real2d phys_tend_micro_qv_gcm  ("phys_tend_micro_qv_gcm",  gcm_nlev,nens);
  real2d phys_tend_micro_qc_gcm  ("phys_tend_micro_qc_gcm",  gcm_nlev,nens);
  real2d phys_tend_micro_qi_gcm  ("phys_tend_micro_qi_gcm",  gcm_nlev,nens);
  real2d phys_tend_micro_qr_gcm  ("phys_tend_micro_qr_gcm",  gcm_nlev,nens);

  real2d phys_tend_dycor_temp_gcm ("phys_tend_dycor_temp_gcm", gcm_nlev,nens);
  real2d phys_tend_dycor_qv_gcm   ("phys_tend_dycor_qv_gcm",   gcm_nlev,nens);
  real2d phys_tend_dycor_qc_gcm   ("phys_tend_dycor_qc_gcm",   gcm_nlev,nens);
  real2d phys_tend_dycor_qi_gcm   ("phys_tend_dycor_qi_gcm",   gcm_nlev,nens);
  real2d phys_tend_dycor_qr_gcm   ("phys_tend_dycor_qr_gcm",   gcm_nlev,nens);

  real2d phys_tend_sponge_temp_gcm("phys_tend_sponge_temp_gcm",gcm_nlev,nens);
  real2d phys_tend_sponge_qv_gcm  ("phys_tend_sponge_qv_gcm",  gcm_nlev,nens);
  real2d phys_tend_sponge_qc_gcm  ("phys_tend_sponge_qc_gcm",  gcm_nlev,nens);
  real2d phys_tend_sponge_qi_gcm  ("phys_tend_sponge_qi_gcm",  gcm_nlev,nens);
  real2d phys_tend_sponge_qr_gcm  ("phys_tend_sponge_qr_gcm",  gcm_nlev,nens);
  // copy variables to equivalent variables on GCM grid
  parallel_for(SimpleBounds<2>(gcm_nlev,nens), YAKL_LAMBDA (int k_gcm, int iens) {
    int k_crm = gcm_nlev-1-k_gcm;
    if (k_crm<crm_nz) {
      liqwp_gcm           (k_gcm,iens) = liqwp           (k_crm,iens);
      icewp_gcm           (k_gcm,iens) = icewp           (k_crm,iens);
      liq_ice_exchange_gcm(k_gcm,iens) = liq_ice_exchange(k_crm,iens);
      vap_liq_exchange_gcm(k_gcm,iens) = vap_liq_exchange(k_crm,iens);
      vap_ice_exchange_gcm(k_gcm,iens) = vap_ice_exchange(k_crm,iens);
      rho_v_forcing_gcm   (k_gcm,iens) = rho_v_forcing   (k_crm,iens);
      rho_l_forcing_gcm   (k_gcm,iens) = rho_l_forcing   (k_crm,iens);
      rho_i_forcing_gcm   (k_gcm,iens) = rho_i_forcing   (k_crm,iens);
      cldfrac_gcm         (k_gcm,iens) = cldfrac         (k_crm,iens);

      phys_tend_sgs_temp_gcm  (k_gcm,iens) = phys_tend_sgs_temp  (k_crm,iens);
      phys_tend_sgs_qv_gcm    (k_gcm,iens) = phys_tend_sgs_qv    (k_crm,iens);
      phys_tend_sgs_qc_gcm    (k_gcm,iens) = phys_tend_sgs_qc    (k_crm,iens);
      phys_tend_sgs_qi_gcm    (k_gcm,iens) = phys_tend_sgs_qi    (k_crm,iens);
      phys_tend_sgs_qr_gcm    (k_gcm,iens) = phys_tend_sgs_qr    (k_crm,iens);

      phys_tend_micro_temp_gcm(k_gcm,iens) = phys_tend_micro_temp(k_crm,iens);
      phys_tend_micro_qv_gcm  (k_gcm,iens) = phys_tend_micro_qv  (k_crm,iens);
      phys_tend_micro_qc_gcm  (k_gcm,iens) = phys_tend_micro_qc  (k_crm,iens);
      phys_tend_micro_qi_gcm  (k_gcm,iens) = phys_tend_micro_qi  (k_crm,iens);
      phys_tend_micro_qr_gcm  (k_gcm,iens) = phys_tend_micro_qr  (k_crm,iens);

      phys_tend_dycor_temp_gcm (k_gcm,iens) = phys_tend_dycor_temp (k_crm,iens);
      phys_tend_dycor_qv_gcm   (k_gcm,iens) = phys_tend_dycor_qv   (k_crm,iens);
      phys_tend_dycor_qc_gcm   (k_gcm,iens) = phys_tend_dycor_qc   (k_crm,iens);
      phys_tend_dycor_qi_gcm   (k_gcm,iens) = phys_tend_dycor_qi   (k_crm,iens);
      phys_tend_dycor_qr_gcm   (k_gcm,iens) = phys_tend_dycor_qr   (k_crm,iens);

      phys_tend_sponge_temp_gcm(k_gcm,iens) = phys_tend_sponge_temp(k_crm,iens);
      phys_tend_sponge_qv_gcm  (k_gcm,iens) = phys_tend_sponge_qv  (k_crm,iens);
      phys_tend_sponge_qc_gcm  (k_gcm,iens) = phys_tend_sponge_qc  (k_crm,iens);
      phys_tend_sponge_qi_gcm  (k_gcm,iens) = phys_tend_sponge_qi  (k_crm,iens);
      phys_tend_sponge_qr_gcm  (k_gcm,iens) = phys_tend_sponge_qr  (k_crm,iens);

    } else {
      liqwp_gcm           (k_gcm,iens) = 0.;
      icewp_gcm           (k_gcm,iens) = 0.;
      liq_ice_exchange_gcm(k_gcm,iens) = 0.;
      vap_liq_exchange_gcm(k_gcm,iens) = 0.;
      vap_ice_exchange_gcm(k_gcm,iens) = 0.;
      rho_v_forcing_gcm   (k_gcm,iens) = 0.;
      rho_l_forcing_gcm   (k_gcm,iens) = 0.;
      rho_i_forcing_gcm   (k_gcm,iens) = 0.;
      cldfrac_gcm         (k_gcm,iens) = 0.;

      phys_tend_sgs_temp_gcm  (k_gcm,iens) = 0.;
      phys_tend_sgs_qv_gcm    (k_gcm,iens) = 0.;
      phys_tend_sgs_qc_gcm    (k_gcm,iens) = 0.;
      phys_tend_sgs_qi_gcm    (k_gcm,iens) = 0.;
      phys_tend_sgs_qr_gcm    (k_gcm,iens) = 0.;

      phys_tend_micro_temp_gcm(k_gcm,iens) = 0.;
      phys_tend_micro_qv_gcm  (k_gcm,iens) = 0.;
      phys_tend_micro_qc_gcm  (k_gcm,iens) = 0.;
      phys_tend_micro_qi_gcm  (k_gcm,iens) = 0.;
      phys_tend_micro_qr_gcm  (k_gcm,iens) = 0.;

      phys_tend_dycor_temp_gcm (k_gcm,iens) = 0.;
      phys_tend_dycor_qv_gcm   (k_gcm,iens) = 0.;
      phys_tend_dycor_qc_gcm   (k_gcm,iens) = 0.;
      phys_tend_dycor_qi_gcm   (k_gcm,iens) = 0.;
      phys_tend_dycor_qr_gcm   (k_gcm,iens) = 0.;

      phys_tend_sponge_temp_gcm(k_gcm,iens) = 0.;
      phys_tend_sponge_qv_gcm  (k_gcm,iens) = 0.;
      phys_tend_sponge_qc_gcm  (k_gcm,iens) = 0.;
      phys_tend_sponge_qi_gcm  (k_gcm,iens) = 0.;
      phys_tend_sponge_qr_gcm  (k_gcm,iens) = 0.;
    }
  });
  //------------------------------------------------------------------------------------------------
  // copy data to host
  auto precip_tot_c_host          = dm_host.get<real,1>("output_precc");
  auto precip_ice_c_host          = dm_host.get<real,1>("output_precsc");
  auto precip_tot_l_host          = dm_host.get<real,1>("output_precl");
  auto precip_ice_l_host          = dm_host.get<real,1>("output_precsl");
  auto liqwp_host                 = dm_host.get<real,2>("output_gliqwp");
  auto icewp_host                 = dm_host.get<real,2>("output_gicewp");
  auto liq_ice_exchange_host      = dm_host.get<real,2>("output_liq_ice_exchange");
  auto vap_liq_exchange_host      = dm_host.get<real,2>("output_vap_liq_exchange");
  auto vap_ice_exchange_host      = dm_host.get<real,2>("output_vap_ice_exchange");
  auto output_rho_v_ls            = dm_host.get<real,2>("output_rho_v_ls");
  auto output_rho_l_ls            = dm_host.get<real,2>("output_rho_l_ls");
  auto output_rho_i_ls            = dm_host.get<real,2>("output_rho_i_ls");
  auto cldfrac_host               = dm_host.get<real,2>("output_cld");
  auto clear_rh_host              = dm_host.get<real,2>("output_clear_rh");

  auto phys_tend_sgs_temp_host    = dm_host.get<real,2>("output_dt_sgs");
  auto phys_tend_sgs_qv_host      = dm_host.get<real,2>("output_dqv_sgs");
  auto phys_tend_sgs_qc_host      = dm_host.get<real,2>("output_dqc_sgs");
  auto phys_tend_sgs_qi_host      = dm_host.get<real,2>("output_dqi_sgs");
  auto phys_tend_sgs_qr_host      = dm_host.get<real,2>("output_dqr_sgs");

  auto phys_tend_micro_temp_host  = dm_host.get<real,2>("output_dt_micro");
  auto phys_tend_micro_qv_host    = dm_host.get<real,2>("output_dqv_micro");
  auto phys_tend_micro_qc_host    = dm_host.get<real,2>("output_dqc_micro");
  auto phys_tend_micro_qi_host    = dm_host.get<real,2>("output_dqi_micro");
  auto phys_tend_micro_qr_host    = dm_host.get<real,2>("output_dqr_micro");

  auto phys_tend_dycor_temp_host  = dm_host.get<real,2>("output_dt_dycor");
  auto phys_tend_dycor_qv_host    = dm_host.get<real,2>("output_dqv_dycor");
  auto phys_tend_dycor_qc_host    = dm_host.get<real,2>("output_dqc_dycor");
  auto phys_tend_dycor_qi_host    = dm_host.get<real,2>("output_dqi_dycor");
  auto phys_tend_dycor_qr_host    = dm_host.get<real,2>("output_dqr_dycor");

  auto phys_tend_sponge_temp_host = dm_host.get<real,2>("output_dt_sponge");
  auto phys_tend_sponge_qv_host   = dm_host.get<real,2>("output_dqv_sponge");
  auto phys_tend_sponge_qc_host   = dm_host.get<real,2>("output_dqc_sponge");
  auto phys_tend_sponge_qi_host   = dm_host.get<real,2>("output_dqi_sponge");
  auto phys_tend_sponge_qr_host   = dm_host.get<real,2>("output_dqr_sponge");

  precip_tot_c            .deep_copy_to(precip_tot_c_host);
  precip_ice_c            .deep_copy_to(precip_ice_c_host);
  precip_tot_l            .deep_copy_to(precip_tot_l_host);
  precip_ice_l            .deep_copy_to(precip_ice_l_host);
  liqwp_gcm               .deep_copy_to(liqwp_host);
  icewp_gcm               .deep_copy_to(icewp_host);
  liq_ice_exchange_gcm    .deep_copy_to(liq_ice_exchange_host);
  vap_liq_exchange_gcm    .deep_copy_to(vap_liq_exchange_host);
  vap_ice_exchange_gcm    .deep_copy_to(vap_ice_exchange_host);
  rho_v_forcing_gcm       .deep_copy_to(output_rho_v_ls);
  rho_l_forcing_gcm       .deep_copy_to(output_rho_l_ls);
  rho_i_forcing_gcm       .deep_copy_to(output_rho_i_ls);
  cldfrac_gcm             .deep_copy_to(cldfrac_host);
  clear_rh                .deep_copy_to(clear_rh_host);

  phys_tend_sgs_temp_gcm  .deep_copy_to(phys_tend_sgs_temp_host);
  phys_tend_sgs_qv_gcm    .deep_copy_to(phys_tend_sgs_qv_host);
  phys_tend_sgs_qc_gcm    .deep_copy_to(phys_tend_sgs_qc_host);
  phys_tend_sgs_qi_gcm    .deep_copy_to(phys_tend_sgs_qi_host);
  phys_tend_sgs_qr_gcm    .deep_copy_to(phys_tend_sgs_qr_host);

  phys_tend_micro_temp_gcm.deep_copy_to(phys_tend_micro_temp_host);
  phys_tend_micro_qv_gcm  .deep_copy_to(phys_tend_micro_qv_host);
  phys_tend_micro_qc_gcm  .deep_copy_to(phys_tend_micro_qc_host);
  phys_tend_micro_qi_gcm  .deep_copy_to(phys_tend_micro_qi_host);
  phys_tend_micro_qr_gcm  .deep_copy_to(phys_tend_micro_qr_host);

  phys_tend_dycor_temp_gcm .deep_copy_to(phys_tend_dycor_temp_host);
  phys_tend_dycor_qv_gcm   .deep_copy_to(phys_tend_dycor_qv_host);
  phys_tend_dycor_qc_gcm   .deep_copy_to(phys_tend_dycor_qc_host);
  phys_tend_dycor_qi_gcm   .deep_copy_to(phys_tend_dycor_qi_host);
  phys_tend_dycor_qr_gcm   .deep_copy_to(phys_tend_dycor_qr_host);
  
  phys_tend_sponge_temp_gcm.deep_copy_to(phys_tend_sponge_temp_host);
  phys_tend_sponge_qv_gcm  .deep_copy_to(phys_tend_sponge_qv_host);
  phys_tend_sponge_qc_gcm  .deep_copy_to(phys_tend_sponge_qc_host);
  phys_tend_sponge_qi_gcm  .deep_copy_to(phys_tend_sponge_qi_host);
  phys_tend_sponge_qr_gcm  .deep_copy_to(phys_tend_sponge_qr_host);
  yakl::fence();
  //------------------------------------------------------------------------------------------------
}

