#ifndef SCREAM_P3_FUNCTIONS_F90_HPP
#define SCREAM_P3_FUNCTIONS_F90_HPP

#include "physics/p3/p3_functions.hpp"
#include "physics/share/physics_test_data.hpp"
#include "share/eamxx_types.hpp"
#include "ekat/util/ekat_file_utils.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

namespace scream {
namespace p3 {

///////////////////////////////////////////////////////////////////////////////

/**
 * Structs for holding data related to specific P3 calls; these are used for
 * the BFB unit tests.
 */

struct LookupIceData
{
  // Inputs
  Real qi, ni, qm, rhop;

  // Outputs
  Int  dumi, dumjj, dumii, dumzz;
  Real dum1, dum4, dum5, dum6;

  PTD_RW_SCALARS_ONLY(8, dumi, dumjj, dumii, dumzz, dum1, dum4, dum5, dum6);
};

///////////////////////////////////////////////////////////////////////////////

struct LookupIceDataB
{
  // Inputs
  Real qr, nr;

  // Outputs
  Int dumj;
  Real dum3;

  PTD_RW_SCALARS_ONLY(2, dumj, dum3);
};

///////////////////////////////////////////////////////////////////////////////

struct AccessLookupTableData
{
  // Inputs
  LookupIceData& lid;
  Int index;

  // Outputs
  Real proc;

  PTD_RW_SCALARS_ONLY(1, proc);
};

///////////////////////////////////////////////////////////////////////////////

struct AccessLookupTableCollData
{
  // Inputs
  LookupIceData& lid;
  LookupIceDataB& lidb;
  Int index;

  // Outputs
  Real proc;

  PTD_RW_SCALARS_ONLY(1, proc);
};

///////////////////////////////////////////////////////////////////////////////

struct BackToCellAverageData
{
  // inputs
  Real cld_frac_l, cld_frac_r, cld_frac_i;

  // in/out
  Real qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend, nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr, qcnuc,
       nc_nuceat_tend, qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend, qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend,
       nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend, ni_selfcollect_tend, qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend,
       qc2qi_berg_tend;
  Real ncheti_cnt=0, qcheti_cnt=0, nicnt=0, qicnt=0, ninuc_cnt=0, qinuc_cnt=0;
  bool use_hetfrz_classnuc=false, context=true;

  // This populates all fields with test data within [0,1].
  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(31, qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend, nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr, qcnuc,
                 nc_nuceat_tend, qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend, qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend,
                 nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend, ni_selfcollect_tend, qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend,
                 qc2qi_berg_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct CloudWaterConservationData
{
  // inputs
  Real qc, dt;

  //output
  Real qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend;
  Real qcheti_cnt=0, qicnt=0;
  bool use_hetfrz_classnuc=false, context=true;

  PTD_RW_SCALARS_ONLY(8, qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend);
};

struct RainWaterConservationData
{
  // inputs
  Real qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt;

  //output
  Real qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend;

  PTD_RW_SCALARS_ONLY(3, qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend);
};

struct IceWaterConservationData
{
  //inputs
  Real qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, dt;
  //output
  Real qi2qv_sublim_tend, qi2qr_melt_tend;

  Real qinuc_cnt=0, qcheti_cnt=0, qicnt=0;
  bool use_hetfrz_classnuc=false, context=true;


  PTD_RW_SCALARS_ONLY(2, qi2qv_sublim_tend, qi2qr_melt_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct CalcRimeDensityData
{
  // inputs
  Real T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend;

  // output
  Real vtrmi1, rho_qm_cloud;

  PTD_RW_SCALARS_ONLY(2, vtrmi1, rho_qm_cloud);
};

///////////////////////////////////////////////////////////////////////////////

struct CldliqImmersionFreezingData
{
  // inputs
  Real T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar;

  // output
  Real qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend;

  PTD_RW_SCALARS_ONLY(2, qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct RainImmersionFreezingData
{
  // inputs
  Real T_atm, lamr, mu_r, cdistr, qr_incld;

  // output
  Real qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend;

  PTD_RW_SCALARS_ONLY(2, qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct DropletSelfCollectionData
{
  // inputs
  Real rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend;

  // output
  Real nc_selfcollect_tend;

  PTD_RW_SCALARS_ONLY(1, nc_selfcollect_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct CloudRainAccretionData
{
  // inputs
  Real rho, inv_rho, qc_incld, nc_incld, qr_incld, inv_qc_relvar;

  // output
  Real qc2qr_accret_tend, nc_accret_tend;

  PTD_RW_SCALARS_ONLY(2, qc2qr_accret_tend, nc_accret_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct CloudWaterAutoconversionData
{
  // inputs
  Real rho, qc_incld, nc_incld, inv_qc_relvar;

  // output
  Real qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr;

  PTD_RW_SCALARS_ONLY(3, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr);
};

///////////////////////////////////////////////////////////////////////////////

struct RainSelfCollectionData
{
  //inputs
  Real rho, qr_incld, nr_incld;

  //output
  Real nr_selfcollect_tend;

  PTD_RW_SCALARS_ONLY(1, nr_selfcollect_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct ImposeMaxTotalNiData{
  // inout
  Real ni_local;

  //input
  Real max_total_ni, inv_rho_local;

  PTD_RW_SCALARS_ONLY(2, ni_local, inv_rho_local);
};

///////////////////////////////////////////////////////////////////////////////

struct IceMeltingData
{
  // inputs
  Real rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld;

  // output
  Real qi2qr_melt_tend,ni2nr_melt_tend;

  PTD_RW_SCALARS_ONLY(2, qi2qr_melt_tend, ni2nr_melt_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct SubgridVarianceScalingData
{
  // inputs
  Real relvar,expon;
  // no outputs - is a function
};

///////////////////////////////////////////////////////////////////////////////

struct GetCloudDsd2Data
{
  // Inputs
  Real qc, rho, nc_in;

  // Outputs
  Real nc_out, mu_c, nu, lamc, cdist, cdist1;

  PTD_RW_SCALARS_ONLY(6, nc_out, mu_c, nu, lamc, cdist, cdist1)
};

//////////////////////////////////////////////////////////////////////////

struct GetRainDsd2Data
{
  // Inputs
  Real qr, nr_in;

  // Outputs
  Real nr_out, lamr, mu_r, cdistr, logn0r;

  PTD_RW_SCALARS_ONLY(5, nr_out, lamr, mu_r, cdistr, logn0r);
};

///////////////////////////////////////////////////////////////////////////////

struct CalcUpwindData : public PhysicsTestData
{
  // Inputs
  Int kts, kte, kdir, kbot, k_qxtop, num_arrays;
  Real dt_sub;
  Real* rho, *inv_rho, *inv_dz;
  Real *vs; // num_arrays x nk

  // In/out
  Real *qnx; // num_arrays x nk

  // Outputs
  Real *fluxes; // num_arrays x nk

  CalcUpwindData(Int kts_, Int kte_, Int kdir_, Int kbot_, Int k_qxtop_, Int num_arrays_, Real dt_sub_);

  PTD_STD_DEF(CalcUpwindData, 7, kts, kte, kdir, kbot, k_qxtop, num_arrays, dt_sub);

  Int nk() const { return (kte - kts) + 1; }

  void convert_to_ptr_arr(std::vector<Real*>& mem_space, Real**& fluxes_, Real**& vs_, Real**& qnx_);
};

///////////////////////////////////////////////////////////////////////////////

struct GenSedData : public CalcUpwindData
{
  // Inputs
  Real Co_max;

  // In/out
  Int k_qxbot;
  Real dt_left, prt_accum;

  GenSedData(Int kts_, Int kte_, Int kdir_, Int k_qxtop_, Int k_qxbot_, Int kbot_, Real Co_max_, Real dt_left_,
             Real prt_accum_, Int num_arrays_);

  PTD_DATA_COPY_CTOR(GenSedData, 10);
  PTD_ASSIGN_OP(GenSedData, 11, kts, kte, kdir, kbot, k_qxtop, num_arrays, dt_sub, Co_max, k_qxbot, dt_left, prt_accum);
  PTD_RW();
  PTD_RW_SCALARS(11, kts, kte, kdir, kbot, k_qxtop, num_arrays, dt_sub, Co_max, k_qxbot, dt_left, prt_accum);
};

///////////////////////////////////////////////////////////////////////////////

struct CloudSedData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 13;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real *qc_incld, *rho, *inv_rho, *cld_frac_l, *acn, *inv_dz;
  Real dt, inv_dt;
  bool do_predict_nc;

  // In/out
  Real *qc, *nc, *nc_incld, *mu_c, *lamc, *qc_tend, *nc_tend;
  Real precip_liq_surf;

  CloudSedData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
               Real dt_, Real inv_dt_, bool do_predict_nc_, Real precip_liq_surf_);

  PTD_STD_DEF(CloudSedData, 9, kts, kte, ktop, kbot, kdir, dt, inv_dt, do_predict_nc, precip_liq_surf);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct IceSedData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 15;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real *rho, *inv_rho, *rhofaci, *cld_frac_i, *inv_dz;
  Real dt, inv_dt;

  // In/out
  Real *qi, *qi_incld, *ni, *ni_incld, *qm, *qm_incld, *bm, *bm_incld, *qi_tend, *ni_tend;
  Real precip_ice_surf;

  IceSedData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
             Real dt_, Real inv_dt_, Real precip_ice_surf_);

  PTD_STD_DEF(IceSedData, 8, kts, kte, ktop, kbot, kdir, dt, inv_dt, precip_ice_surf);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct RainSedData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 14;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real *rho, *inv_rho, *rhofacr, *cld_frac_r, *inv_dz, *qr_incld;
  Real dt, inv_dt;

  // In/out
  Real *qr, *nr, *nr_incld, *mu_r, *lamr, *qr_tend, *nr_tend;
  Real *precip_liq_flux; // has special size (nk+1)
  Real precip_liq_surf;

  RainSedData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
              Real dt_, Real inv_dt_, Real precip_liq_surf_);

  PTD_STD_DEF(RainSedData, 8, kts, kte, ktop, kbot, kdir, dt, inv_dt, precip_liq_surf);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct CalcBulkRhoRimeData
{
  // Inputs
  Real qi_tot;

  // In/out
  Real qi_rim, bi_rim;

  // Outputs
  Real rho_rime;

  PTD_RW_SCALARS_ONLY(3, qi_rim, bi_rim, rho_rime);
};

///////////////////////////////////////////////////////////////////////////////

struct HomogeneousFreezingData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 12;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real* T_atm, *inv_exner, *latent_heat_fusion;

  // In/out
  Real* qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *th_atm;

  HomogeneousFreezingData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_);

  PTD_STD_DEF(HomogeneousFreezingData, 5, kts, kte, ktop, kbot, kdir);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct ComputeRainFallVelocityData
{
  // Inputs
  Real qr_incld, rhofacr;

  // In/out
  Real nr_incld;

  // Outputs
  Real mu_r, lamr, V_qr, V_nr;

  PTD_RW_SCALARS_ONLY(5, nr_incld, mu_r, lamr, V_qr, V_nr);
};

///////////////////////////////////////////////////////////////////////////////

struct GetTimeSpacePhysVarsData
{
  //Inputs
  Real T_atm, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i;

  //Outs
  Real mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii;

  PTD_RW_SCALARS_ONLY(9, mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii);
};

///////////////////////////////////////////////////////////////////////////////

struct P3UpdatePrognosticIceData
{
  // Inputs
  Real qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend, qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend,
    ni2nr_melt_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend, qc2qi_berg_tend, inv_exner, latent_heat_sublim, latent_heat_fusion;
  bool do_predict_nc, log_wetgrowth;
  Real dt, nmltratio, rho_qm_cloud;
  Real ncheti_cnt=0, nicnt=0, ninuc_cnt=0, qcheti_cnt=0, qicnt=0, qinuc_cnt=0;

  // In/outs
  Real th_atm, qv, qi, ni, qm, bm, qc, nc, qr, nr;
  bool use_hetfrz_classnuc=false, context=true;

  PTD_RW_SCALARS_ONLY(10, th_atm, qv, qi, ni, qm, bm, qc, nc, qr, nr);
};

///////////////////////////////////////////////////////////////////////////////

struct EvapRainData
{
  // Inputs
  Real qr_incld, qc_incld, nr_incld, qi_incld, cld_frac_l, cld_frac_r, qv, qv_prev,
    qv_sat_l, qv_sat_i, ab, abi, epsr, epsi_tot, t, t_prev, latent_heat_sublim, dqsdt, dt;

  //Outs
  Real qr2qv_evap_tend, nr_evap_tend;

  PTD_RW_SCALARS_ONLY(2, qr2qv_evap_tend, nr_evap_tend);
};

///////////////////////////////////////////////////////////////////////////////

struct P3UpdatePrognosticLiqData
{
  // Inputs
  Real qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend;

  bool do_predict_nc, do_prescribed_CCN;

  Real inv_rho, inv_exner, latent_heat_vapor, dt;

  // In/outs
  Real th_atm, qv, qc, nc, qr, nr;

  PTD_RW_SCALARS_ONLY(6, th_atm, qv, qc, nc, qr, nr);
};

///////////////////////////////////////////////////////////////////////////////

struct IceDepositionSublimationData
{
  //Inputs
  Real qi_incld, ni_incld, T_atm, qv_sat_l, qv_sat_i, epsi, abi, qv, inv_dt;

  //Outs
  Real qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend;

  // This populates all input fields with test data within [0,1].
  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(4, qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend);
};

struct IceCldliqCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, table_val_qc2qi_collect, qi_incld, qc_incld;
  Real ni_incld, nc_incld;

  // Outputs
  Real qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc;

  PTD_RW_SCALARS_ONLY(4, qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc);
};

struct IceRainCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect, qi_incld;
  Real ni_incld, qr_incld;

  // Outputs
  Real qr2qi_collect_tend, nr_collect_tend;

  PTD_RW_SCALARS_ONLY(2, qr2qi_collect_tend, nr_collect_tend);
};

struct IceSelfCollectionData
{
  // Inputs
  Real rho, rhofaci, table_val_ni_self_collect, eii, qm_incld;
  Real qi_incld, ni_incld;

  // Outputs
  Real ni_selfcollect_tend;

  PTD_RW_SCALARS_ONLY(1, ni_selfcollect_tend);
};

struct IceRelaxationData
{
  // Inputs
  Real rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, mu, sc, qi_incld, ni_incld;

  // Outputs
  Real epsi, epsi_tot;

  PTD_RW_SCALARS_ONLY(2, epsi, epsi_tot);
};

struct CalcLiqRelaxationData
{
  // Inputs
  Real rho, f1r, f2r, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld;

  // Outputs
  Real epsr, epsc;

  // This populates all input fields with test data within [0,1].
  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(2, epsr, epsc);
};

struct IceNucleationData
{
  // Inputs
  Real temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt;

  bool do_predict_nc, do_prescribed_CCN;

  // Outputs
  Real qv2qi_nucleat_tend, ni_nucleat_tend;

  PTD_RW_SCALARS_ONLY(2, qv2qi_nucleat_tend, ni_nucleat_tend);
};

struct IceWetGrowthData
{
  // Inputs
  Real rho, temp, pres, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, kap, mu, sc, qv, qc_incld;
  Real qi_incld, ni_incld, qr_incld;

  // In/Outs
  bool log_wetgrowth;

  Real qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend;

  PTD_RW_SCALARS_ONLY(6, log_wetgrowth, qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend);
};

struct CheckValuesData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 2;

  // Inputs
  Int kts, kte;
  Int timestepcount, source_ind;

  bool force_abort;

  Real *qv, *temp, *col_loc;

  CheckValuesData(Int kts_, Int kte_, Int timestepcount_, Int source_ind_, bool force_abort_);

  // deep copy
  PTD_STD_DEF(CheckValuesData, 5, kts, kte, timestepcount, source_ind, force_abort);

  Int nk() const { return (kte - kts) + 1; }
};

struct IncloudMixingData
{
  // Inputs
  Real qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r;

  // Outputs
  Real qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld;

  PTD_RW_SCALARS_ONLY(8, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld);
};

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart1Data : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 40;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  bool do_predict_nc, do_prescribed_CCN;
  Real dt;
  Real* pres, *dpres, *dz, *nc_nuceat_tend, *inv_exner, *exner, *inv_cld_frac_l, *inv_cld_frac_i, *inv_cld_frac_r, *latent_heat_vapor, *latent_heat_sublim,
     *latent_heat_fusion, *nccn_prescribed;

  // In/out
  Real* T_atm, *rho, *inv_rho, *qv_sat_l, *qv_sat_i, *qv_supersat_i, *rhofacr, *rhofaci,
    *acn, *qv, *th_atm, *qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *qc_incld, *qr_incld, *qi_incld,
    *qm_incld, *nc_incld, *nr_incld, *ni_incld, *bm_incld;

  // Output
  bool is_nucleat_possible, is_hydromet_present;

  P3MainPart1Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                  bool do_predict_nc_, bool do_prescribed_CCN_, Real dt_, bool=false, bool=false);

  PTD_STD_DEF(P3MainPart1Data, 10, kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, is_nucleat_possible, is_hydromet_present);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart2Data : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 76;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  bool do_predict_nc, do_prescribed_CCN;
  Real dt, inv_dt;
  Real* pres, *dpres, *dz, *nc_nuceat_tend, *inv_exner, *exner, *inv_cld_frac_l, *inv_cld_frac_i, *inv_cld_frac_r, *ni_activated, *inv_qc_relvar, *cld_frac_i, *cld_frac_l, *cld_frac_r, *qv_prev, *t_prev;

  // In/out
  Real* T_atm, *rho, *inv_rho, *qv_sat_l, *qv_sat_i, *qv_supersat_i, *rhofacr, *rhofaci, *acn,
    *qv, *th_atm, *qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *latent_heat_vapor, *latent_heat_sublim, *latent_heat_fusion, *qc_incld, *qr_incld,
    *qi_incld, *qm_incld, *nc_incld, *nr_incld, *ni_incld, *bm_incld, *mu_c, *nu, *lamc, *cdist, *cdist1,
    *cdistr, *mu_r, *lamr, *logn0r, *qv2qi_depos_tend, *precip_total_tend, *nevapr, *qr_evap_tend, *vap_liq_exchange,
    *vap_ice_exchange, *liq_ice_exchange, 
    *P3_qr2qv_evap, *P3_qi2qv_sublim, *P3_qc2qr_accret, *P3_qc2qr_autoconv, *P3_qv2qi_vapdep, *P3_qc2qi_berg, *P3_qc2qr_ice_shed, *P3_qc2qi_collect, *P3_qr2qi_collect, *P3_qc2qi_hetero_freeze, *P3_qr2qi_immers_freeze, *P3_qi2qr_melt,
    *pratot, *prctot;

  bool is_hydromet_present;

  P3MainPart2Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                  bool do_predict_nc_, bool do_prescribed_CCN, Real dt_, Real=0., bool=false);

  PTD_STD_DEF(P3MainPart2Data, 10, kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, inv_dt, is_hydromet_present);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart3Data : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 47;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  Real* inv_exner, *cld_frac_l, *cld_frac_r, *cld_frac_i;

  // In/out
  Real* rho, *inv_rho, *rhofaci,
    *qv, *th_atm, *qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *latent_heat_vapor, *latent_heat_sublim,
    *mu_c, *nu, *lamc, *mu_r,
    *lamr, *vap_liq_exchange,
    *ze_rain, *ze_ice, *diag_vm_qi, *diag_eff_radius_qi, *diag_diam_qi, *rho_qi,
    *diag_equiv_reflectivity, *diag_eff_radius_qc, *diag_eff_radius_qr;

  P3MainPart3Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_);

  PTD_STD_DEF(P3MainPart3Data, 5, kts, kte, kbot, ktop, kdir);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct P3MainData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 50;
  static constexpr size_t NUM_INPUT_ARRAYS = 24;

  // Inputs
  Int its, ite, kts, kte, it;
  Real* pres, *dz, *nc_nuceat_tend, *nccn_prescribed, *ni_activated, *dpres, *inv_exner, *cld_frac_i, *cld_frac_l, *cld_frac_r, *inv_qc_relvar, *qv_prev, *t_prev;
  Real dt;
  bool do_predict_nc, do_prescribed_CCN;
  bool use_hetfrz_classnuc=false;
  Real* hetfrz_immersion_nucleation_tend, *hetfrz_contact_nucleation_tend, *hetfrz_deposition_nucleation_tend;

  // In/out
  Real* qc, *nc, *qr, *nr, *qi, *qm, *ni, *bm, *qv, *th_atm;

  // Out
  Real *diag_eff_radius_qc, *diag_eff_radius_qi, *diag_eff_radius_qr, *rho_qi, *mu_c, *lamc, *qv2qi_depos_tend, *precip_total_tend, *nevapr,
       *qr_evap_tend, *liq_ice_exchange, *vap_liq_exchange, *vap_ice_exchange,
       *P3_qr2qv_evap, *P3_qi2qv_sublim, *P3_qc2qr_accret, *P3_qc2qr_autoconv, *P3_qv2qi_vapdep, *P3_qc2qi_berg, *P3_qc2qr_ice_shed, *P3_qc2qi_collect, *P3_qr2qi_collect, *P3_qc2qi_hetero_freeze, *P3_qr2qi_immers_freeze, *P3_qi2qr_melt, *P3_qr_sed, *P3_qc_sed, *P3_qi_sed,
       *precip_liq_flux, *precip_ice_flux, *precip_liq_surf, *precip_ice_surf;
  Real elapsed_s;

  P3MainData(Int its_, Int ite_, Int kts_, Int kte_, Int it_, Real dt_, bool do_predict_nc_, bool do_prescribed_CCN_, Real=0.);

  PTD_STD_DEF(P3MainData, 9, its, ite, kts, kte, it, dt, do_predict_nc, do_prescribed_CCN, elapsed_s);
};

struct IceSupersatConservationData {
  // Inputs
  Real cld_frac_i, qv, qv_sat_i, latent_heat_sublim, t_atm, dt, qi2qv_sublim_tend, qr2qv_evap_tend;

  // Inputs/Outputs
  Real qidep, qinuc;

  bool use_hetfrz_classnuc=false, context=true;
  Real qinuc_cnt=0;

  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(2, qidep, qinuc);
};

struct NcConservationData {
  // Inputs
  Real nc, nc_selfcollect_tend, dt;

  // Inputs/Outputs
  Real nc_collect_tend, nc2ni_immers_freeze_tend, nc_accret_tend, nc2nr_autoconv_tend;
  Real ncheti_cnt=0, nicnt=0;
  bool use_hetfrz_classnuc=false, context=true;

  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(4, nc_collect_tend, nc2ni_immers_freeze_tend, nc_accret_tend, nc2nr_autoconv_tend);
};

struct NrConservationData {
  // Inputs
  Real nr, ni2nr_melt_tend, nr_ice_shed_tend, ncshdc, nc2nr_autoconv_tend, dt, nmltratio;

  // Inputs/Outputs
  Real nr_collect_tend, nr2ni_immers_freeze_tend, nr_selfcollect_tend, nr_evap_tend;

  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(4, nr_collect_tend, nr2ni_immers_freeze_tend, nr_selfcollect_tend, nr_evap_tend);
};

struct NiConservationData {
  // Inputs
  Real ni, ni_nucleat_tend, nr2ni_immers_freeze_tend, nc2ni_immers_freeze_tend, dt;
  Real ncheti_cnt=0, nicnt=0, ninuc_cnt=0;
  bool use_hetfrz_classnuc=false, context=true;

  // Inputs/Outputs
  Real ni2nr_melt_tend, ni_sublim_tend, ni_selfcollect_tend;

  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(3, ni2nr_melt_tend, ni_sublim_tend, ni_selfcollect_tend);
};

struct PreventLiqSupersaturationData {
  // Inputs
  Real pres, t_atm, qv, latent_heat_vapor, latent_heat_sublim, dt, qidep, qinuc;

  // Inputs/Outputs
  Real qi2qv_sublim_tend, qr2qv_evap_tend;

  // This populates all fields with test data within [0,1].
  void randomize(std::mt19937_64& engine);

  PTD_RW_SCALARS_ONLY(2, qi2qv_sublim_tend, qr2qv_evap_tend);
};

/**
 * Convenience functions for calling p3 routines from the host with scalar data.
 * These function will pack your data, sync it to device, call the p3 function,
 * then sync back to host and unpack. These are used by the BFB unit tests.
 */

void calc_first_order_upwind_step_host(
  Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub, Real* rho,
  Real* inv_rho, Real* inv_dz, Int num_arrays, Real** fluxes, Real** vs, Real** qnx);

void generalized_sedimentation_host(Int kts, Int kte, Int kdir, Int k_qxtop, Int *k_qxbot, Int kbot, Real Co_max,
                                 Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
                                 Int num_arrays, Real** vs, Real** fluxes, Real** qnx);

void cloud_sedimentation_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* cld_frac_l, Real* acn, Real* inv_dz,
  Real dt, Real inv_dt, bool do_predict_nc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* precip_liq_surf, Real* qc_tend, Real* nc_tend);

void ice_sedimentation_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* cld_frac_i, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qi, Real* qi_incld, Real* ni, Real* qm, Real* qm_incld, Real* bm, Real* bm_incld,
  Real* ni_incld, Real* precip_ice_surf, Real* qi_tend, Real* ni_tend);

void rain_sedimentation_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* cld_frac_r, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* precip_liq_surf, Real* precip_liq_flux, Real* qr_tend, Real* nr_tend);

void homogeneous_freezing_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* T_atm, Real* inv_exner,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* th_atm);

void check_values_host(Real* Qv, Real* temp, Int kstart, Int kend,
                       Int timestepcount, bool force_abort, Int source_ind, Real* col_loc);

void p3_main_part1_host(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool do_predict_nc, bool do_prescribed_CCN,
  Real dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* nccn_prescribed, Real* inv_exner, Real* exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld,
  Real* qm_incld, Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld,
  bool* is_nucleat_possible, bool* is_hydromet_present);

void p3_main_part2_host(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool do_predict_nc, bool do_prescribed_CCN, Real dt, Real inv_dt,
  const Real *hetfrz_immersion_nucleation_tend, const Real *hetfrz_contact_nucleation_tend, const Real *hetfrz_deposition_nucleation_tend,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* inv_exner, Real* exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i, Real* inv_cld_frac_r, Real* ni_activated, Real* inv_qc_relvar, Real* cld_frac_i, Real* cld_frac_l, Real* cld_frac_r, Real* qv_prev, Real* t_prev,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci, Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni,
  Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld, Real* qm_incld, Real* nc_incld, Real* nr_incld,
  Real* ni_incld, Real* bm_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* qv2qi_depos_tend, Real* precip_total_tend,
  Real* nevapr, Real* qr_evap_tend, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, 
  Real* P3_qr2qv_evap, Real* P3_qi2qv_sublim, Real* P3_qc2qr_accret, Real* P3_qc2qr_autoconv, Real* P3_qv2qi_vapdep, Real* P3_qc2qi_berg, Real* P3_qc2qr_ice_shed, Real* P3_qc2qi_collect, Real* P3_qr2qi_collect, Real* P3_qc2qi_hetero_freeze, Real* P3_qr2qi_immers_freeze, Real* P3_qi2qr_melt,
  Real* pratot,
  Real* prctot, bool* is_hydromet_present);

void p3_main_part3_host(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* inv_exner, Real* cld_frac_l, Real* cld_frac_r, Real* cld_frac_i,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vm_qi, Real* diag_eff_radius_qi, Real* diag_diam_qi, Real* rho_qi, Real* diag_equiv_reflectivity, Real* diag_eff_radius_qc, Real* diag_eff_radius_qr);

Int p3_main_host(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th_atm, Real* qv, Real dt,
  Real* qi, Real* qm, Real* ni, Real* bm, Real* pres, Real* dz,
  Real* nc_nuceat_tend, Real* nccn_prescribed, Real* ni_activated, Real* inv_qc_relvar, Int it, Real* precip_liq_surf,
  Real* precip_ice_surf, Int its, Int ite, Int kts, Int kte, Real* diag_eff_radius_qc,
  Real* diag_eff_radius_qi, Real* diag_eff_radius_qr, Real* rho_qi, bool do_predict_nc, bool do_prescribed_CCN, bool use_hetfrz_classnuc, Real* dpres, Real* inv_exner,
  Real* qv2qi_depos_tend, Real* precip_liq_flux, Real* precip_ice_flux, Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i,
  Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange,
  Real* P3_qr2qv_evap, Real* P3_qi2qv_sublim, Real* P3_qc2qr_accret, Real* P3_qc2qr_autoconv, Real* P3_qv2qi_vapdep, Real* P3_qc2qi_berg, Real* P3_qc2qr_ice_shed, Real* P3_qc2qi_collect, Real* P3_qr2qi_collect, Real* P3_qc2qi_hetero_freeze, Real* P3_qr2qi_immers_freeze, Real* P3_qi2qr_melt, Real* P3_qr_sed, Real* P3_qc_sed, Real* P3_qi_sed,
  Real* qv_prev, Real* t_prev);

}  // namespace p3
}  // namespace scream

#endif
