#ifndef SCREAM_P3_FUNCTIONS_F90_HPP
#define SCREAM_P3_FUNCTIONS_F90_HPP

#include "physics/p3/p3_functions.hpp"
#include "physics/share/physics_test_data.hpp"
#include "share/scream_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

//
// Bridge functions to call fortran version of p3 functions from C++
//

namespace scream {
namespace p3 {

//
// Singleton for holding the same global data that are maintained in
// micro_p3, but for use in C++. This data is necessary to complete
// the "bridge" when calling C++ from micro_p3.
//
struct P3GlobalForFortran
{
  using P3F = Functions<Real, DefaultDevice>;

  using view_1d_table = typename P3F::view_1d_table;
  using view_2d_table = typename P3F::view_2d_table;
  using view_ice_table = typename P3F::view_ice_table;
  using view_collect_table = typename P3F::view_collect_table;
  using view_dnu_table = typename P3F::view_dnu_table;

  // All kokkos views must be destructed before Kokkos::finalize
  static void deinit();

  static const view_1d_table& mu_r_table_vals()   { return get().m_mu_r_table_vals; }
  static const view_2d_table& vn_table_vals()     { return get().m_vn_table_vals; }
  static const view_2d_table& vm_table_vals()     { return get().m_vm_table_vals; }
  static const view_2d_table& revap_table_vals()  { return get().m_revap_table_vals; }
  static const view_ice_table& ice_table_vals()       { return get().m_ice_table_vals; }
  static const view_collect_table& collect_table_vals() { return get().m_collect_table_vals; }
  static const view_dnu_table& dnu()         { return get().m_dnu; }

  P3GlobalForFortran() = delete;
  ~P3GlobalForFortran() = delete;
  P3GlobalForFortran(const P3GlobalForFortran&) = delete;
  P3GlobalForFortran& operator=(const P3GlobalForFortran&) = delete;

 private:
  struct Views {
    view_1d_table m_mu_r_table_vals;
    view_2d_table m_vn_table_vals, m_vm_table_vals, m_revap_table_vals;
    view_ice_table m_ice_table_vals;
    view_collect_table m_collect_table_vals;
    view_dnu_table m_dnu;
  };

  static const Views& get();
  static std::shared_ptr<Views> s_views;
};

///////////////////////////////////////////////////////////////////////////////

struct P3InitAFortranData
{
  // Must use Host as device, f90 code might not be able to use Device memory
  using P3F = Functions<Real, HostDevice>;
  using P3C = typename P3F::P3C;

  using view_ice_table = typename P3F::KT::template lview<Real[P3C::densize][P3C::rimsize][P3C::isize][P3C::ice_table_size]>;
  using view_collect_table = typename P3F::KT::template lview<Real[P3C::densize][P3C::rimsize][P3C::isize][P3C::rcollsize][P3C::collect_table_size]>;

  // Need to be LayoutLeft to be fortran compatible
  view_ice_table ice_table_vals;
  view_collect_table collect_table_vals;

  P3InitAFortranData() :
    ice_table_vals("P3InitAFortranData::ice_table_vals"),
    collect_table_vals("P3InitAFortranData::collect_table_vals")
  {}
};

///////////////////////////////////////////////////////////////////////////////

struct LookupIceData
{
  // Inputs
  Real qi, ni, qm, rhop;

  // Outputs
  Int  dumi, dumjj, dumii, dumzz;
  Real dum1, dum4, dum5, dum6;
};

///////////////////////////////////////////////////////////////////////////////

struct LookupIceDataB
{
  // Inputs
  Real qr, nr;

  // Outputs
  Int dumj;
  Real dum3;
};

///////////////////////////////////////////////////////////////////////////////

struct AccessLookupTableData
{
  // Inputs
  LookupIceData& lid;
  Int index;

  // Outputs
  Real proc;
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

  // This populates all fields with test data within [0,1].
  void randomize(std::mt19937_64& engine);
};

///////////////////////////////////////////////////////////////////////////////

struct CloudWaterConservationData
{
  // inputs
  Real qc, dt;

  //output
  Real qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend;
};

struct RainWaterConservationData
{
  // inputs
  Real qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt;

  //output
  Real qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend;
};

struct IceWaterConservationData
{
  //inputs
  Real qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, dt;

  //output
  Real qi2qv_sublim_tend, qi2qr_melt_tend;
};

///////////////////////////////////////////////////////////////////////////////

struct CalcRimeDensityData
{
  // inputs
  Real T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend;

  // output
  Real vtrmi1, rho_qm_cloud;
};

///////////////////////////////////////////////////////////////////////////////

struct CldliqImmersionFreezingData
{
  // inputs
  Real T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar;

  // output
  Real qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend;
};

///////////////////////////////////////////////////////////////////////////////

struct RainImmersionFreezingData
{
  // inputs
  Real T_atm, lamr, mu_r, cdistr, qr_incld;

  // output
  Real qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend;
};

///////////////////////////////////////////////////////////////////////////////

struct DropletSelfCollectionData
{
  // inputs
  Real rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend;

  // output
  Real nc_selfcollect_tend;
};

///////////////////////////////////////////////////////////////////////////////

struct CloudRainAccretionData
{
  // inputs
  Real rho, inv_rho, qc_incld, nc_incld, qr_incld, inv_qc_relvar;

  // output
  Real qc2qr_accret_tend, nc_accret_tend;
};

///////////////////////////////////////////////////////////////////////////////

struct CloudWaterAutoconversionData
{
  // inputs
  Real rho;
  Real qc_incld;
  Real nc_incld;
  Real inv_qc_relvar;

  // output
  Real qc2qr_autoconv_tend;
  Real nc2nr_autoconv_tend;
  Real ncautr;
};

///////////////////////////////////////////////////////////////////////////////

struct RainSelfCollectionData
{
  //inputs
  Real rho, qr_incld, nr_incld;

  //output
  Real nr_selfcollect_tend;
};

///////////////////////////////////////////////////////////////////////////////

struct ImposeMaxTotalNiData{
  // inout
  Real ni_local;

  //input
  Real max_total_ni, inv_rho_local;
};

///////////////////////////////////////////////////////////////////////////////

struct IceMeltingData
{
  // inputs
  Real rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld;

  // output
  Real qi2qr_melt_tend,ni2nr_melt_tend;
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
};

//////////////////////////////////////////////////////////////////////////

struct GetRainDsd2Data
{
  // Inputs
  Real qr, nr_in;

  // Outputs
  Real nr_out, lamr, mu_r, cdistr, logn0r;
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
};

///////////////////////////////////////////////////////////////////////////////

struct GetTimeSpacePhysVarsData
{
  //Inputs
  Real T_atm, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i;

  //Outs
  Real mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii;
};

///////////////////////////////////////////////////////////////////////////////

struct P3UpdatePrognosticIceData
{
  // Inputs
  Real qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend, qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend,
    ni2nr_melt_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend, qc2qi_berg_tend, inv_exner, latent_heat_sublim, latent_heat_fusion;
  bool do_predict_nc, log_wetgrowth;
  Real dt, nmltratio, rho_qm_cloud;

  // In/outs
  Real th_atm, qv, qi, ni, qm, bm, qc, nc, qr, nr;
};

///////////////////////////////////////////////////////////////////////////////

struct EvapRainData
{
  // Inputs
  Real qr_incld, qc_incld, nr_incld, qi_incld, cld_frac_l, cld_frac_r, qv, qv_prev,
    qv_sat_l, qv_sat_i, ab, abi, epsr, epsi_tot, t, t_prev, latent_heat_sublim, dqsdt, dt;

  //Outs
  Real qr2qv_evap_tend, nr_evap_tend;
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

};

struct IceCldliqCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, table_val_qc2qi_collect, qi_incld, qc_incld;
  Real ni_incld, nc_incld;

  // Outputs
  Real qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc;

};

struct IceRainCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect, qi_incld;
  Real ni_incld, qr_incld;

  // Outputs
  Real qr2qi_collect_tend, nr_collect_tend;

};

struct IceSelfCollectionData
{
  // Inputs
  Real rho, rhofaci, table_val_ni_self_collect, eii, qm_incld;
  Real qi_incld, ni_incld;

  // Outputs
  Real ni_selfcollect_tend;

};

struct IceRelaxationData
{
  // Inputs
  Real rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, mu, sc, qi_incld, ni_incld;

  // Outputs
  Real epsi, epsi_tot;
};

struct CalcLiqRelaxationData
{
  // Inputs
  Real rho, f1r, f2r, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld;

  // Outputs
  Real epsr, epsc;

  // This populates all input fields with test data within [0,1].
  void randomize(std::mt19937_64& engine);
};

struct IceNucleationData
{
  // Inputs
  Real temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt;

  bool do_predict_nc, do_prescribed_CCN;

  // Outputs
  Real qv2qi_nucleat_tend, ni_nucleat_tend;
};

struct IceWetGrowthData
{
  // Inputs
  Real rho, temp, pres, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, kap, mu, sc, qv, qc_incld;
  Real qi_incld, ni_incld, qr_incld;

  // In/Outs
  bool log_wetgrowth;

  Real qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend;
};

struct LatentHeatData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 3;

  // Inputs
  Int its, ite, kts, kte;

  // Outputs
  Real* v, *s, *f;

  LatentHeatData(Int its_, Int ite_, Int kts_, Int kte_);

  PTD_STD_DEF(LatentHeatData, 4, its, ite, kts, kte);
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
                  bool do_predict_nc_, bool do_prescribed_CCN_, Real dt_);

  PTD_STD_DEF(P3MainPart1Data, 8, kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart2Data : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 64;

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
    *vap_ice_exchange, *liq_ice_exchange, *pratot, *prctot;

  bool is_hydromet_present;

  P3MainPart2Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                  bool do_predict_nc_, bool do_prescribed_CCN, Real dt_);

  PTD_DATA_COPY_CTOR(P3MainPart2Data, 8);
  PTD_ASSIGN_OP(P3MainPart2Data, 10, kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, inv_dt, is_hydromet_present);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart3Data : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 33;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  Real* inv_exner, *cld_frac_l, *cld_frac_r, *cld_frac_i;

  // In/out
  Real* rho, *inv_rho, *rhofaci,
    *qv, *th_atm, *qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *latent_heat_vapor, *latent_heat_sublim,
    *mu_c, *nu, *lamc, *mu_r,
    *lamr, *vap_liq_exchange,
    *ze_rain, *ze_ice, *diag_vm_qi, *diag_eff_radius_qi, *diag_diam_qi, *rho_qi, *diag_equiv_reflectivity, *diag_eff_radius_qc;

  P3MainPart3Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_);

  PTD_STD_DEF(P3MainPart3Data, 5, kts, kte, kbot, ktop, kdir);

  Int nk() const { return (kte - kts) + 1; }
};

///////////////////////////////////////////////////////////////////////////////

struct P3MainData : public PhysicsTestData
{
  static constexpr size_t NUM_ARRAYS = 34;
  static constexpr size_t NUM_INPUT_ARRAYS = 24;

  // Inputs
  Int its, ite, kts, kte, it;
  Real* pres, *dz, *nc_nuceat_tend, *nccn_prescribed, *ni_activated, *dpres, *inv_exner, *cld_frac_i, *cld_frac_l, *cld_frac_r, *inv_qc_relvar, *qv_prev, *t_prev;
  Real dt;
  bool do_predict_nc, do_prescribed_CCN;

  // In/out
  Real* qc, *nc, *qr, *nr, *qi, *qm, *ni, *bm, *qv, *th_atm;

  // Out
  Real *diag_eff_radius_qc, *diag_eff_radius_qi, *rho_qi, *mu_c, *lamc, *qv2qi_depos_tend, *precip_total_tend, *nevapr,
       *qr_evap_tend, *liq_ice_exchange, *vap_liq_exchange, *vap_ice_exchange,
       *precip_liq_flux, *precip_ice_flux, *precip_liq_surf, *precip_ice_surf;
  Real elapsed_s;

  P3MainData(Int its_, Int ite_, Int kts_, Int kte_, Int it_, Real dt_, bool do_predict_nc_, bool do_prescribed_CCN_);

  PTD_STD_DEF(P3MainData, 8, its, ite, kts, kte, it, dt, do_predict_nc, do_prescribed_CCN);
};

struct IceSupersatConservationData {
  // Inputs
  Real cld_frac_i, qv, qv_sat_i, latent_heat_sublim, t_atm, dt, qi2qv_sublim_tend, qr2qv_evap_tend;

  // Inputs/Outputs
  Real qidep, qinuc;

  void randomize(std::mt19937_64& engine);
};

struct NcConservationData {
  // Inputs
  Real nc, nc_selfcollect_tend, dt;

  // Inputs/Outputs
  Real nc_collect_tend, nc2ni_immers_freeze_tend, nc_accret_tend, nc2nr_autoconv_tend;

  void randomize(std::mt19937_64& engine);
};

struct NrConservationData {
  // Inputs
  Real nr, ni2nr_melt_tend, nr_ice_shed_tend, ncshdc, nc2nr_autoconv_tend, dt, nmltratio;

  // Inputs/Outputs
  Real nr_collect_tend, nr2ni_immers_freeze_tend, nr_selfcollect_tend, nr_evap_tend;

  void randomize(std::mt19937_64& engine);
};

struct NiConservationData {
  // Inputs
  Real ni, ni_nucleat_tend, nr2ni_immers_freeze_tend, nc2ni_immers_freeze_tend, dt;

  // Inputs/Outputs
  Real ni2nr_melt_tend, ni_sublim_tend, ni_selfcollect_tend;

  void randomize(std::mt19937_64& engine);
};

struct PreventLiqSupersaturationData {
  // Inputs
  Real pres, t_atm, qv, latent_heat_vapor, latent_heat_sublim, dt, qidep, qinuc;

  // Inputs/Outputs
  Real qi2qv_sublim_tend, qr2qv_evap_tend;

  // This populates all fields with test data within [0,1].
  void randomize(std::mt19937_64& engine);
};

// Glue functions to call fortran from from C++ with the Data struct
void p3_init_a(P3InitAFortranData& d);
void find_lookuptable_indices_1a(LookupIceData& d);
void find_lookuptable_indices_1b(LookupIceDataB& d);
void access_lookup_table(AccessLookupTableData& d);
void access_lookup_table_coll(AccessLookupTableCollData& d);
void back_to_cell_average(BackToCellAverageData& d);
void cloud_water_conservation(CloudWaterConservationData& d);
void rain_water_conservation(RainWaterConservationData& d);
void ice_water_conservation(IceWaterConservationData& d);
void calc_rime_density(CalcRimeDensityData& d);
void cldliq_immersion_freezing(CldliqImmersionFreezingData& d);
void rain_immersion_freezing(RainImmersionFreezingData& d);
void droplet_self_collection(DropletSelfCollectionData& d);
void cloud_rain_accretion(CloudRainAccretionData& d);
void cloud_water_autoconversion(CloudWaterAutoconversionData& d);
void rain_self_collection(RainSelfCollectionData& d);
void impose_max_total_ni(ImposeMaxTotalNiData& d);
void ice_melting(IceMeltingData& d);
Real subgrid_variance_scaling(SubgridVarianceScalingData& d);
void get_cloud_dsd2(GetCloudDsd2Data& d);
void get_rain_dsd2(GetRainDsd2Data& d);
void calc_first_order_upwind_step(CalcUpwindData& d);
void generalized_sedimentation(GenSedData& d);
void cloud_sedimentation(CloudSedData& d);
void ice_sedimentation(IceSedData& d);
void rain_sedimentation(RainSedData& d);
void calc_bulk_rho_rime(CalcBulkRhoRimeData& d);
void homogeneous_freezing(HomogeneousFreezingData& d);
void compute_rain_fall_velocity(ComputeRainFallVelocityData& d);
void get_time_space_phys_variables(GetTimeSpacePhysVarsData& d);
void update_prognostic_ice(P3UpdatePrognosticIceData& d);
void evaporate_rain(EvapRainData& d);
void update_prognostic_liquid(P3UpdatePrognosticLiqData& d);
void ice_deposition_sublimation(IceDepositionSublimationData& d);
void ice_cldliq_collection(IceCldliqCollectionData& d);
void ice_rain_collection(IceRainCollectionData& d);
void ice_self_collection(IceSelfCollectionData& d);
void ice_relaxation_timescale(IceRelaxationData& d);
void calc_liq_relaxation_timescale(CalcLiqRelaxationData& d);
void ice_nucleation(IceNucleationData& d);
void ice_cldliq_wet_growth(IceWetGrowthData& d);
void get_latent_heat(LatentHeatData& d);
void check_values(CheckValuesData& d);
void calculate_incloud_mixingratios(IncloudMixingData& d);
void p3_main_part1(P3MainPart1Data& d);
void p3_main_part2(P3MainPart2Data& d);
void p3_main_part3(P3MainPart3Data& d);
void p3_main(P3MainData& d);

void ice_supersat_conservation(IceSupersatConservationData& d);
void nc_conservation(NcConservationData& d);
void nr_conservation(NrConservationData& d);
void ni_conservation(NiConservationData& d);
void prevent_liq_supersaturation(PreventLiqSupersaturationData& d);
extern "C" { // _f function decls

void calc_first_order_upwind_step_f(
  Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub, Real* rho,
  Real* inv_rho, Real* inv_dz, Int num_arrays, Real** fluxes, Real** vs, Real** qnx);

void generalized_sedimentation_f(Int kts, Int kte, Int kdir, Int k_qxtop, Int *k_qxbot, Int kbot, Real Co_max,
                                 Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
                                 Int num_arrays, Real** vs, Real** fluxes, Real** qnx);

void cloud_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* cld_frac_l, Real* acn, Real* inv_dz,
  Real dt, Real inv_dt, bool do_predict_nc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* precip_liq_surf, Real* qc_tend, Real* nc_tend);

void ice_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* cld_frac_i, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qi, Real* qi_incld, Real* ni, Real* qm, Real* qm_incld, Real* bm, Real* bm_incld,
  Real* ni_incld, Real* precip_ice_surf, Real* qi_tend, Real* ni_tend);

void rain_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* cld_frac_r, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* precip_liq_surf, Real* precip_liq_flux, Real* qr_tend, Real* nr_tend);

void homogeneous_freezing_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* T_atm, Real* inv_exner, Real* latent_heat_fusion,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* th_atm);

void get_latent_heat_f(Int its, Int ite, Int kts, Int kte, Real* v, Real* s, Real* f);

void check_values_f(Real* Qv, Real* temp, Int kstart, Int kend,
                    Int timestepcount, bool force_abort, Int source_ind, Real* col_loc);

void p3_main_part1_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool do_predict_nc, bool do_prescribed_CCN,
  Real dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* nccn_prescribed, Real* inv_exner, Real* exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld,
  Real* qm_incld, Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld,
  bool* is_nucleat_possible, bool* is_hydromet_present);

void p3_main_part2_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool do_predict_nc, bool do_prescribed_CCN, Real dt, Real inv_dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* inv_exner, Real* exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i, Real* inv_cld_frac_r, Real* ni_activated, Real* inv_qc_relvar, Real* cld_frac_i, Real* cld_frac_l, Real* cld_frac_r, Real* qv_prev, Real* t_prev,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci, Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni,
  Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion, Real* qc_incld, Real* qr_incld, Real* qi_incld, Real* qm_incld, Real* nc_incld, Real* nr_incld,
  Real* ni_incld, Real* bm_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* qv2qi_depos_tend, Real* precip_total_tend,
  Real* nevapr, Real* qr_evap_tend, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* is_hydromet_present);

void p3_main_part3_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* inv_exner, Real* cld_frac_l, Real* cld_frac_r, Real* cld_frac_i,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vm_qi, Real* diag_eff_radius_qi, Real* diag_diam_qi, Real* rho_qi, Real* diag_equiv_reflectivity, Real* diag_eff_radius_qc);

Int p3_main_f(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th_atm, Real* qv, Real dt,
  Real* qi, Real* qm, Real* ni, Real* bm, Real* pres, Real* dz,
  Real* nc_nuceat_tend, Real* nccn_prescribed, Real* ni_activated, Real* inv_qc_relvar, Int it, Real* precip_liq_surf,
  Real* precip_ice_surf, Int its, Int ite, Int kts, Int kte, Real* diag_eff_radius_qc,
  Real* diag_eff_radius_qi, Real* rho_qi, bool do_predict_nc, bool do_prescribed_CCN, Real* dpres, Real* inv_exner,
  Real* qv2qi_depos_tend, Real* precip_liq_flux, Real* precip_ice_flux, Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i,
  Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* qv_prev, Real* t_prev);

} // end _f function decls

}  // namespace p3
}  // namespace scream

#endif
