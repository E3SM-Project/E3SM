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
  using view_itab_table = typename P3F::view_itab_table;
  using view_itabcol_table = typename P3F::view_itabcol_table;
  using view_dnu_table = typename P3F::view_dnu_table;

  // All kokkos views must be destructed before Kokkos::finalize
  static void deinit();

  static const view_1d_table& mu_r_table()   { return get().m_mu_r_table; }
  static const view_2d_table& vn_table()     { return get().m_vn_table; }
  static const view_2d_table& vm_table()     { return get().m_vm_table; }
  static const view_2d_table& revap_table()  { return get().m_revap_table; }
  static const view_itab_table& itab()       { return get().m_itab; }
  static const view_itabcol_table& itabcol() { return get().m_itabcol; }
  static const view_dnu_table& dnu()         { return get().m_dnu; }

  P3GlobalForFortran() = delete;
  ~P3GlobalForFortran() = delete;
  P3GlobalForFortran(const P3GlobalForFortran&) = delete;
  P3GlobalForFortran& operator=(const P3GlobalForFortran&) = delete;

 private:
  struct Views {
    view_1d_table m_mu_r_table;
    view_2d_table m_vn_table, m_vm_table, m_revap_table;
    view_itab_table m_itab;
    view_itabcol_table m_itabcol;
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

  using view_itab_table = typename P3F::KT::template lview<Real[P3C::densize][P3C::rimsize][P3C::isize][P3C::tabsize]>;
  using view_itabcol_table = typename P3F::KT::template lview<Real[P3C::densize][P3C::rimsize][P3C::isize][P3C::rcollsize][P3C::coltabsize]>;

  // Need to be LayoutLeft to be fortran compatible
  view_itab_table itab;
  view_itabcol_table itabcol;

  P3InitAFortranData() :
    itab("P3InitAFortranData::itab"),
    itabcol("P3InitAFortranData::itabcol")
  {}
};
void p3_init_a(P3InitAFortranData& d);

///////////////////////////////////////////////////////////////////////////////

struct LookupIceData
{
  // Inputs
  Real qi, ni, qm, rhop;

  // Outputs
  Int  dumi, dumjj, dumii, dumzz;
  Real dum1, dum4, dum5, dum6;
};
void find_lookuptable_indices_1a(LookupIceData& d);

extern "C" {

void find_lookuptable_indices_1a_f(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qi, Real ni, Real qm, Real rhop);

}

///////////////////////////////////////////////////////////////////////////////

struct LookupIceDataB
{
  // Inputs
  Real qr, nr;

  // Outputs
  Int dumj;
  Real dum3;
};
void find_lookuptable_indices_1b(LookupIceDataB& d);

extern "C" {

void find_lookuptable_indices_1b_f(Int* dumj, Real* dum3, Real qr, Real nr);

}

///////////////////////////////////////////////////////////////////////////////

struct AccessLookupTableData
{
  // Inputs
  LookupIceData& lid;
  Int index;

  // Outputs
  Real proc;
};
void access_lookup_table(AccessLookupTableData& d);

extern "C" {

void access_lookup_table_f(Int dumjj, Int dumii, Int dumi, Int index,
                           Real dum1, Real dum4, Real dum5, Real* proc);

}

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
void access_lookup_table_coll(AccessLookupTableCollData& d);

extern "C" {

void access_lookup_table_coll_f(Int dumjj, Int dumii, Int dumj, Int dumi, Int index,
                                Real dum1, Real dum3, Real dum4, Real dum5, Real* proc);

}

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
  void randomize();
};

void back_to_cell_average(BackToCellAverageData& d);

extern "C"{
  void back_to_cell_average_f(Real cld_frac_l, Real cld_frac_r, Real cld_frac_i,
                              Real* qc2qr_accret_tend, Real* qr2qv_evap_tend, Real* qc2qr_autoconv_tend,
                              Real* nc_accret_tend, Real* nc_selfcollect_tend, Real* nc2nr_autoconv_tend,
                              Real* nr_selfcollect_tend, Real* nr_evap_tend, Real* ncautr,
                              Real* qi2qv_sublim_tend,
                              Real* nr_ice_shed_tend, Real* qc2qi_hetero_freeze_tend, Real* qr2qi_collect_tend,
                              Real* qc2qr_ice_shed_tend, Real* qi2qr_melt_tend, Real* qc2qi_collect_tend,
                              Real* qr2qi_immers_freeze_tend, Real* ni2nr_melt_tend, Real* nc_collect_tend,
                              Real* ncshdc, Real* nc2ni_immers_freeze_tend, Real* nr_collect_tend,
                              Real* ni_selfcollect_tend, Real* qv2qi_vapdep_tend, Real* nr2ni_immers_freeze_tend,
                              Real* ni_sublim_tend, Real* qv2qi_nucleat_tend, Real* ni_nucleat_tend,
                              Real* qc2qi_berg_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct PreventIceOverdepletionData
{
  // inputs
  Real pres, t, qv, latent_heat_sublim, inv_dt;

  //output
  Real qv2qi_vapdep_tend, qi2qv_sublim_tend;
};

void prevent_ice_overdepletion(PreventIceOverdepletionData& d);

extern "C"{
  void prevent_ice_overdepletion_f(Real pres, Real t, Real qv, Real latent_heat_sublim,
    Real inv_dt, Real* qv2qi_vapdep_tend, Real* qi2qv_sublim_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct CloudWaterConservationData
{
  // inputs
  Real qc, dt;

  //output
  Real qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend;
};

void cloud_water_conservation(CloudWaterConservationData& d);

extern "C"{
  void cloud_water_conservation_f(Real qc, Real dt, Real* qc2qr_autoconv_tend, Real* qc2qr_accret_tend, Real* qc2qi_collect_tend,
    Real* qc2qi_hetero_freeze_tend, Real* qc2qr_ice_shed_tend, Real* qc2qi_berg_tend, Real* qi2qv_sublim_tend, Real* qv2qi_vapdep_tend);
}

struct RainWaterConservationData
{
  // inputs
  Real qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt;

  //output
  Real qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend;
};

void rain_water_conservation(RainWaterConservationData& d);

extern "C"{
  void rain_water_conservation_f(Real qr, Real qc2qr_autoconv_tend, Real qc2qr_accret_tend, Real qi2qr_melt_tend, Real qc2qr_ice_shed_tend,
  Real dt, Real* qr2qv_evap_tend, Real* qr2qi_collect_tend, Real* qr2qi_immers_freeze_tend);
}

struct IceWaterConservationData
{
  //inputs
  Real qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, dt;

  //output
  Real qi2qv_sublim_tend, qi2qr_melt_tend;
};

void ice_water_conservation(IceWaterConservationData& d);

extern "C"{
  void ice_water_conservation_f(Real qi, Real qv2qi_vapdep_tend, Real qv2qi_nucleat_tend, Real qc2qi_berg_tend, Real qr2qi_collect_tend, Real qc2qi_collect_tend,
  Real qr2qi_immers_freeze_tend, Real qc2qi_hetero_freeze_tend, Real dt, Real* qi2qv_sublim_tend, Real* qi2qr_melt_tend);
}
///////////////////////////////////////////////////////////////////////////////

struct CalcRimeDensityData
{
  // inputs
  Real t, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend;

  // output
  Real vtrmi1, rho_qm_cloud;
};

void calc_rime_density(CalcRimeDensityData& d);
extern "C"{
  void calc_rime_density_f(Real t, Real rhofaci, Real table_val_qi_fallspd, Real acn,
    Real lamc, Real mu_c, Real qc_incld, Real qc2qi_collect_tend, Real* vtrmi1,
    Real* rho_qm_cloud);
}

///////////////////////////////////////////////////////////////////////////////

struct CldliqImmersionFreezingData
{
  // inputs
  Real t, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar;

  // output
  Real qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend;
};

void cldliq_immersion_freezing(CldliqImmersionFreezingData& d);
extern "C"{
  void cldliq_immersion_freezing_f(Real t, Real lamc, Real mu_c,
       Real cdist1, Real qc_incld, Real inv_qc_relvar, Real* qc2qi_hetero_freeze_tend, Real* nc2ni_immers_freeze_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct RainImmersionFreezingData
{
  // inputs
  Real t, lamr, mu_r, cdistr, qr_incld;

  // output
  Real qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend;
};

void rain_immersion_freezing(RainImmersionFreezingData& d);
extern "C"{

  void rain_immersion_freezing_f(Real t, Real lamr, Real mu_r,
    Real cdistr, Real qr_incld, Real* qr2qi_immers_freeze_tend, Real* nr2ni_immers_freeze_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct DropletSelfCollectionData
{
  // inputs
  Real rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend;

  // output
  Real nc_selfcollect_tend;
};

void droplet_self_collection(DropletSelfCollectionData& d);
extern "C"{

  void droplet_self_collection_f(Real rho, Real inv_rho, Real qc_incld,
    Real mu_c, Real nu, Real nc2nr_autoconv_tend, Real* nc_selfcollect_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct CloudRainAccretionData
{
  // inputs
  Real rho, inv_rho, qc_incld, nc_incld, qr_incld, inv_qc_relvar;

  // output
  Real qc2qr_accret_tend, nc_accret_tend;
};

void cloud_rain_accretion(CloudRainAccretionData& d);
extern "C"{

  void cloud_rain_accretion_f(Real rho, Real inv_rho, Real qc_incld,
       Real nc_incld, Real qr_incld, Real inv_qc_relvar, Real* qc2qr_accret_tend, Real* nc_accret_tend);
}

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

void cloud_water_autoconversion(CloudWaterAutoconversionData& d);
extern "C"{

  void cloud_water_autoconversion_f(Real rho, Real qc_incld, Real nc_incld, Real inv_qc_relvar,
    Real* qc2qr_autoconv_tend, Real* nc2nr_autoconv_tend, Real* ncautr);
}

///////////////////////////////////////////////////////////////////////////////


struct RainSelfCollectionData
{
  //inputs
  Real rho, qr_incld, nr_incld;

  //output
  Real nr_selfcollect_tend;
};

void rain_self_collection(RainSelfCollectionData& d);
extern "C"{

  void rain_self_collection_f(Real rho, Real qr_incld, Real nr_incld, Real* nr_selfcollect_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct ImposeMaxTotalNiData{
  // inout
  Real ni_local;

  //input
  Real max_total_Ni, inv_rho_local;
};
void impose_max_total_Ni(ImposeMaxTotalNiData& d);
extern "C"{

  void impose_max_total_ni_f(Real* ni_local, Real max_total_Ni, Real inv_rho_local);
}

///////////////////////////////////////////////////////////////////////////////

struct IceMeltingData
{
  // inputs
  Real rho,t,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld;

  // output
  Real qi2qr_melt_tend,ni2nr_melt_tend;
};

void ice_melting(IceMeltingData& d);

extern "C"{
void ice_melting_f(Real rho,Real t,Real pres,Real rhofaci,Real table_val_qi2qr_melting,Real table_val_qi2qr_vent_melt,Real latent_heat_vapor,Real latent_heat_fusion,Real dv,Real sc,Real mu,Real kap,Real qv,Real qi_incld,Real ni_incld,Real* qi2qr_melt_tend,Real* ni2nr_melt_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct SubgridVarianceScalingData
{
  // inputs
  Real relvar,expon;
  // no outputs - is a function
};

Real subgrid_variance_scaling(SubgridVarianceScalingData& d);

extern "C"{
  Real subgrid_variance_scaling_f(Real relvar,Real expon);
}

///////////////////////////////////////////////////////////////////////////////

struct GetCloudDsd2Data
{
  // Inputs
  Real qc, rho, cld_frac_l, nc_in;

  // Outputs
  Real nc_out, mu_c, nu, lamc, cdist, cdist1;
};
void get_cloud_dsd2(GetCloudDsd2Data& d);

extern "C" {

void get_cloud_dsd2_f(Real qc, Real* nc, Real* mu_c, Real rho, Real* nu, Real* lamc,
                      Real* cdist, Real* cdist1, Real cld_frac_l);

}

//////////////////////////////////////////////////////////////////////////

struct GetRainDsd2Data
{
  // Inputs
  Real qr, cld_frac_r, nr_in;

  // Outputs
  Real nr_out, lamr, mu_r, cdistr, logn0r;
};
void get_rain_dsd2(GetRainDsd2Data& d);

extern "C" {

void get_rain_dsd2_f(Real qr, Real* nr, Real* mu_r, Real* lamr, Real* cdistr, Real* logn0r, Real cld_frac_r);

}

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

  // deep copy
  CalcUpwindData(const CalcUpwindData& rhs);

  CalcUpwindData& operator=(const CalcUpwindData& rhs) = delete;

  Int nk() const { return m_nk; }

  void convert_to_ptr_arr(std::vector<Real*>& mem_space, Real**& fluxes_, Real**& vs_, Real**& qnx_);

 private:
  // Internals
  Int m_nk;
};
void calc_first_order_upwind_step(CalcUpwindData& d);

extern "C" {

void calc_first_order_upwind_step_f(Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub, Real* rho,
                                    Real* inv_rho, Real* inv_dz, Int num_arrays, Real** fluxes, Real** vs, Real** qnx);

}

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
};
void generalized_sedimentation(GenSedData& d);

extern "C" {

void generalized_sedimentation_f(Int kts, Int kte, Int kdir, Int k_qxtop, Int *k_qxbot, Int kbot, Real Co_max,
                                 Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
                                 Int num_arrays, Real** vs, Real** fluxes, Real** qnx);

}

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

  // deep copy
  CloudSedData(const CloudSedData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};
void cloud_sedimentation(CloudSedData& d);

extern "C" {

void cloud_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* cld_frac_l, Real* acn, Real* inv_dz,
  Real dt, Real inv_dt, bool do_predict_nc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* precip_liq_surf, Real* qc_tend, Real* nc_tend);

}

///////////////////////////////////////////////////////////////////////////////

struct IceSedData
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
             Real dt_, Real inv_dt_, Real precip_ice_surf_,
             const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  IceSedData(const IceSedData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void ice_sedimentation(IceSedData& d);

extern "C" {

void ice_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* cld_frac_i, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qi, Real* qi_incld, Real* ni, Real* qm, Real* qm_incld, Real* bm, Real* bm_incld,
  Real* ni_incld, Real* precip_ice_surf, Real* qi_tend, Real* ni_tend);

}

///////////////////////////////////////////////////////////////////////////////

struct RainSedData
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
              Real dt_, Real inv_dt_, Real precip_liq_surf_,
              const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  RainSedData(const RainSedData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void rain_sedimentation(RainSedData& d);

extern "C" {

void rain_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* cld_frac_r, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* precip_liq_surf, Real* precip_liq_flux, Real* qr_tend, Real* nr_tend);

}

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
void calc_bulk_rho_rime(CalcBulkRhoRimeData& d);

extern "C" {

void calc_bulk_rho_rime_f(Real qi_tot, Real* qi_rim, Real* bi_rim, Real* rho_rime);

}

///////////////////////////////////////////////////////////////////////////////

struct HomogeneousFreezingData
{
  static constexpr size_t NUM_ARRAYS = 12;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real* t, *exner, *latent_heat_fusion;

  // In/out
  Real* qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *th;

  HomogeneousFreezingData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
                          const std::array<std::pair<Real, Real>, NUM_ARRAYS>& ranges);

  // deep copy
  HomogeneousFreezingData(const HomogeneousFreezingData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;

};
void homogeneous_freezing(HomogeneousFreezingData& d);

extern "C" {

void homogeneous_freezing_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* t, Real* exner, Real* latent_heat_fusion,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* th);

}

///////////////////////////////////////////////////////////////////////////////

struct ComputeRainFallVelocityData
{
  // Inputs
  Real qr_incld, cld_frac_r, rhofacr;

  // In/out
  Real nr_incld;

  // Outputs
  Real mu_r, lamr, V_qr, V_nr;
};
void compute_rain_fall_velocity(ComputeRainFallVelocityData& d);

extern "C" {

void compute_rain_fall_velocity_f(Real qr_incld, Real cld_frac_r, Real rhofacr,
                                  Real* nr_incld, Real* mu_r, Real* lamr, Real* V_qr, Real* V_nr);

}
///////////////////////////////////////////////////////////////////////////////

struct GetTimeSpacePhysVarsData
{
  //Inputs
  Real t, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i;

  //Outs
  Real mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii;
};

void get_time_space_phys_variables(GetTimeSpacePhysVarsData& d);

extern "C"{

void get_time_space_phys_variables_f(Real t, Real pres, Real rho, Real latent_heat_vapor, Real latent_heat_sublim, Real qv_sat_l, Real qv_sat_i,
				     Real* mu, Real* dv, Real* sc, Real* dqsdt, Real* dqsidt, Real* ab,
				     Real* abi, Real* kap, Real* eii);

}

///////////////////////////////////////////////////////////////////////////////

struct P3UpdatePrognosticIceData
{
  // Inputs
  Real qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend, qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend,
    ni2nr_melt_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend, qc2qi_berg_tend, exner, latent_heat_sublim, latent_heat_fusion;
  bool do_predict_nc, log_wetgrowth;
  Real dt, nmltratio, rho_qm_cloud;

  // In/outs
  Real th, qv, qi, ni, qm, bm, qc, nc, qr, nr;
};

void update_prognostic_ice(P3UpdatePrognosticIceData& d);

extern "C"{

void update_prognostic_ice_f( Real qc2qi_hetero_freeze_tend,Real qc2qi_collect_tend, Real qc2qr_ice_shed_tend,  Real nc_collect_tend,  Real nc2ni_immers_freeze_tend, Real ncshdc,
Real qr2qi_collect_tend,  Real nr_collect_tend, Real qr2qi_immers_freeze_tend, Real nr2ni_immers_freeze_tend, Real nr_ice_shed_tend, Real qi2qr_melt_tend, Real ni2nr_melt_tend, Real qi2qv_sublim_tend,
Real qv2qi_vapdep_tend, Real qv2qi_nucleat_tend, Real ni_nucleat_tend, Real ni_selfcollect_tend, Real ni_sublim_tend, Real qc2qi_berg_tend, Real exner, Real latent_heat_sublim, Real latent_heat_fusion,
bool do_predict_nc, bool log_wetgrowth, Real dt, Real nmltratio, Real rho_qm_cloud, Real* th, Real* qv,
Real* qi, Real* ni, Real* qm, Real* bm, Real* qc, Real* nc, Real* qr, Real* nr);

}

///////////////////////////////////////////////////////////////////////////////

struct EvapSublimatePrecipData
{
  // Inputs
  Real qr_incld, qc_incld, nr_incld, qi_incld, cld_frac_l, cld_frac_r, qv_sat_l, ab, epsr, qv;

  //Outs
  Real qr2qv_evap_tend, nr_evap_tend;
};

void evaporate_sublimate_precip(EvapSublimatePrecipData& d);

extern "C"{

void evaporate_sublimate_precip_f( Real qr_incld, Real qc_incld, Real nr_incld, Real qi_incld,
Real cld_frac_l, Real cld_frac_r, Real qv_sat_l, Real ab, Real epsr, Real qv, Real* qr2qv_evap_tend, Real* nr_evap_tend);
}

///////////////////////////////////////////////////////////////////////////////

struct P3UpdatePrognosticLiqData
{
  // Inputs
  Real qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend;

  bool do_predict_nc;

  Real inv_rho, exner, latent_heat_vapor, dt;

  // In/outs
  Real th, qv, qc, nc, qr, nr;
};

void update_prognostic_liquid(P3UpdatePrognosticLiqData& d);

extern "C"{

void update_prognostic_liquid_f( Real qc2qr_accret_tend, Real nc_accret_tend, Real qc2qr_autoconv_tend, Real nc2nr_autoconv_tend, Real ncautr,
Real nc_selfcollect_tend, Real  qr2qv_evap_tend, Real nr_evap_tend, Real nr_selfcollect_tend , bool do_predict_nc,
Real inv_rho, Real exner, Real latent_heat_vapor, Real dt, Real* th, Real* qv,
Real* qc, Real* nc, Real* qr, Real* nr);
}

  ///////////////////////////////////////////////////////////////////////////////

struct IceDepSublimationData
{
  //Inputs
  Real qi_incld, ni_incld, t, qv_sat_l, qv_sat_i, epsi, abi, qv;

  //Outs
  Real qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend;
};

void ice_deposition_sublimation(IceDepSublimationData& d);

extern "C"{
void ice_deposition_sublimation_f( Real qi_incld, Real ni_incld, Real t, Real qv_sat_l, Real qv_sat_i,
Real epsi, Real abi, Real qv, Real* qv2qi_vapdep_tend, Real* qi2qv_sublim_tend, Real* ni_sublim_tend, Real* qc2qi_berg_tend);
}

struct IceCldliqCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, table_val_qc2qi_collect, qi_incld, qc_incld;
  Real ni_incld, nc_incld;

  // Outputs
  Real qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc;

};
void ice_cldliq_collection(IceCldliqCollectionData& d);

extern "C" {

void ice_cldliq_collection_f(Real rho, Real temp, Real rhofaci, Real table_val_qc2qi_collect,
                             Real qi_incld,Real qc_incld, Real ni_incld, Real nc_incld,
                             Real* qc2qi_collect_tend, Real* nc_collect_tend, Real* qc2qr_ice_shed_tend, Real* ncshdc);
}

struct IceRainCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect, qi_incld;
  Real ni_incld, qr_incld;

  // Outputs
  Real qr2qi_collect_tend, nr_collect_tend;

};
void ice_rain_collection(IceRainCollectionData& d);

extern "C" {

void ice_rain_collection_f(Real rho, Real temp, Real rhofaci, Real logn0r, Real table_val_nr_collect, Real table_val_qr2qi_collect,
                         Real qi_incld, Real ni_incld, Real qr_incld, Real* qr2qi_collect_tend, Real* nr_collect_tend);

}

struct IceSelfCollectionData
{
  // Inputs
  Real rho, rhofaci, table_val_ni_self_collect, eii, qm_incld;
  Real qi_incld, ni_incld;

  // Outputs
  Real ni_selfcollect_tend;

};
void ice_self_collection(IceSelfCollectionData& d);

extern "C" {

void ice_self_collection_f(Real rho, Real rhofaci, Real table_val_ni_self_collect, Real eii,
                           Real qm_incld, Real qi_incld, Real ni_incld, Real* ni_selfcollect_tend);

}

struct IceRelaxationData
{
  // Inputs
  Real rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, mu, sc, qi_incld, ni_incld;

  // Outputs
  Real epsi, epsi_tot;
};
void ice_relaxation_timescale(IceRelaxationData& d);

extern "C" {
 void ice_relaxation_timescale_f(Real rho, Real temp, Real rhofaci, Real table_val_qi2qr_melting, Real table_val_qi2qr_vent_melt,
                                 Real dv, Real mu, Real sc, Real qi_incld, Real ni_incld,
                                 Real* epsi, Real* epsi_tot);
}

struct CalcLiqRelaxationData
{
  // Inputs
  Real rho, f1r, f2r, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld;

  // Outputs
  Real epsr, epsc;

  // This populates all input fields with test data within [0,1].
  void randomize();
};
void calc_liq_relaxation_timescale(CalcLiqRelaxationData& d);

extern "C" {
void calc_liq_relaxation_timescale_f(Real rho, Real f1r, Real f2r, Real dv,
                                     Real mu, Real sc, Real mu_r, Real lamr,
                                     Real cdistr, Real cdist, Real qr_incld,
                                     Real qc_incld, Real* epsr, Real* epsc);
}

struct IceNucleationData
{
  // Inputs
  Real temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt;

  bool do_predict_nc;

  // Outputs
  Real qv2qi_nucleat_tend, ni_nucleat_tend;
};
void ice_nucleation(IceNucleationData& d);

extern "C" {
void ice_nucleation_f(Real temp, Real inv_rho, Real ni, Real ni_activated,
                      Real qv_supersat_i, Real inv_dt, bool do_predict_nc,
                      Real* qv2qi_nucleat_tend, Real* ni_nucleat_tend);
}

struct IceWetGrowthData
{
  // Inputs
  Real rho, temp, pres, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, kap, mu, sc, qv, qc_incld;
  Real qi_incld, ni_incld, qr_incld;

  // In/Outs
  bool log_wetgrowth;

  Real qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend;
};
void ice_cldliq_wet_growth(IceWetGrowthData& d);

extern "C" {
void ice_cldliq_wet_growth_f(Real rho, Real temp, Real pres, Real rhofaci, Real table_val_qi2qr_melting,
                             Real table_val_qi2qr_vent_melt, Real latent_heat_vapor, Real latent_heat_fusion, Real dv,
                             Real kap, Real mu, Real sc, Real qv, Real qc_incld,
                             Real qi_incld, Real ni_incld, Real qr_incld, bool* log_wetgrowth,
                             Real* qr2qi_collect_tend, Real* qc2qi_collect_tend, Real* qc_growth_rate, Real* nr_ice_shed_tend, Real* qc2qr_ice_shed_tend);
}

struct LatentHeatData
{
  static constexpr size_t NUM_ARRAYS = 3;

  // Inputs
  Int its, ite, kts, kte;

  // Outputs
  Real* v, *s, *f;

  LatentHeatData(Int its_, Int ite_, Int kts_, Int kte_);
  LatentHeatData(const LatentHeatData& rhs);
  LatentHeatData& operator=(const LatentHeatData& rhs) = delete;

  void transpose();
  void init_ptrs();

  // Internals
  Int m_ni, m_nk, m_total;
  std::vector<Real> m_data;
};
void get_latent_heat(LatentHeatData& d);

extern "C" {
void get_latent_heat_f(Int its, Int ite, Int kts, Int kte, Real* v, Real* s, Real* f);
}

struct CheckValuesData
{
  static constexpr size_t NUM_ARRAYS = 2;

  // Inputs
   Int kts, kte;
   Int timestepcount, source_ind;

   bool force_abort;

   Real *qv, *temp, *col_loc;

   CheckValuesData(Int kts_, Int kte_, Int timestepcount_, Int source_ind_, bool force_abort_,
                   const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy

  CheckValuesData(const CheckValuesData& rhs);

  Int nk() const { return m_nk; }

  private:
  // Internals

  Int m_nk;
  std::vector<Real> m_data;
};
void check_values(CheckValuesData& d);

extern "C" {
void check_values_f(Real* Qv, Real* temp, Int kstart, Int kend,
                    Int timestepcount, bool force_abort, Int source_ind, Real* col_loc);
}

struct IncloudMixingData
{
  // Inputs
  Real qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r;

  // Outputs
  Real qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld;
};
void calculate_incloud_mixingratios(IncloudMixingData& d);

extern "C" {
void calculate_incloud_mixingratios_f(Real qc, Real qr, Real qi, Real qm, Real nc, Real nr, Real ni, Real bm,
                                      Real inv_cld_frac_l, Real inv_cld_frac_i, Real inv_cld_frac_r,
                                      Real* qc_incld, Real* qr_incld, Real* qi_incld, Real* qm_incld,
                                      Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld);
}

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart1Data
{
  static constexpr size_t NUM_ARRAYS = 39;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  bool do_predict_nc;
  Real dt;
  Real* pres, *dpres, *dz, *nc_nuceat_tend, *exner, *inv_exner, *inv_cld_frac_l, *inv_cld_frac_i, *inv_cld_frac_r, *latent_heat_vapor, *latent_heat_sublim, *latent_heat_fusion;

  // In/out
  Real* t, *rho, *inv_rho, *qv_sat_l, *qv_sat_i, *qv_supersat_i, *rhofacr, *rhofaci,
    *acn, *qv, *th, *qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *qc_incld, *qr_incld, *qi_incld,
    *qm_incld, *nc_incld, *nr_incld, *ni_incld, *bm_incld;

  // Output
  bool is_nucleat_possible, is_hydromet_present;

  P3MainPart1Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                    bool do_predict_nc_, Real dt_,
                    const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  P3MainPart1Data(const P3MainPart1Data& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void p3_main_part1(P3MainPart1Data& d);

extern "C" {

void p3_main_part1_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool do_predict_nc,
  Real dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* exner, Real* inv_exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i, Real* inv_cld_frac_r, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion,
  Real* t, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld,
  Real* qm_incld, Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld,
  bool* is_nucleat_possible, bool* is_hydromet_present);

}

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart2Data
{
  static constexpr size_t NUM_ARRAYS = 62;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  bool do_predict_nc;
  Real dt, inv_dt;
  Real* pres, *dpres, *dz, *nc_nuceat_tend, *exner, *inv_exner, *inv_cld_frac_l, *inv_cld_frac_i, *inv_cld_frac_r, *ni_activated, *inv_qc_relvar, *cld_frac_i, *cld_frac_l, *cld_frac_r;

  // In/out
  Real* t, *rho, *inv_rho, *qv_sat_l, *qv_sat_i, *qv_supersat_i, *rhofacr, *rhofaci, *acn,
    *qv, *th, *qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *latent_heat_vapor, *latent_heat_sublim, *latent_heat_fusion, *qc_incld, *qr_incld,
    *qi_incld, *qm_incld, *nc_incld, *nr_incld, *ni_incld, *bm_incld, *mu_c, *nu, *lamc, *cdist, *cdist1,
    *cdistr, *mu_r, *lamr, *logn0r, *cmeiout, *precip_total_tend, *nevapr, *qr_evap_tend, *vap_liq_exchange,
    *vap_ice_exchange, *liq_ice_exchange, *pratot, *prctot;

  bool is_hydromet_present;

  P3MainPart2Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                     bool do_predict_nc_, Real dt_,
                     const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  P3MainPart2Data(const P3MainPart2Data& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void p3_main_part2(P3MainPart2Data& d);

extern "C" {

void p3_main_part2_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool do_predict_nc, Real dt, Real inv_dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* exner, Real* inv_exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i, Real* inv_cld_frac_r, Real* ni_activated, Real* inv_qc_relvar, Real* cld_frac_i, Real* cld_frac_l, Real* cld_frac_r,
  Real* t, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci, Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni,
  Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion, Real* qc_incld, Real* qr_incld, Real* qi_incld, Real* qm_incld, Real* nc_incld, Real* nr_incld,
  Real* ni_incld, Real* bm_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* cmeiout, Real* precip_total_tend,
  Real* nevapr, Real* qr_evap_tend, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* is_hydromet_present);

}

///////////////////////////////////////////////////////////////////////////////

struct P3MainPart3Data
{
  static constexpr size_t NUM_ARRAYS = 32;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  Real* exner, *cld_frac_l, *cld_frac_r;

  // In/out
  Real* rho, *inv_rho, *rhofaci,
    *qv, *th, *qc, *nc, *qr, *nr, *qi, *ni, *qm, *bm, *latent_heat_vapor, *latent_heat_sublim,
    *mu_c, *nu, *lamc, *mu_r,
    *lamr, *vap_liq_exchange,
    *ze_rain, *ze_ice, *diag_vmi, *diag_effi, *diag_di, *rho_qi, *diag_ze, *diag_effc;

  P3MainPart3Data(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                     const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  P3MainPart3Data(const P3MainPart3Data& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void p3_main_part3(P3MainPart3Data& d);

extern "C" {

void p3_main_part3_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* exner, Real* cld_frac_l, Real* cld_frac_r,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vmi, Real* diag_effi, Real* diag_di, Real* rho_qi, Real* diag_ze, Real* diag_effc);
}

///////////////////////////////////////////////////////////////////////////////

struct P3MainData
{
  static constexpr size_t NUM_ARRAYS = 36;
  static constexpr size_t NUM_INPUT_ARRAYS = 20;

  // Inputs
  Int its, ite, kts, kte, it;
  Real* pres, *dz, *nc_nuceat_tend, *ni_activated, *dpres, *exner, *cld_frac_i, *cld_frac_l, *cld_frac_r, *inv_qc_relvar;
  Real dt;
  bool do_predict_nc;

  // In/out
  Real* qc, *nc, *qr, *nr, *qi, *qm, *ni, *bm, *qv, *th;

  // Out
  Real *diag_effc, *diag_effi, *rho_qi, *mu_c, *lamc, *cmeiout, *precip_total_tend, *nevapr,
       *qr_evap_tend, *liq_ice_exchange, *vap_liq_exchange, *vap_ice_exchange,
       *precip_liq_flux, *precip_ice_flux, *precip_liq_surf, *precip_ice_surf;

  P3MainData(Int its_, Int ite_, Int kts_, Int kte_, Int it_, Real dt_, bool do_predict_nc_,
             const std::array< std::pair<Real, Real>, NUM_INPUT_ARRAYS >& ranges);

  template <ekat::util::TransposeDirection::Enum D>
  void transpose()
  {
    using ekat::util::transpose;

    P3MainData d_trans(*this);

    transpose<D>(pres, d_trans.pres, m_ni, m_nk);
    transpose<D>(dz, d_trans.dz, m_ni, m_nk);
    transpose<D>(nc_nuceat_tend, d_trans.nc_nuceat_tend, m_ni, m_nk);
    transpose<D>(ni_activated, d_trans.ni_activated, m_ni, m_nk);
    transpose<D>(dpres, d_trans.dpres, m_ni, m_nk);
    transpose<D>(exner, d_trans.exner, m_ni, m_nk);
    transpose<D>(cld_frac_i, d_trans.cld_frac_i, m_ni, m_nk);
    transpose<D>(cld_frac_l, d_trans.cld_frac_l, m_ni, m_nk);
    transpose<D>(cld_frac_r, d_trans.cld_frac_r, m_ni, m_nk);
    transpose<D>(inv_qc_relvar, d_trans.inv_qc_relvar, m_ni, m_nk);
    transpose<D>(qc, d_trans.qc, m_ni, m_nk);
    transpose<D>(nc, d_trans.nc, m_ni, m_nk);
    transpose<D>(qr, d_trans.qr, m_ni, m_nk);
    transpose<D>(nr, d_trans.nr, m_ni, m_nk);
    transpose<D>(qi, d_trans.qi, m_ni, m_nk);
    transpose<D>(qm, d_trans.qm, m_ni, m_nk);
    transpose<D>(ni, d_trans.ni, m_ni, m_nk);
    transpose<D>(bm, d_trans.bm, m_ni, m_nk);
    transpose<D>(qv, d_trans.qv, m_ni, m_nk);
    transpose<D>(th, d_trans.th, m_ni, m_nk);
    transpose<D>(diag_effc, d_trans.diag_effc, m_ni, m_nk);
    transpose<D>(diag_effi, d_trans.diag_effi, m_ni, m_nk);
    transpose<D>(rho_qi, d_trans.rho_qi, m_ni, m_nk);
    transpose<D>(mu_c, d_trans.mu_c, m_ni, m_nk);
    transpose<D>(lamc, d_trans.lamc, m_ni, m_nk);
    transpose<D>(cmeiout, d_trans.cmeiout, m_ni, m_nk);
    transpose<D>(precip_total_tend, d_trans.precip_total_tend, m_ni, m_nk);
    transpose<D>(nevapr, d_trans.nevapr, m_ni, m_nk);
    transpose<D>(qr_evap_tend, d_trans.qr_evap_tend, m_ni, m_nk);
    transpose<D>(liq_ice_exchange, d_trans.liq_ice_exchange, m_ni, m_nk);
    transpose<D>(vap_liq_exchange, d_trans.vap_liq_exchange, m_ni, m_nk);
    transpose<D>(vap_ice_exchange, d_trans.vap_ice_exchange, m_ni, m_nk);
    transpose<D>(precip_liq_flux, d_trans.precip_liq_flux, m_ni, m_nk+1);
    transpose<D>(precip_ice_flux, d_trans.precip_ice_flux, m_ni, m_nk+1);

    *this = std::move(d_trans);
  }

  // deep copy
  P3MainData(const P3MainData& rhs);

  P3MainData& operator=(P3MainData&&) = default;

  Int nt() const { return m_nt; }

 private:
  // Internals
  Int m_ni, m_nk, m_nt;
  std::vector<Real> m_data;
};

void p3_main(P3MainData& d);

extern "C" {

void p3_main_f(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th, Real* qv, Real dt,
  Real* qi, Real* qm, Real* ni, Real* bm, Real* pres, Real* dz,
  Real* nc_nuceat_tend, Real* ni_activated, Real* inv_qc_relvar, Int it, Real* precip_liq_surf,
  Real* precip_ice_surf, Int its, Int ite, Int kts, Int kte, Real* diag_effc,
  Real* diag_effi, Real* rho_qi, bool do_predict_nc, Real* dpres, Real* exner,
  Real* cmeiout, Real* precip_total_tend, Real* nevapr, Real* qr_evap_tend, Real* precip_liq_flux,
  Real* precip_ice_flux, Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i, Real* mu_c, Real* lamc,
  Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange);
}

}  // namespace p3
}  // namespace scream

#endif
