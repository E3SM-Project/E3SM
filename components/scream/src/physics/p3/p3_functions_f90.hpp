#ifndef SCREAM_P3_FUNCTIONS_F90_HPP
#define SCREAM_P3_FUNCTIONS_F90_HPP

#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

#include "p3_functions.hpp"

#include <array>
#include <utility>

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
  Real qitot, nitot, qirim, rhop;

  // Outputs
  Int  dumi, dumjj, dumii, dumzz;
  Real dum1, dum4, dum5, dum6;
};
void find_lookuptable_indices_1a(LookupIceData& d);

extern "C" {

void find_lookuptable_indices_1a_f(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qitot, Real nitot, Real qirim, Real rhop);

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
  Real lcldm, rcldm, icldm;

  // in/out
  Real qcacc, qrevp, qcaut, ncacc, ncslf, ncautc, nrslf, nrevp, ncautr, qcnuc,
       ncnuc, qisub, nrshdr, qcheti, qrcol, qcshd, qimlt, qccol, qrheti, nimlt,
       nccol, ncshdc, ncheti, nrcol, nislf, qidep, nrheti, nisub, qinuc, ninuc,
       qiberg;

  // This populates all fields with test data within [0,1].
  void randomize();
};

void back_to_cell_average(BackToCellAverageData& d);

extern "C"{
  void back_to_cell_average_f(Real lcldm, Real rcldm, Real icldm,
                              Real* qcacc, Real* qrevp, Real* qcaut,
                              Real* ncacc, Real* ncslf, Real* ncautc,
                              Real* nrslf, Real* nrevp, Real* ncautr,
                              Real* qisub,
                              Real* nrshdr, Real* qcheti, Real* qrcol,
                              Real* qcshd, Real* qimlt, Real* qccol,
                              Real* qrheti, Real* nimlt, Real* nccol,
                              Real* ncshdc, Real* ncheti, Real* nrcol,
                              Real* nislf, Real* qidep, Real* nrheti,
                              Real* nisub, Real* qinuc, Real* ninuc,
                              Real* qiberg);
}

///////////////////////////////////////////////////////////////////////////////

struct PreventIceOverdepletionData
{
  // inputs
  Real pres, t, qv, xxls, odt;

  //output
  Real qidep, qisub;
};

void prevent_ice_overdepletion(PreventIceOverdepletionData& d);

extern "C"{
  void prevent_ice_overdepletion_f(Real pres, Real t, Real qv, Real xxls,
    Real odt, Real* qidep, Real* qisub);
}

///////////////////////////////////////////////////////////////////////////////

struct CloudWaterConservationData
{
  // inputs
  Real qc, dt;

  //output
  Real qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep;
};

void cloud_water_conservation(CloudWaterConservationData& d);

extern "C"{
  void cloud_water_conservation_f(Real qc, Real dt, Real* qcaut, Real* qcacc, Real* qccol,
    Real* qcheti, Real* qcshd, Real* qiberg, Real* qisub, Real* qidep);
}

struct RainWaterConservationData
{
  // inputs
  Real qr, qcaut, qcacc, qimlt, qcshd, dt;

  //output
  Real qrevp, qrcol, qrheti;
};

void rain_water_conservation(RainWaterConservationData& d);

extern "C"{
  void rain_water_conservation_f(Real qr, Real qcaut, Real qcacc, Real qimlt, Real qcshd,
  Real dt, Real* qrevp, Real* qrcol, Real* qrheti);
}

struct IceWaterConservationData
{
  //inputs
  Real qitot, qidep, qinuc, qiberg, qrcol, qccol, qrheti, qcheti, dt;

  //output
  Real qisub, qimlt;
};

void ice_water_conservation(IceWaterConservationData& d);

extern "C"{
  void ice_water_conservation_f(Real qitot, Real qidep, Real qinuc, Real qiberg, Real qrcol, Real qccol,
  Real qrheti, Real qcheti, Real dt, Real* qisub, Real* qimlt);
}
///////////////////////////////////////////////////////////////////////////////

struct CalcRimeDensityData
{
  // inputs
  Real t, rhofaci, f1pr02, acn, lamc, mu_c, qc_incld, qccol;

  // output
  Real vtrmi1, rhorime_c;
};

void calc_rime_density(CalcRimeDensityData& d);
extern "C"{
  void calc_rime_density_f(Real t, Real rhofaci, Real f1pr02, Real acn,
    Real lamc, Real mu_c, Real qc_incld, Real qccol, Real* vtrmi1,
    Real* rhorime_c);
}

///////////////////////////////////////////////////////////////////////////////

struct CldliqImmersionFreezingData
{
  // inputs
  Real t, lamc, mu_c, cdist1, qc_incld, qc_relvar;

  // output
  Real qcheti, ncheti;
};

void cldliq_immersion_freezing(CldliqImmersionFreezingData& d);
extern "C"{
  void cldliq_immersion_freezing_f(Real t, Real lamc, Real mu_c,
       Real cdist1, Real qc_incld, Real qc_relvar, Real* qcheti, Real* ncheti);
}

///////////////////////////////////////////////////////////////////////////////

struct RainImmersionFreezingData
{
  // inputs
  Real t, lamr, mu_r, cdistr, qr_incld;

  // output
  Real qrheti, nrheti;
};

void rain_immersion_freezing(RainImmersionFreezingData& d);
extern "C"{

  void rain_immersion_freezing_f(Real t, Real lamr, Real mu_r,
    Real cdistr, Real qr_incld, Real* qrheti, Real* nrheti);
}

///////////////////////////////////////////////////////////////////////////////

struct DropletSelfCollectionData
{
  // inputs
  Real rho, inv_rho, qc_incld, mu_c, nu, ncautc;

  // output
  Real ncslf;
};

void droplet_self_collection(DropletSelfCollectionData& d);
extern "C"{

  void droplet_self_collection_f(Real rho, Real inv_rho, Real qc_incld,
    Real mu_c, Real nu, Real ncautc, Real* ncslf);
}

///////////////////////////////////////////////////////////////////////////////

struct CloudRainAccretionData
{
  // inputs
  Real rho, inv_rho, qc_incld, nc_incld, qr_incld, qc_relvar;

  // output
  Real qcacc, ncacc;
};

void cloud_rain_accretion(CloudRainAccretionData& d);
extern "C"{

  void cloud_rain_accretion_f(Real rho, Real inv_rho, Real qc_incld,
       Real nc_incld, Real qr_incld, Real qc_relvar, Real* qcacc, Real* ncacc);
}

///////////////////////////////////////////////////////////////////////////////

struct CloudWaterAutoconversionData
{
  // inputs
  Real rho;
  Real qc_incld;
  Real nc_incld;
  Real qc_relvar;

  // output
  Real qcaut;
  Real ncautc;
  Real ncautr;
};

void cloud_water_autoconversion(CloudWaterAutoconversionData& d);
extern "C"{

  void cloud_water_autoconversion_f(Real rho, Real qc_incld, Real nc_incld, Real qc_relvar,
    Real* qcaut, Real* ncautc, Real* ncautr);
}

///////////////////////////////////////////////////////////////////////////////


struct RainSelfCollectionData
{
  //inputs
  Real rho, qr_incld, nr_incld;

  //output
  Real nrslf;
};

void rain_self_collection(RainSelfCollectionData& d);
extern "C"{

  void rain_self_collection_f(Real rho, Real qr_incld, Real nr_incld, Real* nrslf);
}

///////////////////////////////////////////////////////////////////////////////

struct ImposeMaxTotalNiData{
  // inout
  Real nitot_local;

  //input
  Real max_total_Ni, inv_rho_local;
};
void impose_max_total_Ni(ImposeMaxTotalNiData& d);
extern "C"{

  void impose_max_total_ni_f(Real* nitot_local, Real max_total_Ni, Real inv_rho_local);
}

///////////////////////////////////////////////////////////////////////////////

struct IceMeltingData
{
  // inputs
  Real rho,t,pres,rhofaci,f1pr05,f1pr14,xxlv,xlf,dv,sc,mu,kap,qv,qitot_incld,nitot_incld;

  // output
  Real qimlt,nimlt;
};

void ice_melting(IceMeltingData& d);

extern "C"{
void ice_melting_f(Real rho,Real t,Real pres,Real rhofaci,Real f1pr05,Real f1pr14,Real xxlv,Real xlf,Real dv,Real sc,Real mu,Real kap,Real qv,Real qitot_incld,Real nitot_incld,Real* qimlt,Real* nimlt);
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
  Real qc, rho, lcldm, nc_in;

  // Outputs
  Real nc_out, mu_c, nu, lamc, cdist, cdist1;
};
void get_cloud_dsd2(GetCloudDsd2Data& d);

extern "C" {

void get_cloud_dsd2_f(Real qc, Real* nc, Real* mu_c, Real rho, Real* nu, Real* lamc,
                      Real* cdist, Real* cdist1, Real lcldm);

}

//////////////////////////////////////////////////////////////////////////

struct GetRainDsd2Data
{
  // Inputs
  Real qr, rcldm, nr_in;

  // Outputs
  Real nr_out, lamr, mu_r, cdistr, logn0r;
};
void get_rain_dsd2(GetRainDsd2Data& d);

extern "C" {

void get_rain_dsd2_f(Real qr, Real* nr, Real* mu_r, Real* lamr, Real* cdistr, Real* logn0r, Real rcldm);

}

///////////////////////////////////////////////////////////////////////////////

struct CalcUpwindData
{
  // Inputs
  Int kts, kte, kdir, kbot, k_qxtop, num_arrays;
  Real dt_sub;
  Real* rho, *inv_rho, *inv_dzq;
  Real **vs;

  // In/out
  Real **qnx;

  // Outputs
  Real** fluxes;

  CalcUpwindData(Int kts_, Int kte_, Int kdir_, Int kbot_, Int k_qxtop_, Int num_arrays_, Real dt_sub_,
                 std::pair<Real, Real> rho_range, std::pair<Real, Real> inv_dzq_range,
                 std::pair<Real, Real> vs_range, std::pair<Real, Real> qnx_range);

  // deep copy
  CalcUpwindData(const CalcUpwindData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
  std::vector<Real*> m_ptr_data;
};
void calc_first_order_upwind_step(CalcUpwindData& d);

extern "C" {

void calc_first_order_upwind_step_f(Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub, Real* rho,
                                    Real* inv_rho, Real* inv_dzq, Int num_arrays, Real** fluxes, Real** vs, Real** qnx);

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
             Real prt_accum_, Int num_arrays_,
             std::pair<Real, Real> rho_range, std::pair<Real, Real> inv_dzq_range,
             std::pair<Real, Real> vs_range, std::pair<Real, Real> qnx_range);

};
void generalized_sedimentation(GenSedData& d);

extern "C" {

void generalized_sedimentation_f(Int kts, Int kte, Int kdir, Int k_qxtop, Int *k_qxbot, Int kbot, Real Co_max,
                                 Real* dt_left, Real* prt_accum, Real* inv_dzq, Real* inv_rho, Real* rho,
                                 Int num_arrays, Real** vs, Real** fluxes, Real** qnx);

}

///////////////////////////////////////////////////////////////////////////////

struct CloudSedData
{
  static constexpr size_t NUM_ARRAYS = 13;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real *qc_incld, *rho, *inv_rho, *lcldm, *acn, *inv_dzq;
  Real dt, odt;
  bool log_predictNc;

  // In/out
  Real *qc, *nc, *nc_incld, *mu_c, *lamc, *qc_tend, *nc_tend;
  Real prt_liq;

  CloudSedData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
               Real dt_, Real odt_, bool log_predictNc_, Real prt_liq_,
               const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

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
  Real* qc_incld, Real* rho, Real* inv_rho, Real* lcldm, Real* acn, Real* inv_dzq,
  Real dt, Real odt, bool log_predictNc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* prt_liq, Real* qc_tend, Real* nc_tend);

}

///////////////////////////////////////////////////////////////////////////////

struct IceSedData
{
  static constexpr size_t NUM_ARRAYS = 15;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real *rho, *inv_rho, *rhofaci, *icldm, *inv_dzq;
  Real dt, odt;

  // In/out
  Real *qitot, *qitot_incld, *nitot, *nitot_incld, *qirim, *qirim_incld, *birim, *birim_incld, *qi_tend, *ni_tend;
  Real prt_sol;

  IceSedData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
             Real dt_, Real odt_, Real prt_sol_,
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
  Real* rho, Real* inv_rho, Real* rhofaci, Real* icldm, Real* inv_dzq,
  Real dt, Real odt,
  Real* qitot, Real* qitot_incld, Real* nitot, Real* qirim, Real* qirim_incld, Real* birim, Real* birim_incld,
  Real* nitot_incld, Real* prt_sol, Real* qi_tend, Real* ni_tend);

}

///////////////////////////////////////////////////////////////////////////////

struct RainSedData
{
  static constexpr size_t NUM_ARRAYS = 14;

  // Inputs
  Int kts, kte, ktop, kbot, kdir;
  Real *rho, *inv_rho, *rhofacr, *rcldm, *inv_dzq, *qr_incld;
  Real dt, odt;

  // In/out
  Real *qr, *nr, *nr_incld, *mu_r, *lamr, *qr_tend, *nr_tend;
  Real *rflx; // has special size (nk+1)
  Real prt_liq;

  RainSedData(Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
              Real dt_, Real odt_, Real prt_liq_,
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
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* rcldm, Real* inv_dzq,
  Real dt, Real odt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* prt_liq, Real* rflx, Real* qr_tend, Real* nr_tend);

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
  Real* t, *exner, *xlf;

  // In/out
  Real* qc, *nc, *qr, *nr, *qitot, *nitot, *qirim, *birim, *th;

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
  Real* t, Real* exner, Real* xlf,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* th);

}

///////////////////////////////////////////////////////////////////////////////

struct ComputeRainFallVelocityData
{
  // Inputs
  Real qr_incld, rcldm, rhofacr;

  // In/out
  Real nr, nr_incld;

  // Outputs
  Real mu_r, lamr, V_qr, V_nr;
};
void compute_rain_fall_velocity(ComputeRainFallVelocityData& d);

extern "C" {

void compute_rain_fall_velocity_f(Real qr_incld, Real rcldm, Real rhofacr,
                                  Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* V_qr, Real* V_nr);

}
///////////////////////////////////////////////////////////////////////////////

struct GetTimeSpacePhysVarsData
{
  //Inputs
  Real t, pres, rho, xxlv, xxls, qvs, qvi;

  //Outs
  Real mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii;
};

void get_time_space_phys_variables(GetTimeSpacePhysVarsData& d);

extern "C"{

void get_time_space_phys_variables_f(Real t, Real pres, Real rho, Real xxlv, Real xxls, Real qvs, Real qvi,
				     Real* mu, Real* dv, Real* sc, Real* dqsdt, Real* dqsidt, Real* ab,
				     Real* abi, Real* kap, Real* eii);

}

///////////////////////////////////////////////////////////////////////////////

struct P3UpdatePrognosticIceData
{
  // Inputs
  Real qcheti, qccol, qcshd, nccol, ncheti, ncshdc, qrcol, nrcol, qrheti, nrheti, nrshdr, qimlt,
    nimlt, qisub, qidep, qinuc, ninuc, nislf, nisub, qiberg, exner, xxls, xlf;
  bool log_predictNc, log_wetgrowth;
  Real dt, nmltratio, rhorime_c;

  // In/outs
  Real th, qv, qitot, nitot, qirim, birim, qc, nc, qr, nr;
};

void update_prognostic_ice(P3UpdatePrognosticIceData& d);

extern "C"{

void update_prognostic_ice_f( Real qcheti,Real qccol, Real qcshd,  Real nccol,  Real ncheti, Real ncshdc,
Real qrcol,  Real nrcol, Real qrheti, Real nrheti, Real nrshdr, Real qimlt, Real nimlt, Real qisub,
Real qidep, Real qinuc, Real ninuc, Real nislf, Real nisub, Real qiberg, Real exner, Real xxls, Real xlf,
bool log_predictNc, bool log_wetgrowth, Real dt, Real nmltratio, Real rhorime_c, Real* th, Real* qv,
Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* qc, Real* nc, Real* qr, Real* nr);

}

///////////////////////////////////////////////////////////////////////////////

struct EvapSublimatePrecipData
{
  // Inputs
  Real qr_incld, qc_incld, nr_incld, qitot_incld, lcldm, rcldm, qvs, ab, epsr, qv;

  //Outs
  Real qrevp, nrevp;
};

void evaporate_sublimate_precip(EvapSublimatePrecipData& d);

extern "C"{

void evaporate_sublimate_precip_f( Real qr_incld, Real qc_incld, Real nr_incld, Real qitot_incld,
Real lcldm, Real rcldm, Real qvs, Real ab, Real epsr, Real qv, Real* qrevp, Real* nrevp);
}

///////////////////////////////////////////////////////////////////////////////

struct P3UpdatePrognosticLiqData
{
  // Inputs
  Real qcacc, ncacc, qcaut, ncautc, ncautr, ncslf, qrevp, nrevp, nrslf;

  bool log_predictNc;

  Real inv_rho, exner, xxlv, dt;

  // In/outs
  Real th, qv, qc, nc, qr, nr;
};

void update_prognostic_liquid(P3UpdatePrognosticLiqData& d);

extern "C"{

void update_prognostic_liquid_f( Real qcacc, Real ncacc, Real qcaut, Real ncautc, Real ncautr,
Real ncslf, Real  qrevp, Real nrevp, Real nrslf , bool log_predictNc,
Real inv_rho, Real exner, Real xxlv, Real dt, Real* th, Real* qv,
Real* qc, Real* nc, Real* qr, Real* nr);
}

  ///////////////////////////////////////////////////////////////////////////////

struct IceDepSublimationData
{
  //Inputs
  Real qitot_incld, nitot_incld, t, qvs, qvi, epsi, abi, qv;

  //Outs
  Real qidep, qisub, nisub, qiberg;
};

void ice_deposition_sublimation(IceDepSublimationData& d);

extern "C"{
void ice_deposition_sublimation_f( Real qitot_incld, Real nitot_incld, Real t, Real qvs, Real qvi,
Real epsi, Real abi, Real qv, Real* qidep, Real* qisub, Real* nisub, Real* qiberg);
}

struct IceCldliqCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, f1pr04, qitot_incld, qc_incld;
  Real nitot_incld, nc_incld;

  // Outputs
  Real qccol, nccol, qcshd, ncshdc;

};
void ice_cldliq_collection(IceCldliqCollectionData& d);

extern "C" {

void ice_cldliq_collection_f(Real rho, Real temp, Real rhofaci, Real f1pr04,
                             Real qitot_incld,Real qc_incld, Real nitot_incld, Real nc_incld,
                             Real* qccol, Real* nccol, Real* qcshd, Real* ncshdc);
}

struct IceRainCollectionData
{
  // Inputs
  Real rho, temp, rhofaci, logn0r, f1pr07, f1pr08, qitot_incld;
  Real nitot_incld, qr_incld;

  // Outputs
  Real qrcol, nrcol;

};
void ice_rain_collection(IceRainCollectionData& d);

extern "C" {

void ice_rain_collection_f(Real rho, Real temp, Real rhofaci, Real logn0r, Real f1pr07, Real f1pr08,
                         Real qitot_incld, Real nitot_incld, Real qr_incld, Real* qrcol, Real* nrcol);

}

struct IceSelfCollectionData
{
  // Inputs
  Real rho, rhofaci, f1pr03, eii, qirim_incld;
  Real qitot_incld, nitot_incld;

  // Outputs
  Real nislf;

};
void ice_self_collection(IceSelfCollectionData& d);

extern "C" {

void ice_self_collection_f(Real rho, Real rhofaci, Real f1pr03, Real eii,
                           Real qirim_incld, Real qitot_incld, Real nitot_incld, Real* nislf);

}

struct IceRelaxationData
{
  // Inputs
  Real rho, temp, rhofaci, f1pr05, f1pr14, dv, mu, sc, qitot_incld, nitot_incld;

  // Outputs
  Real epsi, epsi_tot;
};
void ice_relaxation_timescale(IceRelaxationData& d);

extern "C" {
 void ice_relaxation_timescale_f(Real rho, Real temp, Real rhofaci, Real f1pr05, Real f1pr14,
                                 Real dv, Real mu, Real sc, Real qitot_incld, Real nitot_incld,
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
  Real temp, inv_rho, nitot, naai, supi, odt;

  bool log_predictNc;

  // Outputs
  Real qinuc, ninuc;
};
void ice_nucleation(IceNucleationData& d);

extern "C" {
void ice_nucleation_f(Real temp, Real inv_rho, Real nitot, Real naai,
                      Real supi, Real odt, bool log_predictNc,
                      Real* qinuc, Real* ninuc);
}

struct IceWetGrowthData
{
  // Inputs
  Real rho, temp, pres, rhofaci, f1pr05, f1pr14, xxlv, xlf, dv, kap, mu, sc, qv, qc_incld;
  Real qitot_incld, nitot_incld, qr_incld;

  // In/Outs
  bool log_wetgrowth;

  Real qrcol, qccol, qwgrth, nrshdr, qcshd;
};
void ice_cldliq_wet_growth(IceWetGrowthData& d);

extern "C" {
void ice_cldliq_wet_growth_f(Real rho, Real temp, Real pres, Real rhofaci, Real f1pr05,
                             Real f1pr14, Real xxlv, Real xlf, Real dv,
                             Real kap, Real mu, Real sc, Real qv, Real qc_incld,
                             Real qitot_incld, Real nitot_incld, Real qr_incld, bool* log_wetgrowth,
                             Real* qrcol, Real* qccol, Real* qwgrth, Real* nrshdr, Real* qcshd);
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
  LatentHeatData& operator=(const LatentHeatData& rhs);

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
  Real qc, qr, qitot, qirim, nc, nr, nitot, birim, inv_lcldm, inv_icldm, inv_rcldm;

  // Outputs
  Real qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld;
};
void calculate_incloud_mixingratios(IncloudMixingData& d);

extern "C" {
void calculate_incloud_mixingratios_f(Real qc, Real qr, Real qitot, Real qirim, Real nc, Real nr, Real nitot, Real birim,
                                      Real inv_lcldm, Real inv_icldm, Real inv_rcldm,
                                      Real* qc_incld, Real* qr_incld, Real* qitot_incld, Real* qirim_incld,
                                      Real* nc_incld, Real* nr_incld, Real* nitot_incld, Real* birim_incld);
}

///////////////////////////////////////////////////////////////////////////////

struct P3MainPreLoopData
{
  static constexpr size_t NUM_ARRAYS = 39;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  bool log_predictNc;
  Real dt;
  Real* pres, *pdel, *dzq, *ncnuc, *exner, *inv_exner, *inv_lcldm, *inv_icldm, *inv_rcldm, *xxlv, *xxls, *xlf;

  // In/out
  Real* t, *rho, *inv_rho, *qvs, *qvi, *supi, *rhofacr, *rhofaci,
    *acn, *qv, *th, *qc, *nc, *qr, *nr, *qitot, *nitot, *qirim, *birim, *qc_incld, *qr_incld, *qitot_incld,
    *qirim_incld, *nc_incld, *nr_incld, *nitot_incld, *birim_incld;

  // Output
  bool log_nucleationPossible, log_hydrometeorsPresent;

  P3MainPreLoopData(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                    bool log_predictNc_, Real dt_,
                    const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  P3MainPreLoopData(const P3MainPreLoopData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void p3_main_pre_main_loop(P3MainPreLoopData& d);

extern "C" {

void p3_main_pre_main_loop_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool log_predictNc,
  Real dt,
  Real* pres, Real* pdel, Real* dzq, Real* ncnuc, Real* exner, Real* inv_exner, Real* inv_lcldm, Real* inv_icldm, Real* inv_rcldm, Real* xxlv, Real* xxls, Real* xlf,
  Real* t, Real* rho, Real* inv_rho, Real* qvs, Real* qvi, Real* supi, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* qc_incld, Real* qr_incld, Real* qitot_incld,
  Real* qirim_incld, Real* nc_incld, Real* nr_incld, Real* nitot_incld, Real* birim_incld,
  bool* log_nucleationPossible, bool* log_hydrometeorsPresent);

}

///////////////////////////////////////////////////////////////////////////////

struct P3MainLoopData
{
  static constexpr size_t NUM_ARRAYS = 62;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  bool log_predictNc;
  Real dt, odt;
  Real* pres, *pdel, *dzq, *ncnuc, *exner, *inv_exner, *inv_lcldm, *inv_icldm, *inv_rcldm, *naai, *qc_relvar, *icldm, *lcldm, *rcldm;

  // In/out
  Real* t, *rho, *inv_rho, *qvs, *qvi, *supi, *rhofacr, *rhofaci, *acn,
    *qv, *th, *qc, *nc, *qr, *nr, *qitot, *nitot, *qirim, *birim, *xxlv, *xxls, *xlf, *qc_incld, *qr_incld,
    *qitot_incld, *qirim_incld, *nc_incld, *nr_incld, *nitot_incld, *birim_incld, *mu_c, *nu, *lamc, *cdist, *cdist1,
    *cdistr, *mu_r, *lamr, *logn0r, *cmeiout, *prain, *nevapr, *prer_evap, *vap_liq_exchange,
    *vap_ice_exchange, *liq_ice_exchange, *pratot, *prctot;

  bool log_hydrometeorsPresent;

  P3MainLoopData(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                     bool log_predictNc_, Real dt_,
                     const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  P3MainLoopData(const P3MainLoopData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void p3_main_main_loop(P3MainLoopData& d);

extern "C" {

void p3_main_main_loop_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool log_predictNc, Real dt, Real odt,
  Real* pres, Real* pdel, Real* dzq, Real* ncnuc, Real* exner, Real* inv_exner, Real* inv_lcldm, Real* inv_icldm, Real* inv_rcldm, Real* naai, Real* qc_relvar, Real* icldm, Real* lcldm, Real* rcldm,
  Real* t, Real* rho, Real* inv_rho, Real* qvs, Real* qvi, Real* supi, Real* rhofacr, Real* rhofaci, Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot,
  Real* qirim, Real* birim, Real* xxlv, Real* xxls, Real* xlf, Real* qc_incld, Real* qr_incld, Real* qitot_incld, Real* qirim_incld, Real* nc_incld, Real* nr_incld,
  Real* nitot_incld, Real* birim_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* cmeiout, Real* prain,
  Real* nevapr, Real* prer_evap, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* log_hydrometeorsPresent);

}

///////////////////////////////////////////////////////////////////////////////

struct P3MainPostLoopData
{
  static constexpr size_t NUM_ARRAYS = 32;

  // Inputs
  Int kts, kte, kbot, ktop, kdir;
  Real* exner, *lcldm, *rcldm;

  // In/out
  Real* rho, *inv_rho, *rhofaci,
    *qv, *th, *qc, *nc, *qr, *nr, *qitot, *nitot, *qirim, *birim, *xxlv, *xxls,
    *mu_c, *nu, *lamc, *mu_r,
    *lamr, *vap_liq_exchange,
    *ze_rain, *ze_ice, *diag_vmi, *diag_effi, *diag_di, *diag_rhoi, *diag_ze, *diag_effc;

  P3MainPostLoopData(Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
                     const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges);

  // deep copy
  P3MainPostLoopData(const P3MainPostLoopData& rhs);

  Int nk() const { return m_nk; }

 private:
  // Internals
  Int m_nk;
  std::vector<Real> m_data;
};

void p3_main_post_main_loop(P3MainPostLoopData& d);

extern "C" {

void p3_main_post_main_loop_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* exner, Real* lcldm, Real* rcldm,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* xxlv, Real* xxls,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vmi, Real* diag_effi, Real* diag_di, Real* diag_rhoi, Real* diag_ze, Real* diag_effc);
}

///////////////////////////////////////////////////////////////////////////////
// BFB math stuff
///////////////////////////////////////////////////////////////////////////////

extern "C" {

Real cxx_pow(Real base, Real exp);
Real cxx_sqrt(Real base);
Real cxx_cbrt(Real base);
Real cxx_gamma(Real input);
Real cxx_log(Real input);
Real cxx_log10(Real input);
Real cxx_exp(Real input);

}

}  // namespace p3
}  // namespace scream

#endif
