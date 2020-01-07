#ifndef SCREAM_P3_FUNCTIONS_F90_HPP
#define SCREAM_P3_FUNCTIONS_F90_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

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
    view_2d_table m_vn_table, m_vm_table;
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

struct CloudWaterConservationData
{
  // inputs
  Real qc;
  Real qcnuc;
  Real dt;

  //output
  Real qcaut;
  Real qcacc;
  Real qccol;
  Real qcheti;
  Real qcshd;
  Real qiberg;
  Real qisub;
  Real qidep;
};

void cloud_water_conservation(CloudWaterConservationData& d);

extern "C"{
  void cloud_water_conservation_f(Real qc, Real qcnuc, Real dt, Real* qcaut, Real* qcacc, Real* qccol,
    Real* qcheti, Real* qcshd, Real* qiberg, Real* qisub, Real* qidep);
}

struct RainWaterConservationData
{
  // inputs
  Real qr;
  Real qcaut;
  Real qcacc;
  Real qimlt;
  Real qcshd;
  Real dt;

  //output
  Real qrevp;
  Real qrcol;
  Real qrheti;
};

void rain_water_conservation(RainWaterConservationData& d);

extern "C"{
  void rain_water_conservation_f(Real qr, Real qcaut, Real qcacc, Real qimlt, Real qcshd,
  Real dt, Real* qrevp, Real* qrcol, Real* qrheti);
}

struct IceWaterConservationData
{
  //inputs
  Real qitot;
  Real qidep;
  Real qinuc;
  Real qiberg;
  Real qrcol;
  Real qccol;
  Real qrheti;
  Real qcheti;
  Real dt;

  //output
  Real qisub;
  Real qimlt;

};

void ice_water_conservation(IceWaterConservationData& d);

extern "C"{
  void ice_water_conservation_f(Real qitot, Real qidep, Real qinuc, Real qiberg, Real qrcol, Real qccol,
  Real qrheti, Real qcheti, Real dt, Real* qisub, Real* qimlt);
}
///////////////////////////////////////////////////////////////////////////////

struct CloudWaterAutoconversionData
{
  // inputs
  Real rho;
  Real qc_incld;
  Real nc_incld;

  // output
  Real qcaut;
  Real ncautc;
  Real ncautr;
};

void cloud_water_autoconversion(CloudWaterAutoconversionData& d);
extern "C"{

  void cloud_water_autoconversion_f(Real rho, Real qc_incld, Real nc_incld, 
    Real* qcaut, Real* ncautc, Real* ncautr);
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

///////////////////////////////////////////////////////////////////////////////

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
// BFB math stuff
///////////////////////////////////////////////////////////////////////////////

extern "C" {

Real cxx_pow(Real base, Real exp);
Real cxx_cbrt(Real base);
Real cxx_gamma(Real input);
Real cxx_log(Real input);
Real cxx_log10(Real input);
Real cxx_exp(Real input);

}

}  // namespace p3
}  // namespace scream

#endif
