#ifndef SCREAM_P3_FUNCTIONS_F90_HPP
#define SCREAM_P3_FUNCTIONS_F90_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

#include "p3_functions.hpp"

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
  using P3F = Functions<Real, HostDevice>;

  using view_itab_table = typename P3F::view_itab_table;
  using view_itabcol_table = typename P3F::view_itabcol_table;

  // All kokkos views must be destructed before Kokkos::finalize
  static void deinit();

  static const view_itab_table& itab()       { return get().m_itab; }
  static const view_itabcol_table& itabcol() { return get().m_itabcol; }

  P3GlobalForFortran() = delete;
  ~P3GlobalForFortran() = delete;
  P3GlobalForFortran(const P3GlobalForFortran&) = delete;
  P3GlobalForFortran& operator=(const P3GlobalForFortran&) = delete;

 private:
  struct Views {
    view_itab_table m_itab;
    view_itabcol_table m_itabcol;
  };

  static const Views& get();
  static std::shared_ptr<Views> s_views;
};

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

}  // namespace p3
}  // namespace scream

#endif
