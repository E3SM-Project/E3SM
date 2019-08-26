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

struct P3InitAFortranData
{
  // Must use Host as device, f90 code might not be able to use Device memory
  using P3F = Functions<Real, HostDevice>;
  using view_itab_table = typename P3F::KT::template lview<Real[P3F::DENSIZE][P3F::RIMSIZE][P3F::ISIZE][P3F::TABSIZE]>;
  using view_itabcol_table = typename P3F::KT::template lview<Real[P3F::DENSIZE][P3F::RIMSIZE][P3F::ISIZE][P3F::RCOLLSIZE][P3F::COLTABSIZE]>;

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

  LookupIceData(Real qitot_, Real nitot_, Real qirim_, Real rhop_) :
    qitot(qitot_), nitot(nitot_), qirim(qirim_), rhop(rhop_) {}
};
void find_lookuptable_indices_1a(LookupIceData& d);

extern "C" {

void find_lookuptable_indices_1a_f_(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                    Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                    Real* qitot, Real* nitot, Real* qirim, Real* rhop);

}

struct LookupIceDataB
{
  // Inputs
  Real qr, nr;

  // Outputs
  Int dumj;
  Real dum3;

  LookupIceDataB(Real qr_, Real nr_) : qr(qr_), nr(nr_) {}
};
void find_lookuptable_indices_1b(LookupIceDataB& d);

struct AccessLookupTableData
{
  // Inputs
  LookupIceData& lid;
  Int index;

  // Outputs
  Real proc;

  AccessLookupTableData(LookupIceData& lid_, Int index_) : lid(lid_), index(index_) {}
};
void access_lookup_table(AccessLookupTableData& d);

struct AccessLookupTableCollData
{
  // Inputs
  LookupIceData& lid;
  LookupIceDataB& lidb;
  Int index;

  // Outputs
  Real proc;

  AccessLookupTableCollData(LookupIceData& lid_, LookupIceDataB& lidb_, Int index_) : lid(lid_), lidb(lidb_), index(index_) {}
};
void access_lookup_table_coll(AccessLookupTableCollData& d);

}  // namespace p3
}  // namespace scream

#endif
