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

}  // namespace p3
}  // namespace scream

#endif
