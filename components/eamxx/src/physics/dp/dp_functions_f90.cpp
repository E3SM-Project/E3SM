#include "dp_functions_f90.hpp"

#include "dp_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include "share/util/scream_deep_copy.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C interface to DP fortran calls. The stubs below will link to fortran definitions in dp_iso_c.f90
//

extern "C" {

void advance_iop_forcing_c(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* t_phys_frc, Real* u_update, Real* v_update, Real* t_update, Real* q_update);
} // extern "C" : end _c decls

namespace scream {
namespace dp {

//
// Glue functions to call fortran from from C++ with the Data struct
//

void advance_iop_forcing(AdvanceIopForcingData& d)
{
  dp_init(d.plev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  advance_iop_forcing_c(d.plev, d.pcnst, d.scm_dt, d.ps_in, d.u_in, d.v_in, d.t_in, d.q_in, d.t_phys_frc, d.u_update, d.v_update, d.t_update, d.q_update);
  d.transpose<ekat::TransposeDirection::f2c>();
}


// end _c impls

//
// _f function definitions. These expect data in C layout
//

void advance_iop_forcing_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* t_phys_frc, Real* u_update, Real* v_update, Real* t_update, Real* q_update)
{
  // TODO
}
} // namespace dp
} // namespace scream
