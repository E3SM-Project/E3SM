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
void advance_iop_nudging_c(Int plev, Real scm_dt, Real ps_in, Real* t_in, Real* q_in, Real* t_update, Real* q_update, Real* relaxt, Real* relaxq);
void advance_iop_subsidence_c(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* u_update, Real* v_update, Real* t_update, Real* q_update);
void iop_setinitial_c(Int nelemd, scream::dp::element_t* elem);
void iop_broadcast_c();
void apply_iop_forcing_c(Int nelemd, scream::dp::element_t* elem, scream::dp::hvcoord_t* hvcoord, scream::dp::hybrid_t* hybrid, scream::dp::timelevel_t* tl, Int n, bool t_before_advance, Int nets, Int nete);
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


void advance_iop_nudging(AdvanceIopNudgingData& d)
{
  dp_init(d.plev, true);
  advance_iop_nudging_c(d.plev, d.scm_dt, d.ps_in, d.t_in, d.q_in, d.t_update, d.q_update, d.relaxt, d.relaxq);
}

void advance_iop_subsidence(AdvanceIopSubsidenceData& d)
{
  dp_init(d.plev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  advance_iop_subsidence_c(d.plev, d.pcnst, d.scm_dt, d.ps_in, d.u_in, d.v_in, d.t_in, d.q_in, d.u_update, d.v_update, d.t_update, d.q_update);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void iop_setinitial(IopSetinitialData& d)
{
  dp_init(d.plev, true);
  iop_setinitial_c(d.nelemd, d.elem);
}

void iop_broadcast(IopBroadcastData& d)
{
  dp_init(d.plev, true);
  iop_broadcast_c();
}

void apply_iop_forcing(ApplyIopForcingData& d)
{
  dp_init(d.plev, true);
  apply_iop_forcing_c(d.nelemd, d.elem, &d.hvcoord, &d.hybrid, &d.tl, d.n, d.t_before_advance, d.nets, d.nete);
}

// end _c impls

//
// _f function definitions. These expect data in C layout
//

void advance_iop_forcing_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* t_phys_frc, Real* u_update, Real* v_update, Real* t_update, Real* q_update)
{
  // TODO
}
void advance_iop_nudging_f(Int plev, Real scm_dt, Real ps_in, Real* t_in, Real* q_in, Real* t_update, Real* q_update, Real* relaxt, Real* relaxq)
{
  // TODO
}
void advance_iop_subsidence_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* u_update, Real* v_update, Real* t_update, Real* q_update)
{
  // TODO
}
void iop_setinitial_f(Int nelemd, element_t* elem)
{
  // TODO
}
void iop_broadcast_f()
{
#if 0
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    PF::iop_broadcast();
  });
#endif

}
void apply_iop_forcing_f(Int nelemd, element_t* elem, hvcoord_t* hvcoord, hybrid_t hybrid, timelevel_t tl, Int n, bool t_before_advance, Int nets, Int nete)
{
  // TODO
}
} // namespace dp
} // namespace scream
