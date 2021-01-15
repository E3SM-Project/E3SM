#ifndef INCLUDE_COMPOSE_CEDR_SL_HPP
#define INCLUDE_COMPOSE_CEDR_SL_HPP

#include "compose_homme.hpp"

namespace homme {
namespace sl {

using homme::Int;
using homme::Real;
using homme::FA2;
using homme::FA4;
using homme::FA5;

// Following are naming conventions in element_state and sl_advection:
//     elem(ie)%state%Q(:,:,k,q) is tracer mixing ratio.
//     elem(ie)%state%dp3d(:,:,k,tl%np1) is essentially total density.
//     elem(ie)%state%Qdp(:,:,k,q,n0_qdp) is Q*dp3d.
//     rho(:,:,k,ie) is spheremp*dp3d, essentially total mass at a GLL point.
//     Hence Q*rho = Q*spheremp*dp3d is tracer mass at a GLL point.
// We need to get pointers to some of these; elem can't be given the bind(C)
// attribute, so we can't take the elem array directly. We get these quantities
// at a mix of previous and current time steps.
//   In the code that follows, _p is previous and _c is current time step. Q is
// renamed to q, and Q is tracer mass at a GLL point.
struct Data {
  typedef std::shared_ptr<Data> Ptr;

  const TracerArrays<ko::MachineTraits>::Ptr ta;

  struct Check {
    Kokkos::View<Real**, Kokkos::Serial>
      mass_p, mass_c, mass_lo, mass_hi,
      q_lo, q_hi, q_min_l, q_max_l, qd_lo, qd_hi;
    Check (const Int nlev, const Int qsize)
      : mass_p("mass_p", nlev, qsize), mass_c("mass_c", nlev, qsize),
        mass_lo("mass_lo", nlev, qsize), mass_hi("mass_hi", nlev, qsize),
        q_lo("q_lo", nlev, qsize), q_hi("q_hi", nlev, qsize),
        q_min_l("q_min_l", nlev, qsize), q_max_l("q_max_l", nlev, qsize),
        qd_lo("qd_lo", nlev, qsize), qd_hi("qd_hi", nlev, qsize)
    {}
  };
  std::shared_ptr<Check> check;

  Data (const TracerArrays<ko::MachineTraits>::Ptr& tracer_arrays)
    : ta(tracer_arrays)
  {}
};

template <typename MT>
void run_global(CDR<MT>& cdr, const Data& d, Real* q_min_r, const Real* q_max_r,
                const Int nets, const Int nete);

template <typename MT>
void run_local(CDR<MT>& cdr, const Data& d, Real* q_min_r, const Real* q_max_r,
               const Int nets, const Int nete, const bool scalar_bounds,
               const Int limiter_option);

template <typename MT>
void check(CDR<MT>& cdr, Data& d, const Real* q_min_r, const Real* q_max_r,
           const Int nets, const Int nete);

} // namespace sl
} // namespace homme

#endif
