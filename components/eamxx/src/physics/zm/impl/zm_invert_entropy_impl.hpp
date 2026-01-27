#ifndef ZM_INVERT_ENTROPY_IMPL_HPP
#define ZM_INVERT_ENTROPY_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm invert_entropy. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::invert_entropy(
  // Inputs
  const MemberType& team,
  const Real& s,
  const Real& p,
  const Real& qt,
  const Real& tfg,
  // Outputs
  Real& t,
  Real& qst)
{
  //----------------------------------------------------------------------------
  // Local variables
  bool converged;   // flag for convergence
  Real a, b, c, d, ebr, fa, fb, fc, pbr, qbr, rbr, sbr, xm;

  //----------------------------------------------------------------------------
  // initialize variables
  converged = false;
  t = tfg;            // first guess based on input temperature
  a = tfg-10;         // low bracket
  b = tfg+10;         // high bracket

  fa = entropy(a, p, qt) - s;
  fb = entropy(b, p, qt) - s;

  c = b;
  fc = fb;

  static constexpr Real half = 0.5;

  //----------------------------------------------------------------------------
  //
  for (int i = 0; i <= ZMC::LOOPMAX; ++i) {

    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c   = a;
      d   = b-a;
      fc  = fa;
      ebr = d;
    }

    if (std::abs(fc) < std::abs(fb)) {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    Real tolerance = 2*ZMC::tol_eps*std::abs(b) + half*ZMC::tol_coeff;
    xm = half*(c-b);

    converged = (std::abs(xm) <= tolerance || fb == 0);
    if (converged) break;

    if (std::abs(ebr) >= tolerance && std::abs(fa) > std::abs(fb)) {
      sbr=fb/fa;
      if (a == c) {
        pbr = 2*xm*sbr;
        qbr = 1 - sbr;
      }
      else {
        qbr = fa/fc;
        rbr = fb/fc;
        pbr = sbr*(2*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1));
        qbr = (qbr-1)*(rbr-1)*(sbr-1);
      }

      if (pbr > 0) qbr=-qbr;

      pbr=std::abs(pbr);
      if (2*pbr < ekat::impl::min(3*xm*qbr-std::abs(tolerance*qbr), std::abs(ebr*qbr))) {
        ebr = d;
        d = pbr/qbr;
      }
      else {
        d = xm;
        ebr = d;
      }
    }
    else {
      d = xm;
      ebr = d;
    }

    a = b;
    fa = fb;
    if (std::abs(d) > tolerance) {
      b += d;
    }
    else {
      b += Kokkos::copysign(tolerance, xm);
    }

    fb = entropy(b, p, qt) - s;
  }

  t = b;
  Real est; // saturation vapor pressure
  qsat_hPa(t, p, est, qst);

  EKAT_KERNEL_REQUIRE_MSG(converged, "ZM_CONV: INVERT_ENTROPY: Failed to converge");
}

} // namespace zm
} // namespace scream

#endif
