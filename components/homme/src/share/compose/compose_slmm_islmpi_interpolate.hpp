#ifndef INCLUDE_COMPOSE_SLMM_ISLMPI_INTERPOLATE_HPP
#define INCLUDE_COMPOSE_SLMM_ISLMPI_INTERPOLATE_HPP

#include "compose.hpp"
#include "compose_slmm.hpp"
#include "compose_slmm_islmpi.hpp"

namespace slmm {

static constexpr Real sqrt5 = 2.23606797749978969641; // std::sqrt(5.0);
static constexpr Real oosqrt5 = 1.0 / sqrt5;

SLMM_KIF void gll_np4_eval (const Real x, Real y[4]) {
  static constexpr Real oo8 = 1.0/8.0;
  const Real x2 = x*x;
  y[0] = (1.0 - x)*(5.0*x2 - 1.0)*oo8;
  y[1] = -sqrt5*oo8*(sqrt5 - 5.0*x)*(x2 - 1.0);
  y[2] = -sqrt5*oo8*(sqrt5 + 5.0*x)*(x2 - 1.0);
  y[3] = (1.0 + x)*(5.0*x2 - 1.0)*oo8;
}

// Linear interp in each region.
SLMM_KIF void gll_np4_subgrid_eval_impl (const Real& x, Real y[4]) {
  if (x < -oosqrt5) {
    const Real alpha = (x + 1)/(1 - oosqrt5);
    y[0] = 1 - alpha;
    y[1] = alpha;
    y[2] = 0;
    y[3] = 0;
  } else {
    const Real alpha = (x + oosqrt5)/(2*oosqrt5);
    y[0] = 0;
    y[1] = 1 - alpha;
    y[2] = alpha;
    y[3] = 0;
  }
}

SLMM_KIF void gll_np4_subgrid_eval (const Real& x, Real y[4]) {
  if (x > 0) {
    gll_np4_subgrid_eval_impl(-x, y);
    ko::swap(y[0], y[3]);
    ko::swap(y[1], y[2]);    
    return;
  }
  gll_np4_subgrid_eval_impl(x, y);
}

// Quadratic interpolant across nodes 1,2,3 -- i.e., excluding node 0 -- of the
// np=4 reference element.
SLMM_KIF void outer_eval (const Real& x, Real v[4]) {
  static constexpr Real
    xbar = (2*oosqrt5) / (1 + oosqrt5),
    ooxbar = 1 / xbar,
    ybar = 1 / (xbar - 1);
  const Real xn = (x + oosqrt5) / (1 + oosqrt5);
  v[0] = 0;
  v[1] = 1 + ybar*xn*((1 - ooxbar)*xn + ooxbar - xbar);
  v[2] = ybar*ooxbar*xn*(xn - 1);
  v[3] = ybar*xn*(xbar - xn);
}

// In the middle region, use the standard GLL np=4 interpolant; in the two outer
// regions, use an order-reduced interpolant that stabilizes the method.
SLMM_KIF void gll_np4_subgrid_exp_eval (const Real& x, Real y[4]) {
  static constexpr Real
    alpha = 0.5527864045000416708,
    v = 0.427*(1 + alpha),
    x2 = 0.4472135954999579277,
    x3 = 1 - x2,
    det = x2*x3*(x2 - x3),
    y2 = alpha,
    y3 = v,
    c1 = (x3*y2 - x2*y3)/det,
    c2 = (-x3*x3*y2 + x2*x2*y3)/det;
  if (x < -oosqrt5 || x > oosqrt5) {
    if (x < -oosqrt5) {
      outer_eval(-x, y);
      ko::swap(y[0], y[3]);
      ko::swap(y[1], y[2]);
    } else
      outer_eval(x, y);
    Real y4[4];
    gll_np4_eval(x, y4);
    const Real x0 = 1 - std::abs(x);
    const Real a = (c1*x0 + c2)*x0;
    for (int i = 0; i < 4; ++i)
      y[i] = a*y[i] + (1 - a)*y4[i];
  } else
    gll_np4_eval(x, y);
}

} // namespace slmm

namespace homme {
namespace islmpi {

template <typename MT>
SLMM_KIF void interpolate (const typename IslMpi<MT>::Advecter::Alg::Enum& alg,
                           const Real ref_coord[2], Real rx[4], Real ry[4]) {
  typedef typename IslMpi<MT>::Advecter::Alg Alg;
  switch (alg) {
  case Alg::csl_gll:
    slmm::gll_np4_eval(ref_coord[0], rx);
    slmm::gll_np4_eval(ref_coord[1], ry);
    break;
  case Alg::csl_gll_subgrid:
    slmm::gll_np4_subgrid_eval(ref_coord[0], rx);
    slmm::gll_np4_subgrid_eval(ref_coord[1], ry);
    break;
  case Alg::csl_gll_exp:
    slmm::gll_np4_subgrid_exp_eval(ref_coord[0], rx);
    slmm::gll_np4_subgrid_exp_eval(ref_coord[1], ry);
    break;
  default:
    slmm_kernel_assert(0);
  }  
}

SLMM_KIF Real calc_q_tgt (const Real rx[4], const Real ry[4], const Real qs[16]) {
  return (ry[0]*(rx[0]*qs[ 0] + rx[1]*qs[ 1] + rx[2]*qs[ 2] + rx[3]*qs[ 3]) +
          ry[1]*(rx[0]*qs[ 4] + rx[1]*qs[ 5] + rx[2]*qs[ 6] + rx[3]*qs[ 7]) +
          ry[2]*(rx[0]*qs[ 8] + rx[1]*qs[ 9] + rx[2]*qs[10] + rx[3]*qs[11]) +
          ry[3]*(rx[0]*qs[12] + rx[1]*qs[13] + rx[2]*qs[14] + rx[3]*qs[15]));
}

SLMM_KIF Real calc_q_tgt (const Real rx[4], const Real ry[4], const Real qdp[16],
                          const Real dp[16]) {
  return (ry[0]*(rx[0]*(qdp[ 0]/dp[ 0]) + rx[1]*(qdp[ 1]/dp[ 1])  +
                 rx[2]*(qdp[ 2]/dp[ 2]) + rx[3]*(qdp[ 3]/dp[ 3])) +
          ry[1]*(rx[0]*(qdp[ 4]/dp[ 4]) + rx[1]*(qdp[ 5]/dp[ 5])  +
                 rx[2]*(qdp[ 6]/dp[ 6]) + rx[3]*(qdp[ 7]/dp[ 7])) +
          ry[2]*(rx[0]*(qdp[ 8]/dp[ 8]) + rx[1]*(qdp[ 9]/dp[ 9])  +
                 rx[2]*(qdp[10]/dp[10]) + rx[3]*(qdp[11]/dp[11])) +
          ry[3]*(rx[0]*(qdp[12]/dp[12]) + rx[1]*(qdp[13]/dp[13])  +
                 rx[2]*(qdp[14]/dp[14]) + rx[3]*(qdp[15]/dp[15])));
}

} // namespace islmpi
} // namespace homme

#endif
