#ifndef ZM_ZM_CALC_FRACTIONAL_ENTRAINMENT_IMPL_HPP
#define ZM_ZM_CALC_FRACTIONAL_ENTRAINMENT_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_calc_fractional_entrainment. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_calc_fractional_entrainment(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& jb, // updraft base level
  const Int& jt, // updraft top level
  // Inputs/Outputs
  Int& j0, // level where updraft begins detraining
  // Inputs
  const uview_1d<const Real>& z_mid, // env altitude at mid-point
  const uview_1d<const Real>& z_int, // env altitude at interface
  const uview_1d<const Real>& dz, // layer thickness
  const uview_1d<const Real>& h_env, // env moist stat energy
  const uview_1d<const Real>& h_env_sat, // env saturated moist stat energy
  // Inputs/Outputs
  Real& h_env_min, // mid-tropospheric MSE minimum
  // Outputs
  const uview_1d<Real>& lambda, // fractional entrainment
  Real& lambda_max) // fractional entrainment maximum
{
  //----------------------------------------------------------------------------
  // Purpose: Determine properties of ZM updrafts and downdrafts
  //----------------------------------------------------------------------------
  // Local variables
  // variables used for Taylor series expansion when solving eq (4.78) for lamda_D
  // Note: k1, i2, i3, i4 are maintained as running accumulators (downward
  // k+1→k direction). lambda_tmp holds the computed lambda values before the
  // monotonicity pass and the final write into lambda.
  constexpr Real lambda_limit_min = 0.0;       // limiter
  constexpr Real lambda_limit_max = 0.0002;    // limiter
  constexpr Real lambda_threshold = 1.e-6;     // threshold for moving detrainment level

  uview_1d<Real> lambda_tmp, k1, i2, i3, i4;
  workspace.template take_many_contiguous_unsafe<5>(
    {"lambda_tmp", "k1", "i2", "i3", "i4"},
    {&lambda_tmp, &k1, &i2, &i3, &i4});

  //----------------------------------------------------------------------------
  // initialize variables
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&](const Int& k) {
    lambda(k)     = 0.0;
    lambda_tmp(k) = 0.0;
    k1(k) = 0.0;
    i2(k) = 0.0;
    i3(k) = 0.0;
    i4(k) = 0.0;
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // re-initialize minimum MSE for ensuing calculation
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    h_env_min = 1.e6;
  });
  team.team_barrier();

  Real h_env_min_new;
  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, msg, pver),
    [&](const Int& k, Real& lmin) {
      if (k >= j0 && k <= jb) {
        lmin = ekat::impl::min(lmin, h_env(k));
      }
    }, Kokkos::Min<Real>(h_env_min_new));
  team.team_barrier();

  Kokkos::single(Kokkos::PerTeam(team), [&] {
    h_env_min = ekat::impl::min(h_env_min, h_env_min_new);
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // compute taylor series for approximate lambda(z) below
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    for (Int k = pver-2; k >= msg; --k) {
      if (k < jb && k >= jt) {
        k1(k) = k1(k+1) + (h_env(jb)-h_env(k))*dz(k);
        const Real ihat = 0.5* (k1(k+1)+k1(k));
        i2(k) = i2(k+1) + ihat*dz(k);
        const Real idag = 0.5* (i2(k+1)+i2(k));
        i3(k) = i3(k+1) + idag*dz(k);
        const Real iprm = 0.5* (i3(k+1)+i3(k));
        i4(k) = i4(k+1) + iprm*dz(k);
      }
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // compute approximate lambda(z) using above taylor series - see eq (A6) in ZM95
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    for (Int k = msg+1; k < pver; ++k) {
      // --- lambda computation (Fortran: k from msg+2 to pver 1-based = k>=msg+1 0-based) ---
      Real expnum = 0.0;  // term for Taylor series
      if (k < jt || k >= jb) {
        expnum = 0.0;
      } else {
        expnum = h_env(jb) - (h_env_sat(k-1)*(z_int(k) - z_mid(k)) +
                              h_env_sat(k)*(z_mid(k-1) - z_int(k))) /
          (z_mid(k-1) - z_mid(k));
      }
      if ((h_env(jb) - h_env_min > 100.0 && expnum > 0.0) && k1(k) > expnum*dz(k)) {
        const Real tmp = expnum / k1(k);  // term for Taylor series
        lambda_tmp(k) = tmp +
          i2(k)/k1(k) * tmp*tmp +
          (2.0*i2(k)*i2(k) - k1(k)*i3(k))/(k1(k)*k1(k)) * tmp*tmp*tmp +
          (-5.0*k1(k)*i2(k)*i3(k) + 5.0*i2(k)*i2(k)*i2(k) + k1(k)*k1(k)*i4(k))/
          (k1(k)*k1(k)*k1(k)) * tmp*tmp*tmp*tmp;
        lambda_tmp(k) = ekat::impl::max(lambda_tmp(k), lambda_limit_min);
        lambda_tmp(k) = ekat::impl::min(lambda_tmp(k), lambda_limit_max);
      }
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // move detrainment level downward if fractional entrainment is too low
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    if (j0 < jb) {
      if (lambda_tmp(j0) < lambda_threshold && lambda_tmp(j0+1) > lambda_tmp(j0)) {
        j0 = j0 + 1;
      }
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // ensure that entrainment does not increase above the level that detrainment starts
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    for (Int k = msg+1; k < pver; ++k) {
      if (k >= jt && k <= j0) {
        lambda_tmp(k) = ekat::impl::max(lambda_tmp(k), lambda_tmp(k-1));
      }
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // specify maximum fractional entrainment
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    lambda_max = lambda_tmp(j0);
    lambda(jb) = lambda_max;
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // The modification below comes from:
  //   Rasch, P. J., J. E. Kristjánsson, A comparison of the CCM3 model climate
  //   using diagnosed and predicted condensate parameterizations, J. Clim., 1997.
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&](const Int& k) {
    if (k >= j0 && k <= jb) lambda(k) = lambda_tmp(j0);
    if (k < j0 && k >= jt)  lambda(k) = lambda_tmp(k);
  });
  team.team_barrier();

  workspace.template release_many_contiguous<5>({&lambda_tmp, &k1, &i2, &i3, &i4});
}

} // namespace zm
} // namespace scream

#endif
