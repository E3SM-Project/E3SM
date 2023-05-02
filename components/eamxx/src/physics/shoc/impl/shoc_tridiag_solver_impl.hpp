#ifndef SHOC_TRIDIAG_SOLVER_IMPL_HPP
#define SHOC_TRIDIAG_SOLVER_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "ekat/util/ekat_tridiag.hpp"

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_shoc_decomp(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& kv_term,
  const uview_1d<const Spack>& tmpi,
  const uview_1d<const Spack>& rdp_zt,
  const Scalar&                dtime,
  const Scalar&                flux,
  const uview_1d<Scalar>&       du,
  const uview_1d<Scalar>&       dl,
  const uview_1d<Scalar>&       d)
{
  const auto ggr = C::gravit;

  const auto skv_term = scalarize(kv_term);
  const auto stmpi = scalarize(tmpi);

  const Int nlev_pack = ekat::npack<Spack>(nlev);

  // Compute entries of the tridiagonal system
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

    // Compute local diagonals
    Spack du_k, dl_k, d_k;

    // Compute shift of kv_term and tmpi
    Spack kv_term_k, kv_term_kp1, tmpi_k, tmpi_kp1;
    const auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);

    // Original code was: auto shift_range = range_pack; but that caused mysterious test
    // failures on blake.
    IntSmallPack shift_range;
    vector_simd
    for (int s = 0; s < Spack::n; ++s) {
      shift_range[s] = range_pack[s];
    }

    shift_range.set(range_pack > nlev-1, 1); // don't calculate shift above nlev-1
    ekat::index_and_shift<1>(skv_term, shift_range, kv_term_k, kv_term_kp1);
    ekat::index_and_shift<1>(stmpi, shift_range, tmpi_k, tmpi_kp1);

    // Determine superdiagonal (du) and subdiagonal (dl) coeffs of the
    // tridiagonal diffusion matrix.
    du_k = -kv_term_kp1*tmpi_kp1*rdp_zt(k);
    dl_k = -kv_term(k)*tmpi(k)*rdp_zt(k);

    // The bottom element of the superdiagonal (du) and the top element of
    // the subdiagonal (dl) is set to zero (not included in linear system).
    du_k.set(range_pack == nlev-1, 0);
    dl_k.set(range_pack == 0, 0);

    // The diagonal elements are a combination of du and dl (d=1-du-dl). Surface
    // fluxes are applied explicitly in the diagonal at the top level.
    d_k = 1 - du_k - dl_k;
    d_k.set(range_pack == nlev-1, d_k + flux*dtime*ggr*rdp_zt(k));

    // Diagonals must be scalar in nlev.
    for (Int p=0; p<Spack::n && range_pack[p]<nlev; ++p) {
      du(range_pack[p]) = du_k[p];
      dl(range_pack[p]) = dl_k[p];
      d (range_pack[p]) = d_k [p];
    }
  });
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_shoc_solve(
  const MemberType&      team,
  const uview_1d<Scalar>& du,
  const uview_1d<Scalar>& dl,
  const uview_1d<Scalar>& d,
  const uview_2d<Spack>&  var)
{
#ifdef EKAT_DEFAULT_BFB
  ekat::tridiag::bfb(team, dl, d, du, var);
#else
#ifdef EAMXX_ENABLE_GPU
  ekat::tridiag::cr(team, dl, d, du, ekat::scalarize(var));
#else
  const auto f = [&] () { ekat::tridiag::thomas(dl, d, du, var); };
  Kokkos::single(Kokkos::PerTeam(team), f);
#endif
#endif
}

} // namespace shoc
} // namespace scream

#endif
