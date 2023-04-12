#ifndef SHOC_PBLINTD_HEIGHT_IMPL_HPP
#define SHOC_PBLINTD_HEIGHT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd_height. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

// PBL height calculation: Scan upward until the Richardson number between
// the first level and the current level exceeds the "critical" value.

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd_height(
  const MemberType& team,
  const Int& nlev,
  const Int& npbl,
  const uview_1d<const Spack>& z,
  const uview_1d<const Spack>& u,
  const uview_1d<const Spack>& v,
  const Scalar& ustar,
  const uview_1d<const Spack>& thv,
  const Scalar& thv_ref,
  Scalar& pblh,
  const uview_1d<Spack>& rino,
  bool& check)
{
  if (!check)
    return;

  const auto ggr = C::gravit;
  const Scalar tiny = 1e-36;
  const Scalar fac = 100;
  const Scalar ricr = 0.3;

  // Scalarize views for single entry access
  const auto s_u    = ekat::scalarize(u);
  const auto s_v    = ekat::scalarize(v);
  const auto s_z    = ekat::scalarize(z);
  const auto s_thv  = ekat::scalarize(thv);
  const auto s_rino = ekat::scalarize(rino);

  // Compute range for reduction. Only run the
  // reduction if there is a valid range
  const Int lower_indx = nlev-npbl;
  const Int upper_indx = nlev-1;
  if (lower_indx >= upper_indx) {
    return;
  }
  const Int lower_pack_indx = lower_indx/Spack::n;
  const Int upper_pack_indx = upper_indx/Spack::n;

  // Compute rino values and find max index k s.t. rino(k) >= ricr
  Int max_indx = Kokkos::reduction_identity<Int>::max();
  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, lower_pack_indx, upper_pack_indx+1),
                          [&] (const Int& k, Int& local_max) {
    auto indices_pack = ekat::range<IntSmallPack>(k*Spack::n);
    const auto in_range = (indices_pack < nlev-1 && indices_pack >= nlev-npbl);

    Spack vvk(0);
    vvk.set(in_range,
            ekat::max(tiny,
                      ekat::square(u(k) - s_u(nlev-1)) +
                      ekat::square(v(k) - s_v(nlev-1)) +
                      fac*(ustar*ustar)));

    if (in_range.any()) {
      rino(k).set(in_range,
                  ggr*(thv(k) - thv_ref)*(z(k) - s_z(nlev-1))/(s_thv(nlev-1)*vvk));
    }

    // Set indices_pack entry to -1 if rino(k)<ricr or
    // if global index is not in [nlev-npbl, nlev-1)
    indices_pack.set(rino(k)<ricr || !in_range, Kokkos::reduction_identity<Int>::max());

    // Local max index s.t. rino(k)>=ricr and k in [nlev-npbl, nlev-1)
    if (local_max < ekat::max(indices_pack))
      local_max = ekat::max(indices_pack);

  }, Kokkos::Max<Int>(max_indx));

  // Set check=false and compute pblh only if
  // there was an index s.t. rino(k)>=ricr.
  if (max_indx != Kokkos::reduction_identity<Int>::max()) {
    pblh = s_z(max_indx+1) +
          (ricr - s_rino(max_indx+1))/
          (s_rino(max_indx) - s_rino(max_indx+1))*
          (s_z(max_indx) - s_z(max_indx+1));
    check = false;
  }
}
} // namespace shoc
} // namespace scream

#endif
