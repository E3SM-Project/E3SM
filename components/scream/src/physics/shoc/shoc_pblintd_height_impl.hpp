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

  const auto nlev_v = (nlev-1)/Spack::n;
  const auto nlev_p = (nlev-1)%Spack::n;

  // Compute rino values and find max index s.t. rino(k) >= ricr
  Int max_indx;

  const Int lower_pack_indx = (nlev-npbl)/Spack::n;
  const Int upper_pack_indx = (nlev-1)/Spack::n+1;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, lower_pack_indx, upper_pack_indx),
                          [&] (const Int& k, Int& local_max) {
    auto indices_pack = ekat::range<IntSmallPack>(k*Spack::n);
    const auto in_range = (indices_pack < nlev-1 && indices_pack >= nlev-npbl);

    Spack vvk(0);
    vvk.set(in_range,
            ekat::max(tiny,
                      ekat::square(u(k) - u(nlev_v)[nlev_p]) +
                      ekat::square(v(k) - v(nlev_v)[nlev_p]) +
                      fac*(ustar*ustar)));
    rino(k).set(in_range,
                ggr*(thv(k) - thv_ref)*(z(k) - z(nlev_v)[nlev_p])/(thv(nlev_v)[nlev_p]*vvk));

    // Set indices_pack entry to -1 if rino(k)<ricr or
    // if global index is not in [nlev-npbl, nlev-1)
    indices_pack.set(rino(k)<ricr || !in_range, -1);

    // Local max index s.t. rino(k)>=ricr and k in [nlev-npbl, nlev-1)
    if (local_max < ekat::max(indices_pack))
      local_max = ekat::max(indices_pack);

  }, Kokkos::Max<Int>(max_indx));

  // Set check=false and compute pblh only if
  // there was an index s.t. rino(k)>=ricr.
  // If no index was found, set max_index=nlev-npbl.
  if (max_indx != -1) {
    const auto max_indx_v = max_indx/Spack::n;
    const auto max_indx_p = max_indx%Spack::n;
    const auto max_indx_p1_v = (max_indx+1)/Spack::n;
    const auto max_indx_p1_p = (max_indx+1)%Spack::n;

    pblh = z(max_indx_p1_v)[max_indx_p1_p] +
           (ricr - rino(max_indx_p1_v)[max_indx_p1_p])/
           (rino(max_indx_v)[max_indx_p] - rino(max_indx_p1_v)[max_indx_p1_p])*
           (z(max_indx_v)[max_indx_p] - z(max_indx_p1_v)[max_indx_p1_p]);

    check = false;
  } else {
    max_indx = nlev-npbl;
  }
}
} // namespace shoc
} // namespace scream

#endif
