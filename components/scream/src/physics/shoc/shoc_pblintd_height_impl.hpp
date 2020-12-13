#ifndef SHOC_PBLINTD_HEIGHT_IMPL_HPP
#define SHOC_PBLINTD_HEIGHT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd_height. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd_height(const MemberType& team, const Int& nlev, const Int& npbl, 
        const uview_1d<const Spack>& z, const uview_1d<const Spack>& u, 
        const uview_1d<const Spack>& v, const Scalar& ustar, const uview_1d<const Spack>& thv, 
        const Scalar& thv_ref, Scalar& pblh, const uview_1d<Spack>& rino, bool& check)
{
  //
  // PBL height calculation:  Scan upward until the Richardson number between
  // the first level and the current level exceeds the "critical" value.
  //
  const auto tiny = 1.e-36;
  const auto fac  = 100.;
  const auto ricr  =  0.3;
  const auto ggr = C::gravit;

  const auto sz = scalarize(z);
  
  Spack z_plus, z_mid, rino_plus, rino_mid;
  Spack vvk, pblh_sp, rino_sp;

  Int index = 0;
  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k, Int& local_max) {

    const Int i = team.league_rank();

    local_max = 0;
    auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
    auto range_pack2 = range_pack1;
    range_pack2.set(range_pack1 > nlev, nlev); 

    const auto rmask = range_pack1 < nlev-1;
    const auto lmask = range_pack1 >= nlev-npbl;

    auto index_max_pack = ekat::range<IntSmallPack>(k*Spack::n);
    index_max_pack.set(!(lmask && rmask), 0);

    const int nlev_indx  = (nlev-1)%Spack::n;
    const auto u_nlev    = u(nlev_pack-1)[nlev_indx];
    const auto v_nlev    = v(nlev_pack-1)[nlev_indx];
    const auto thv_nlev  = thv(nlev_pack-1)[nlev_indx];
    const auto z_nlev    = z(nlev_pack-1)[nlev_indx];

    if ((lmask && rmask).any()) {
      vvk.set(lmask && rmask, 
              pow((u(k) - u_nlev), 2) + pow((v(k) - v_nlev),2) + fac*pow(ustar,2));
      vvk.set(lmask && rmask,  ekat::max(vvk, tiny));

      rino_sp.set(lmask && rmask,  
                  ggr*(thv(k)-thv_ref)*(z(k)-z_nlev)/(thv_nlev*vvk));

      const auto rino_lt_ricr = rino_sp < ricr;

      index_max_pack.set(rino_lt_ricr, 0);
      const Int max_indx = ekat::max(index_max_pack);

      if (max_indx > local_max) {
         local_max = max_indx;
      }
    }
  }, Kokkos::Max<int>(index));

  team.team_barrier();

  // calculate the pblh
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {
     auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
     auto range_pack2 = range_pack1;
     range_pack2.set(range_pack1 > nlev, nlev);

     const auto lmask = range_pack2 < nlev-1;
     const auto rmask = range_pack2 >= nlev-npbl;
     auto max_index = (range_pack2 >= index) && lmask && rmask;

     if (max_index.any() && check) {
       const int nlev_indx  = (nlev-1)%Spack::n;
       const auto u_nlev    = u(nlev_pack-1)[nlev_indx];
       const auto v_nlev    = v(nlev_pack-1)[nlev_indx];
       const auto thv_nlev  = thv(nlev_pack-1)[nlev_indx];
       const auto z_nlev    = z(nlev_pack-1)[nlev_indx];

       vvk.set(max_index,
               pow((u(k) - u_nlev), 2) + pow((v(k) - v_nlev),2) + fac*pow(ustar,2));
       vvk.set(max_index,  ekat::max(vvk, tiny));

       rino(k).set(max_index,
                   ggr*(thv(k)-thv_ref)*(z(k)-z_nlev)/(thv_nlev*vvk));

       auto srino = scalarize(rino);
       ekat::index_and_shift<1>(sz, range_pack2, z_mid, z_plus);
       ekat::index_and_shift<1>(srino, range_pack2, rino_mid, rino_plus);
       pblh_sp.set(max_index, z_plus+(ricr-rino_plus)/(rino_mid-rino_plus)*(z_mid-z_plus));

       if (max_index.any() && (k == index/Spack::n)) {
         auto pblh_indx = index - k*Spack::n;
         const auto rino_ge_ricr = rino(k) >= ricr;
         if (rino_ge_ricr[pblh_indx]) {
           pblh = pblh_sp[pblh_indx];
           check = false;
         }
       }
     }
  });
}
} // namespace shoc
} // namespace scream

#endif
