#ifndef SHOC_COMPUTE_CONV_VEL_SHOC_LENGTH_HPP
#define SHOC_COMPUTE_CONV_VEL_SHOC_LENGTH_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::compute_conv_vel_shoc_length(
  const MemberType&            team,
  const Int&                   nlev,
  const Scalar&                pblh,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& thv,
  const uview_1d<const Spack>& wthv_sec,
  Scalar&                      conv_vel)
{
  const auto ggr = C::gravit;

  // Loop over levels runs inverted, i.e., k=nlev-1,...,0
  // We invert the packs for BFB matching which changes the
  // indexing of the view reduction.
  const int inv_begin = (nlev%Spack::n == 0 ? 0 : Spack::n - nlev%Spack::n);
  const int inv_end   = nlev + inv_begin;
  ekat::ExeSpaceUtils<typename KT::ExeSpace>::view_reduction(team,inv_begin,inv_end,
                                                             [&] (const int k) -> Spack {
      Spack
        return_pack(0),
        inv_zt_grid,
        inv_dz_zt,
        inv_thv,
        inv_wthv_sec;

     // Consider k-index in reverse order. For each k-index,
     // flip the pack.
     const int inv_k_indx = ekat::npack<Spack>(nlev) - (k+1);
     for (int p=0; p<Spack::n; ++p) {
       const int inv_p_indx = Spack::n-(p+1);

       inv_zt_grid[p] = zt_grid(inv_k_indx)[inv_p_indx];
       inv_dz_zt[p] = dz_zt(inv_k_indx)[inv_p_indx];
       inv_thv[p] = thv(inv_k_indx)[inv_p_indx];
       inv_wthv_sec[p] = wthv_sec(inv_k_indx)[inv_p_indx];
     }

    return_pack.set(inv_zt_grid < pblh,
                    sp(2.5)*inv_dz_zt*(ggr/inv_thv)*inv_wthv_sec);

    return return_pack;
  }, conv_vel);
}

} // namespace shoc
} // namespace scream

#endif
