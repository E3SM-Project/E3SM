#include "share/util/scream_vertical_interpolation.hpp"

namespace scream {
namespace vinterp {

template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  perform_vertical_interpolation_impl_2d(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, msk_val);
}

template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  const view_2d<Mask<N>> mask("",x_src.extent(0),x_tgt.extent(0));
  perform_vertical_interpolation_impl_2d(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, msk_val);
}

template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt)
{
  perform_vertical_interpolation_impl_2d(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, masked_val);
}

template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const int nlevs_src,
  const int nlevs_tgt)
{
  const view_2d<Mask<N>> mask("",x_src.extent(0),x_tgt.extent(0));
  perform_vertical_interpolation_impl_2d(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, masked_val);
}

template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation_impl_2d(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  const int ncols = x_src.extent(0);
  //Do a bunch of checks to make sure that Pack<T,N>s are consistent
  auto npacks_src = ekat::PackInfo<Pack<T,N>::n>::num_packs(nlevs_src);
  auto npacks_tgt = ekat::PackInfo<Pack<T,N>::n>::num_packs(nlevs_tgt);
  EKAT_REQUIRE(x_src.extent(0)==input.extent(0));
  EKAT_REQUIRE(x_src.extent(1)==input.extent(1));
  EKAT_REQUIRE(x_src.extent_int(1)*N >= nlevs_src);
  EKAT_REQUIRE(x_tgt.extent(0)==output.extent(1));
  EKAT_REQUIRE(x_tgt.extent_int(0)==npacks_tgt);
  
  LIV<T,N> vert_interp(ncols,nlevs_src,nlevs_tgt);

  const int num_vert_packs = x_tgt.extent(0);
  const auto policy = ESU::get_default_team_policy(ncols, num_vert_packs);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
			 
    const int icol = team.league_rank();
    auto x1=ekat::subview(x_src, icol);
    auto in=ekat::subview(input, icol);
    auto out=ekat::subview(output, icol);
    auto msk=ekat::subview(mask, icol);
    
    vert_interp.setup(team, x1, x_tgt);
    vert_interp.lin_interp(team, x1, x_tgt, in, out, icol);
    const auto x_src_s = ekat::scalarize(x1);
    const auto x_tgt_s = ekat::scalarize(x_tgt);
    const auto range = Kokkos::TeamThreadRange(team, x_tgt.extent(0));
    //Mask out values above (below) maximum (minimum) source grid
    Kokkos::parallel_for(range, [&] (const Int & k) {
      const auto above_max = x_tgt[k] > x_src_s[nlevs_src-1];
      const auto below_min = x_tgt[k] < x_src_s[0];
      const auto combined_m = above_max || below_min;
      msk(k) = combined_m;
      out(k).set(combined_m,msk_val);
    });
    team.team_barrier();
  });
  Kokkos::fence();   
}

template<typename Src, typename Tgt, typename Input, typename T, int N> 
KOKKOS_FUNCTION
void perform_vertical_interpolation_impl_1d(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_1d<Pack<T,N>>& output,
  const view_1d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const int icol,
  const Real msk_val,
  const MemberType& team,
  const LIV<T,N>& vert_interp)
{
  //Setup linear interpolation
  vert_interp.setup(team, x_src, x_tgt);
  //Run linear interpolation
  vert_interp.lin_interp(team, x_src, x_tgt, input, output, icol);
  const auto x_src_s = ekat::scalarize(x_src);
  const auto x_tgt_s = ekat::scalarize(x_tgt);
  const auto range = Kokkos::TeamThreadRange(team, x_tgt.extent(0));
  //Mask out values above (below) maximum (minimum) source grid
  Kokkos::parallel_for(range, [&] (const Int & k) {
    const auto above_max = x_tgt[k] > x_src_s[nlevs_src-1];
    const auto below_min = x_tgt[k] < x_src_s[0];
    const auto combined_m = above_max || below_min;
    mask(k) = combined_m;
    output(k).set(combined_m,msk_val);
  });
  team.team_barrier();
}
  
} // namespace vinterp
} // namespace scream

