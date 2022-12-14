#include "share/util/scream_vertical_interpolation.hpp"

namespace scream {
namespace vinterp {

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,N>& x_src,
  const view_1d<Pack<T,P>>&   x_tgt,
  const view_Nd<Pack<T,P>,N>& input,
  const view_Nd<Pack<T,P>,N>& output,
  const view_Nd<Mask<P>,N>&   mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  perform_vertical_interpolation_impl_Nd<Real,P,N>(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, msk_val);
}

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,N>& x_src,
  const view_1d<Pack<T,P>>&   x_tgt,
  const view_Nd<Pack<T,P>,N>& input,
  const view_Nd<Pack<T,P>,N>& output,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  std::vector<int> extents;
  for (int ii=0;ii<output.rank;ii++) {
    extents.pushback(output.extent_int(ii));
  }
  const auto mask = allocate_mask(extents);
  perform_vertical_interpolation_impl_Nd<Real,P,N>(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, msk_val);
}

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,N>& x_src,
  const view_1d<Pack<T,P>>&   x_tgt,
  const view_Nd<Pack<T,P>,N>& input,
  const view_Nd<Pack<T,P>,N>& output,
  const view_Nd<Mask<P>,N>&   mask,
  const int nlevs_src,
  const int nlevs_tgt)
{
  perform_vertical_interpolation_impl_Nd<Real,P,N>(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, masked_val);
}

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,N>& x_src,
  const view_1d<Pack<T,P>>&   x_tgt,
  const view_Nd<Pack<T,P>,N>& input,
  const view_Nd<Pack<T,P>,N>& output,
  const int nlevs_src,
  const int nlevs_tgt)
{
  std::vector<int> extents;
  for (int ii=0;ii<output.rank;ii++) {
    extents.pushback(output.extent_int(ii));
  }
  const auto mask = allocate_mask(extents);
  perform_vertical_interpolation_impl_Nd<Real,P,N>(x_src, x_tgt, input, output, mask,
                                         nlevs_src, nlevs_tgt, masked_val);
}

template<int P, int N>
view_Nd<Mask<P>,N> allocate_mask(const std::vector<int>& extents)
{
  switch(extents.size()) {
    case 1:
      return view_Nd<Mask<P>,N>("",extents[0]);
    case 2:
      return view_Nd<Mask<P>,N>("",extents[0],extents[1]);
    case 3:
      return view_Nd<Mask<P>,N>("",extents[0],extents[1],extents[2]);
    case 4:
      return view_Nd<Mask<P>,N>("",extents[0],extents[1],extents[2],extents[3]);
    default:
      EKAT_ERROR_MSG("vertical_remap::allocate_mask only supports a rank of <= 4, received rank = " + std::to_string(N));
  }
}

template<typename T, int P> 
KOKKOS_FUNCTION
void perform_vertical_interpolation_impl_1d(
  const view_1d<Pack<T,P>>& x_src,
  const view_1d<Pack<T,P>>& x_tgt,
  const view_1d<Pack<T,P>>& input,
  const view_1d<Pack<T,P>>& output,
  const view_1d<Mask<P>>& mask,
  const int nlevs_src,
  const int icol,
  const Real msk_val,
  const MemberType& team,
  const LIV<T,P>& vert_interp)
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

template<typename T, int P, int N> 
void perform_vertical_interpolation_impl_Nd(
  const view_Nd<const Pack<T,P>,2>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<Pack<T,P>,N>&       output,
  const view_Nd<Mask<P>,N>&         mask,
  const int                         nlevs_src,
  const int                         nlevs_tgt,
  const Real                        msk_val)
{
  auto rank = input.rank;
  EKAT_REQUIRE_MSG (rank<=4,"Error::scream_vertical_interpolation, passed view of rank (" + std::to_string(rank) +"), only support rank <= 4\n");
  const int ncols = rank > 1 ? x_src.extent(0) : 1;
  auto npacks_src = ekat::PackInfo<Pack<T,P>::n>::num_packs(nlevs_src);
  auto npacks_tgt = ekat::PackInfo<Pack<T,P>::n>::num_packs(nlevs_tgt);
  // Do a series of checks to make sure that Pack<T,P>s are consistent:
  // Check 1: That the source and input share the same size
  EKAT_REQUIRE(x_src.extent(0)==input.extent(0));
  if (rank>1) {
    EKAT_REQUIRE(x_src.extent(1)==input.extent(rank-1));
  }
  // Check 2: That unpacked there is more data than the expected number of levels.
  EKAT_REQUIRE(x_src.extent_int(1)*P >= nlevs_src);
  // Check 3: That the size of the output matches the target interpolation extent.
  EKAT_REQUIRE(x_tgt.extent(0)==output.extent(1));
  EKAT_REQUIRE(x_tgt.extent_int(0)==npacks_tgt);
  // Check 4: Source and target views should be of the same rank
  EKAT_REQUIRE(rank==output.rank);
  // Check 5: If rank=1, then the source x space needs to also be rank 1
  if (rank==1) {
    EKAT_REQUIRE(x_tgt.rank==1);
  }
  
  LIV<T,P> vert_interp(ncols,nlevs_src,nlevs_tgt);
  apply_interpolation(nlevs_src, msk_val, vert_interp, x_src, x_tgt, input, output, mask);
}

template<typename T, int P> 
void apply_interpolation(
  const                        int  num_levs,
  const                          T  mask_val,
  const                   LIV<T,P>& vert_interp,
  const view_Nd<const Pack<T,P>,2>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,2>& input,
  const view_Nd<      Pack<T,P>,2>& output,
  const view_Nd<        Mask<P>,2>& mask_out)
{
  const int d_0      = input.extent_int(0);
  const int npacks   = output.extent_int(output.rank-1);
  const auto policy = ESU::get_default_team_policy(d_0, npacks);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
        		 
    const int icol  = team.league_rank();
    const auto x1   = ekat::subview(x_src,  icol);
    const auto x1_s = ekat::scalarize(x1);
    const auto in   = ekat::subview(input,  icol);
    const auto out  = ekat::subview(output, icol);
    const auto mask = ekat::subview(mask_out, icol);
    
    vert_interp.setup(team, x1, x_tgt);
    vert_interp.lin_interp(team, x1, x_tgt, in, out, icol);

    const int ivec     = num_levs % P;
    apply_masking<T,P>(team, mask_val, x1_s[0], x1_s[num_levs-1], x_tgt, out, mask);   
    team.team_barrier();
  });
  Kokkos::fence();
}

template<typename T, int P> 
void apply_interpolation(
  const                        int  num_levs,
  const                          T  mask_val,
  const                   LIV<T,P>& vert_interp,
  const view_Nd<const Pack<T,P>,2>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,3>& input,
  const view_Nd<      Pack<T,P>,3>& output,
  const view_Nd<        Mask<P>,3>& mask_out)
{
  const int d_0      = input.extent_int(0);
  const int num_vars = input.extent_int(1);
  const int npacks   = output.extent_int(output.rank-1);
  const auto policy = ESU::get_default_team_policy(d_0*num_vars, npacks);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
        		 
    const int icol  = team.league_rank() / num_vars;
    const int ivar  = team.league_rank() % num_vars;
    const auto x1   = ekat::subview(x_src,  icol);
    const auto x1_s = ekat::scalarize(x1);
    const auto in   = ekat::subview(input,  icol, ivar);
    const auto out  = ekat::subview(output, icol, ivar);
    const auto mask = ekat::subview(mask_out, icol, ivar);
    
    vert_interp.setup(team, x1, x_tgt);
    vert_interp.lin_interp(team, x1, x_tgt, in, out, icol);

    const int ivec     = num_levs % P;
    apply_masking<T,P>(team, mask_val, x1_s[0], x1_s[num_levs-1], x_tgt, out, mask);   
    team.team_barrier();
  });
  Kokkos::fence();   
}

template<typename T, int P> 
void apply_interpolation(
  const                        int  num_levs,
  const                          T  mask_val,
  const                   LIV<T,P>& vert_interp,
  const view_Nd<const Pack<T,P>,2>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,4>& input,
  const view_Nd<      Pack<T,P>,4>& output,
  const view_Nd<        Mask<P>,4>& mask_out)
{
  const int d_0 = input.extent_int(0);
  const int d_1 = input.extent_int(1);
  const int d_2 = input.extent_int(2);
  const int npacks = output.extent_int(output.rank-1);
  const int num_vars   = d_1*d_2;
  const auto policy = ESU::get_default_team_policy(d_0*d_1*d_2, npacks);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
        		 
    const int icol   = team.league_rank() / num_vars;
    const int ivar   = team.league_rank() % num_vars;
    const int islc_1 = ivar / d_2;
    const int islc_2 = ivar % d_2;
    const auto x1    = ekat::subview(x_src,  icol);
    const auto x1_s  = ekat::scalarize(x1);
    const auto in    = ekat::subview(input,  icol, islc_1, islc_2);
    const auto out   = ekat::subview(output, icol, islc_1, islc_2);
    const auto mask  = ekat::subview(mask_out, icol, islc_1, islc_2);
    
    vert_interp.setup(team, x1, x_tgt);
    vert_interp.lin_interp(team, x1, x_tgt, in, out, icol);

    const int ivec     = num_levs % P;
    apply_masking<T,P>(team, mask_val, x1_s[0], x1_s[num_levs-1], x_tgt, out, mask);   
    team.team_barrier();
  });
  Kokkos::fence();   
}
  
template<typename T, int P>
KOKKOS_FUNCTION
void apply_masking(
  const MemberType& team,
  const T           mask_val,
  const T           min_val,
  const T           max_val,
  const view_1d<const Pack<T,P>>& x_tgt,
  const view_1d<      Pack<T,P>>& out,
  const view_1d<      Mask<P>>&   mask)
{
  const auto range = Kokkos::TeamThreadRange(team, x_tgt.extent(0));
  Kokkos::parallel_for(range, [&] (const Int & k) {
    const auto above_max = x_tgt(k) > max_val;
    const auto below_min = x_tgt(k) < min_val;
    const auto combined_m = above_max || below_min;
    mask(k) = combined_m;
    out(k).set(combined_m,mask_val);
  });
  team.team_barrier();
}

} // namespace vinterp
} // namespace scream

