#include "share/util/scream_vertical_interpolation.hpp"

namespace scream {
namespace vinterp {

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
    default:
      EKAT_ERROR_MSG("vertical_remap::allocate_mask only supports a rank of 2 or 3, received rank = " + std::to_string(N));
  }
}

template<typename T, int P> 
KOKKOS_FUNCTION
void apply_interpolation_impl_1d(
  const view_1d<const Pack<T,P>>& x_src_col,
  const view_1d<const Pack<T,P>>& x_tgt_col,
  const view_1d<const Pack<T,P>>& input_col,
  const view_1d<      Pack<T,P>>& output_col,
  const view_1d<      Mask<P>>&   mask_col,
  const int nlevs_src,
  const int nlevs_tgt,
  const int icol,
  const T msk_val,
  const MemberType& team,
  const LIV<T,P>& vert_interp)
{
  // Recast source views to support different packsizes
  using PackInfo = ekat::PackInfo<P>;
  const int num_src_packs = PackInfo::num_packs(nlevs_src);
  const int num_tgt_packs = PackInfo::num_packs(nlevs_tgt);

  auto x_src  = Kokkos::subview(x_src_col,Kokkos::pair<int,int>(0,num_src_packs));
  auto input  = Kokkos::subview(input_col,Kokkos::pair<int,int>(0,num_src_packs));

  auto x_tgt  = Kokkos::subview(x_tgt_col, Kokkos::pair<int,int>(0,num_tgt_packs));
  auto output = Kokkos::subview(output_col,Kokkos::pair<int,int>(0,num_tgt_packs));
  auto mask   = Kokkos::subview(mask_col,  Kokkos::pair<int,int>(0,num_tgt_packs));

  // The input/output data and x_src/x_tgt data should match in the appropriate size, respectively.
  EKAT_KERNEL_REQUIRE_MSG(x_tgt.size() == output.size(), "Error! vertical_interpolation::apply_interpolation_imple_1d - target pressure level size does not match the size of the target data output.");
  EKAT_KERNEL_REQUIRE_MSG(x_src.size() == input.size() , "Error! vertical_interpolation::apply_interpolation_imple_1d - source pressure level size does not match the size of the source data input.");

  //Setup linear interpolation
  vert_interp.setup(team, x_src, x_tgt);
  //Run linear interpolation
  vert_interp.lin_interp(team, x_src, x_tgt, input, output, icol);
  const auto x_src_s = ekat::scalarize(x_src);
  const auto x_tgt_s = ekat::scalarize(x_tgt);
  const auto range = Kokkos::TeamVectorRange(team, x_tgt.extent(0));
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

/* ---------------------------------------------------------------------- 
 * Versions where x_tgt is a 2-D view
 * ---------------------------------------------------------------------- */
template<typename T, int P, int N> 
void perform_checks(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_2d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int                         nlevs_src,
  const int                         nlevs_tgt)
{
  auto rank = N;
  EKAT_REQUIRE_MSG (rank>1 &&rank<=3,"Error::scream_vertical_interpolation, passed view of rank (" + std::to_string(rank) +"), only support ranks 2 or 3\n");

  // The input data and x_src data should match in the appropriate size
  EKAT_REQUIRE(x_src.extent_int(0) == input.extent_int(0));
  EKAT_REQUIRE(x_src.extent_int(1) == input.extent_int(input.rank-1));
  // The output data and x_tgt data should match in the appropriate size
  EKAT_REQUIRE(x_tgt.extent_int(0) == output.extent_int(0));
  EKAT_REQUIRE(x_tgt.extent_int(1) == output.extent_int(input.rank-1));

  // The output data and the input data should match in all sizes except the last one
  for (int ii=0;ii<rank-1;ii++) {
    EKAT_REQUIRE(input.extent_int(ii)==output.extent_int(ii));
  }

  // The nlevs_src and nlevs_tgt values should at least be <= to the extent of x_src and x_tgt * packsize
  EKAT_REQUIRE(nlevs_src <= x_src.extent_int(1)*P);
  EKAT_REQUIRE(nlevs_tgt <= x_tgt.extent_int(1)*P);

}

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_2d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const view_Nd<      Mask<P>,N>&   mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  int ndofs = x_src.extent(0);
  for (int ii=1; ii<N-1; ii++) {
    ndofs *= input.extent_int(ii);
  }
  perform_checks<T,P,N>(x_src, x_tgt, input, output, nlevs_src, nlevs_tgt);
  LIV<T,P> vert_interp(ndofs,nlevs_src,nlevs_tgt);
  apply_interpolation(nlevs_src, nlevs_tgt, msk_val, vert_interp, x_src, x_tgt, input, output, mask);
}

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_2d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  int ndofs = x_src.extent(0);
  for (int ii=1; ii<N-1; ii++) {
    ndofs *= input.extent_int(ii);
  }
  perform_checks<T,P,N>(x_src, x_tgt, input, output, nlevs_src, nlevs_tgt);

  std::vector<int> extents;
  for (int ii=0;ii<output.rank;ii++) {
    extents.push_back(output.extent_int(ii));
  }
  const auto mask = allocate_mask<P,N>(extents);

  LIV<T,P> vert_interp(ndofs,nlevs_src,nlevs_tgt);
  for (int ii=1; ii<N-1; ii++) {
    ndofs *= input.extent_int(ii);
  }
  apply_interpolation(nlevs_src, nlevs_tgt, msk_val, vert_interp, x_src, x_tgt, input, output, mask);
}

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_2d<const Pack<T,P>>& x_tgt,
  const view_2d<const Pack<T,P>>& input,
  const view_2d<      Pack<T,P>>& output,
  const view_2d<        Mask<P>>& mask_out)
{
  const int d_0      = input.extent_int(0);
  const int npacks   = output.extent_int(output.rank-1);
  const auto policy = ESU::get_default_team_policy(d_0, npacks);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
        		 
    const int  icol  = team.league_rank();
    const auto x1   = ekat::subview(x_src,  icol);
    const auto xt   = ekat::subview(x_tgt,  icol);
    const auto in   = ekat::subview(input,  icol);
    const auto out  = ekat::subview(output, icol);
    const auto mask = ekat::subview(mask_out, icol);
    
    apply_interpolation_impl_1d<T,P>(x1,xt,in,out,mask,num_levs_src,num_levs_tgt,icol,mask_val,team,vert_interp);
  });
  Kokkos::fence();
}

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_2d<const Pack<T,P>>& x_tgt,
  const view_3d<const Pack<T,P>>& input,
  const view_3d<      Pack<T,P>>& output,
  const view_3d<        Mask<P>>& mask_out)
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
    const auto xt   = ekat::subview(x_tgt,  icol);
    const auto in   = ekat::subview(input,  icol, ivar);
    const auto out  = ekat::subview(output, icol, ivar);
    const auto mask = ekat::subview(mask_out, icol, ivar);

    apply_interpolation_impl_1d<T,P>(x1,xt,in,out,mask,num_levs_src,num_levs_tgt,icol,mask_val,team,vert_interp);
  });
  Kokkos::fence();   
}

/* ---------------------------------------------------------------------- 
 * Versions where x_tgt is a single 1-D vertical profile
 * ---------------------------------------------------------------------- */
template<typename T, int P, int N> 
void perform_checks(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int                         nlevs_src,
  const int                         nlevs_tgt)
{
  const int rank = input.rank;
  EKAT_REQUIRE_MSG (rank>1 &&rank<=3,"Error::scream_vertical_interpolation, passed view of rank (" + std::to_string(rank) +"), only support ranks 2 or 3\n");

  // The input data and x_src data should match in the appropriate size
  EKAT_REQUIRE(x_src.extent_int(0) == input.extent_int(0));


  // The output and input data should match in rank
  EKAT_REQUIRE(static_cast<int>(output.rank)==rank);
  // The output data and the input data should match in all sizes except the last one
  for (int ii=0;ii<rank-1;ii++) {
    EKAT_REQUIRE(input.extent_int(ii)==output.extent_int(ii));
  }

  // The nlevs_src and nlevs_tgt values should at least be <= to the extent of x_src and x_tgt * packsize
  EKAT_REQUIRE(nlevs_src <= x_src.extent_int(1)*P);
  EKAT_REQUIRE(nlevs_tgt <= x_tgt.extent_int(0)*P);

}

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const view_Nd<      Mask<P>,N>&   mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  int ndofs = x_src.extent(0);
  for (int ii=1; ii<N-1; ii++) {
    ndofs *= input.extent_int(ii);
  }
  perform_checks<T,P,N>(x_src, x_tgt, input, output, nlevs_src, nlevs_tgt);
  LIV<T,P> vert_interp(ndofs,nlevs_src,nlevs_tgt);
  apply_interpolation(nlevs_src, nlevs_tgt, msk_val, vert_interp, x_src, x_tgt, input, output, mask);
}

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val)
{
  int ndofs = x_src.extent(0);
  for (int ii=1; ii<N-1; ii++) {
    ndofs *= input.extent_int(ii);
  }
  perform_checks<T,P,N>(x_src, x_tgt, input, output, nlevs_src, nlevs_tgt);

  std::vector<int> extents;
  for (int ii=0;ii<output.rank;ii++) {
    extents.push_back(output.extent_int(ii));
  }
  const auto mask = allocate_mask<P,N>(extents);

  LIV<T,P> vert_interp(ndofs,nlevs_src,nlevs_tgt);
  apply_interpolation(nlevs_src, nlevs_tgt, msk_val, vert_interp, x_src, x_tgt, input, output, mask);
}

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_1d<const Pack<T,P>>& x_tgt,
  const view_2d<const Pack<T,P>>& input,
  const view_2d<      Pack<T,P>>& output,
  const view_2d<        Mask<P>>& mask_out)
{
  const int d_0      = input.extent_int(0);
  const int npacks   = output.extent_int(output.rank-1);
  const auto policy = ESU::get_default_team_policy(d_0, npacks);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
        		 
    const int  icol  = team.league_rank();
    const auto x1   = ekat::subview(x_src,  icol);
    const auto in   = ekat::subview(input,  icol);
    const auto out  = ekat::subview(output, icol);
    const auto mask = ekat::subview(mask_out, icol);
    
    apply_interpolation_impl_1d<T,P>(x1,x_tgt,in,out,mask,num_levs_src,num_levs_tgt,icol,mask_val,team,vert_interp);
  });
  Kokkos::fence();
}

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_1d<const Pack<T,P>>& x_tgt,
  const view_3d<const Pack<T,P>>& input,
  const view_3d<      Pack<T,P>>& output,
  const view_3d<        Mask<P>>& mask_out)
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
    const auto in   = ekat::subview(input,  icol, ivar);
    const auto out  = ekat::subview(output, icol, ivar);
    const auto mask = ekat::subview(mask_out, icol, ivar);

    apply_interpolation_impl_1d<T,P>(x1,x_tgt,in,out,mask,num_levs_src,num_levs_tgt,team.league_rank(),mask_val,team,vert_interp);
  });
  Kokkos::fence();   
}
  
} // namespace vinterp
} // namespace scream

