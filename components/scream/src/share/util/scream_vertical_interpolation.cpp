#include "share/util/scream_vertical_interpolation.hpp"

namespace scream {

  
void perform_vertical_interpolation(const view_2d<const Spack>& x_src,
				    const view_1d<const Spack>& x_tgt,
				    const view_2d<const Spack>& input,
				    const view_2d<      Spack>& output,
				    const int                 & nlevs_src,
				    const int                 & nlevs_tgt
				    )
{
const view_2d<Smask> mask("",x_src.extent(0),x_tgt.extent(0));
perform_vertical_interpolation(x_src, x_tgt, input, output, mask,
			       nlevs_src, nlevs_tgt, masked_val);
}

void perform_vertical_interpolation(const view_2d<const Spack>& x_src,
				    const view_1d<const Spack>& x_tgt,
				    const view_2d<const Spack>& input,
				    const view_2d<      Spack>& output,
				    const view_2d<      Smask>& mask,
				    const int                 & nlevs_src,
				    const int                 & nlevs_tgt
				    )
{
perform_vertical_interpolation(x_src, x_tgt, input, output, mask,
			       nlevs_src, nlevs_tgt, masked_val);
}
  
void perform_vertical_interpolation(const view_2d<const Spack>& x_src,
				    const view_1d<const Spack>& x_tgt,
				    const view_2d<const Spack>& input,
				    const view_2d<      Spack>& output,
				    const view_2d<      Smask>& mask,
				    const int                 & nlevs_src,
				    const int                 & nlevs_tgt,
				    const Real                & masked_val
				    )
{
  const int ncols = x_src.extent(0);
  //Do a bunch of checks to make sure that Spacks are consistent
  auto npacks_src = ekat::PackInfo<Spack::n>::num_packs(nlevs_src);
  auto npacks_tgt = ekat::PackInfo<Spack::n>::num_packs(nlevs_tgt);
  EKAT_REQUIRE(x_src.extent(0)==input.extent(0));
  EKAT_REQUIRE(x_src.extent(1)==input.extent(1));
  EKAT_REQUIRE(x_src.extent(1)==npacks_src);
  EKAT_REQUIRE(input.extent(1)==npacks_src);
  EKAT_REQUIRE(x_tgt.extent(0)==output.extent(1));
  EKAT_REQUIRE(x_tgt.extent(0)==npacks_tgt); 
  
  LIV vert_interp(ncols,nlevs_src,nlevs_tgt);

  const int num_vert_packs = x_tgt.extent(0);
  const auto policy = ESU::get_default_team_policy(ncols, num_vert_packs);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
			 
    const int icol = team.league_rank();

    auto x1=ekat::subview(x_src, icol);
    auto in=ekat::subview(input, icol);
    auto out=ekat::subview(output, icol);
    auto msk=ekat::subview(mask, icol);
    
    perform_vertical_interpolation(x1, x_tgt, in, out, msk, nlevs_src, nlevs_tgt,
				   icol, masked_val, team, vert_interp);

  });
  Kokkos::fence();   
}

void perform_vertical_interpolation(const view_1d<const Spack>& x_src,
				    const view_1d<const Spack>& x_tgt,
				    const view_1d<const Spack>& input,
				    const view_1d<      Spack>& output,
				    const view_1d<      Smask>& mask,
				    const int                 & nlevs_src,
				    const int                 & nlevs_tgt,
				    const int                 & icol,
				    const Real                & masked_val,
				    const LIV::MemberType     & team,
				    const LIV                 & vert_interp
				    )
{
  //Setup linear interpolation
  vert_interp.setup(team, x_src, x_tgt);
  //Run linear interpolation
  vert_interp.lin_interp(team, x_src, x_tgt, input, output, icol);
  const auto x_src_s = ekat::scalarize(x_src);
  const auto x_tgt_s = ekat::scalarize(x_tgt);
  const auto range_boundary = KT::RangePolicy(0, x_tgt.extent(0));
  //Mask out values above (below) maximum (minimum) source grid
  Kokkos::parallel_for(range_boundary, [&] (const int & k) {
    const auto above_max = x_tgt[k] > x_src_s[nlevs_src-1];
    const auto below_min = x_tgt[k] < x_src_s[0];
    const auto combined_m = above_max || below_min;
    mask(k) = combined_m;
    output(k).set(combined_m,masked_val);
  });
  Kokkos::fence();
}

} // namespace scream
