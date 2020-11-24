#ifndef SHOC_INTEG_COLUMN_STABILITY_IMPL_HPP
#define SHOC_INTEG_COLUMN_STABILITY_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::integ_column_stability(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& brunt,
  Scalar&                      brunt_int)
{
  static constexpr auto troppres = 80000;

  using ExeSpaceUtils = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  // Compute se_int
  ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {

    auto press_gt_troppress = (pres(k) > troppres);

    Spack my_result(0);
    my_result.set(press_gt_troppress, dz_zt(k) * brunt(k));

    return my_result ;
  }, brunt_int);

}

} // namespace shoc
} // namespace scream

#endif
