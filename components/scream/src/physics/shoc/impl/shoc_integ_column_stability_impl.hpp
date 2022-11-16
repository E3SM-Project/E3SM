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
  //Lower troposphere pressure [Pa]
  static constexpr Scalar troppres = 80000;

  using ExeSpaceUtils = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  brunt_int = ExeSpaceUtils::view_reduction(team,0,nlev,
                                            [&] (const int k) -> Spack {

    //calculate only when pressure is > troposphere pressure
    auto press_gt_troppress = (pres(k) > troppres);

    Spack return_val(0); //initialize return value for brunt_int
    return_val.set(press_gt_troppress, dz_zt(k) * brunt(k));// compute brunt_int for each column

    return return_val ;

    });
}

} // namespace shoc
} // namespace scream

#endif
