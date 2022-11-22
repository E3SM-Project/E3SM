#ifndef SHOC_DIAG_OBKLEN_IMPL_HPP
#define SHOC_DIAG_OBKLEN_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_diag_obklen. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_diag_obklen(
  const Scalar& uw_sfc,
  const Scalar& vw_sfc,
  const Scalar& wthl_sfc,
  const Scalar& wqw_sfc,
  const Scalar& thl_sfc,
  const Scalar& cldliq_sfc,
  const Scalar& qv_sfc,
  Scalar&       ustar,
  Scalar&       kbfs,
  Scalar&       obklen)
{
  // Constants
  const auto cp = C::CP;
  const auto lcond = C::LatVap;
  const auto eps = C::ZVIR;
  const auto ustar_min = SC::ustar_min;
  const auto ggr = C::gravit;
  const auto vk = C::Karman;

  // Local variables
  const Scalar th_sfc = thl_sfc + (lcond/cp)*cldliq_sfc;
  const Scalar thv_sfc = th_sfc*(1 + eps*qv_sfc - cldliq_sfc);
  const Scalar ustar_val = std::sqrt(uw_sfc*uw_sfc + vw_sfc*vw_sfc);

  // Outputs
  ustar = ekat::impl::max(ustar_min, ustar_val);
  kbfs = wthl_sfc + eps*th_sfc*wqw_sfc;
  const Scalar sign_val = (kbfs>=0 ? 1e-10 : -1e-10);
  obklen = -thv_sfc*(ustar*ustar*ustar)/(ggr*vk*(kbfs + sign_val));
}

} // namespace shoc
} // namespace scream

#endif
