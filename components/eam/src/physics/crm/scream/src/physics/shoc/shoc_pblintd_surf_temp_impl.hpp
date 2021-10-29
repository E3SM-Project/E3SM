#ifndef SHOC_PBLINTD_SURF_TEMP_IMPL_HPP
#define SHOC_PBLINTD_SURF_TEMP_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd_surf_temp. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd_surf_temp(const Int& nlev, const Int& nlevi, const Int& npbl,
      const uview_1d<const Spack>& z, const Scalar& ustar,
      const Scalar& obklen, const Scalar& kbfs,
      const uview_1d<const Spack>& thv, Scalar& tlv,
      Scalar& pblh, bool& check, const uview_1d<Spack>& rino)
{
  // const parameter for Diagnosis of PBL depth
  const Scalar fak   =  8.5;
  const Scalar betam = 15.0;
  const Scalar sffrac=  0.1;
  const Scalar binm  = betam*sffrac;

  // Scalarize views for single entry access
  const auto s_z    = ekat::scalarize(z);
  const auto s_thv  = ekat::scalarize(thv);
  const auto s_rino = ekat::scalarize(rino);

  // Estimate an effective surface temperature to account for surface
  // fluctuations
  if (check) {
    pblh = s_z(nlevi-npbl-1);
  }
  check = (kbfs > 0);

  tlv = s_thv(nlev-1);
  if (check) {
    s_rino(nlev-1) = 0;
    const Scalar phiminv = std::cbrt(1-binm*pblh/obklen);
    tlv += kbfs*fak/(ustar*phiminv);
  }
}

} // namespace shoc
} // namespace scream

#endif
