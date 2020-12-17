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
  const Scalar onet  = 1./3.;
  const Scalar fak   =  8.5;
  const Scalar betam = 15.0;
  const Scalar sffrac=  0.1;
  const Scalar binm  = betam*sffrac;

  //
  // Estimate an effective surface temperature to account for surface
  // fluctuations
  //
  const auto zs    = scalarize(z);
  const auto thvs  = scalarize(thv);
  const auto rinos = scalarize(rino);

  if (check) pblh = zs(nlevi-npbl-1);
  check  = kbfs > 0.;
  tlv    = thvs(nlev-1);
  if (check) {
     const auto phiminv = std::pow((1 - binm*pblh/obklen), onet);
     rinos(nlev-1) = 0.0;
     tlv  = thvs(nlev-1) + kbfs*fak/(ustar*phiminv);
  }
}

} // namespace shoc
} // namespace scream

#endif
