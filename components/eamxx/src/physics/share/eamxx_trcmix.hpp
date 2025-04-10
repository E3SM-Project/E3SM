#ifndef SCREAM_TRCMIX_HPP
#define SCREAM_TRCMIX_HPP

#include "share/eamxx_types.hpp"

//
// Bridge function to call trcmix from CXX
//

namespace scream {
namespace physics {

template<typename Scalar>
using trcmix_view1d = typename ekat::KokkosTypes<DefaultDevice>::template view_1d<Scalar>;

template<typename Scalar>
using trcmix_view2d = typename ekat::KokkosTypes<DefaultDevice>::template view_2d<Scalar>;

// NOTE: nlevs be smaller than pmid/q extent, in case one/both of them are padded
void trcmix(
  const std::string& name,      // constituent name
  const int nlevs,              // number of physical levels
  trcmix_view1d<const Real> const& clat,  // latitude for columns in degrees
  trcmix_view2d<const Real> const& pmid,  // model pressures
  trcmix_view2d<Real>            & q,     // constituent mass mixing ratio (output)
  // items below likely come from namelists/yaml
  const Real co2vmr,     // co2 volume mixing ratio
  const Real n2ovmr,     // n2o volume mixing ratio
  const Real ch4vmr,     // ch4 volume mixing ratio
  const Real f11vmr,     // cfc11 volume mixing ratio
  const Real f12vmr);    // cfc12 volume mixing ratio

} // namespace physics
} // namespace scream

#endif
