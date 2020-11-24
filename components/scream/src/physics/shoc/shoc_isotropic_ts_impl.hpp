#ifndef SHOC_ISOTROPIC_TS_IMPL_HPP
#define SHOC_ISOTROPIC_TS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc isotropic_ts. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::isotropic_ts(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   shcol,
  const uview_1d<const Spack>& brunt_int,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& a_diss,
  const uview_1d<const Spack>& brunt,
  const uview_1d<Spack>&       isotropy)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
