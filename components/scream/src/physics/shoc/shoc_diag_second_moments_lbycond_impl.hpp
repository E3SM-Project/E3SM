#ifndef SHOC_DIAG_SECOND_MOMENTS_LBYCOND_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_LBYCOND_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc diag_second_moments_lbycond. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::diag_second_moments_lbycond(const Int& shcol, const uview_1d<const Spack>& wthl_sfc, const uview_1d<const Spack>& wqw_sfc, const uview_1d<const Spack>& uw_sfc, const uview_1d<const Spack>& vw_sfc, const uview_1d<const Spack>& ustar2, const uview_1d<const Spack>& wstar, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& thl_sec, const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& qwthl_sec)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace shoc
} // namespace scream

#endif
