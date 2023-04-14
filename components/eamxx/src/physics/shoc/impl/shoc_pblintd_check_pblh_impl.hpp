#ifndef SHOC_PBLINTD_CHECK_PBLH_IMPL_HPP
#define SHOC_PBLINTD_CHECK_PBLH_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd_check_pblh. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd_check_pblh(const Int& nlevi, const Int& npbl, 
       const uview_1d<const Spack>& z, const Scalar& ustar, const bool& check, Scalar& pblh)
{
   // PBL height must be greater than some minimum mechanical mixing depth
   // Several investigators have proposed minimum mechanical mixing depth
   // relationships as a function of the local friction velocity, u*.  We
   // make use of a linear relationship of the form h = c u* where c=700.
   // The scaling arguments that give rise to this relationship most often
   // represent the coefficient c as some constant over the local coriolis
   // parameter.  Here we make use of the experimental results of Koracin
   // and Berkowicz (1988) [BLM, Vol 43] for which they recommend 0.07/f
   // where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
   // latitude value for f so that c = 0.07/f = 700.  Also, do not allow
   // PBL to exceed some maximum (npbl) number of allowable points
   //
   const auto s_z = ekat::scalarize(z);
   if (check) pblh = s_z(nlevi-npbl-1);
   pblh = ekat::impl::max(pblh, 700*ustar);
}

} // namespace shoc
} // namespace scream

#endif
