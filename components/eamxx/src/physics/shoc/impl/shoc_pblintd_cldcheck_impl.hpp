
#ifndef SHOC_PBLINTD_CLDCHECK_IMPL_HPP
#define SHOC_PBLINTD_CLDCHECK_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_pblintd_cldcheck(
    const Scalar& zi, const Scalar& cldn, 
    Scalar& pblh)
{
    //
    // Final requirement on PBL heightis that it must be greater than the depth
    // of the lowest model level if there is any cloud diagnosed in
    // the lowest model level.  This is to deal with the inadequacies of the
    // current "dry" formulation of the boundary layer, where this test is
    // used to identify circumstances where there is marine stratus in the
    // lowest level, and to provide a weak ventilation of the layer to avoid
    // a pathology in the cloud scheme (locking in low-level stratiform cloud)
    // If  any cloud is diagnosed in the lowest level, set pblh to 50 meters 
    // higher than top interface of lowest level
    //

    if (cldn >= 0) pblh = ekat::impl::max<Scalar>(pblh, zi + sp(50));
}
} // namespace shoc
} // namespace scream

#endif
