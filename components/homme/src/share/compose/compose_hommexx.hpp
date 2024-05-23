#ifndef INCLUDE_COMPOSE_HOMMEXX_HPP
#define INCLUDE_COMPOSE_HOMMEXX_HPP

#include "compose_homme.hpp" // To get Real
#include <Kokkos_Core.hpp>

namespace homme {
namespace compose {

template <typename DataType>
using SetView = Kokkos::View<DataType, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>;

void set_views(const SetView<Real***>& spheremp,
               const SetView<Real****>& dp, const SetView<Real*****>& dp3d,
               const SetView<Real******>& qdp, const SetView<Real*****>& q,
               const SetView<Real*****>& dep_points);

void advect(const int np1, const int n0_qdp, const int np1_qdp);

void set_dp3d_np1(const int np1);
bool property_preserve_global();
bool property_preserve_local(const int limiter_option);
void property_preserve_check();

void finalize();

} // namespace compose
} // namespace homme

#endif
