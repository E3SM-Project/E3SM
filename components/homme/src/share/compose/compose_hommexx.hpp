#ifndef INCLUDE_COMPOSE_HOMMEXX_HPP
#define INCLUDE_COMPOSE_HOMMEXX_HPP

#include <Kokkos_Core.hpp>

namespace homme {
namespace compose {

typedef double HommexxReal;

template <typename DataType>
using SetView = Kokkos::View<DataType, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>;

typedef SetView<HommexxReal*****> SetView5;

void set_views(const SetView<HommexxReal***>& spheremp,
               const SetView<HommexxReal****>& dp, const SetView5& dp3d,
               const SetView<HommexxReal******>& qdp, const SetView5& q,
               const SetView5& dep_points, const SetView5& vnode,
               const SetView5& vdep, const int trajectory_ndim);

void set_hvcoord(const HommexxReal etai_beg, const HommexxReal etai_end,
                 const HommexxReal* etam);
void calc_v_departure(const int step, const HommexxReal dtsub);

void advect(const int np1, const int n0_qdp, const int np1_qdp);

void set_dp3d_np1(const int np1);
bool property_preserve_global();
bool property_preserve_local(const int limiter_option);
void property_preserve_check();

void finalize();

} // namespace compose
} // namespace homme

#endif
