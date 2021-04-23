#ifndef COMMON_PHYSICS_FUNCTIONS_HPP
#define COMMON_PHYSICS_FUNCTIONS_HPP

#include "physics_constants.hpp"

namespace scream {
namespace physics {

template <typename DeviceT>
struct PhysicsFunctions
{

  using Device = DeviceT;
  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;

  //--------------------------------------------------------------------//
  //  Functions
  //--------------------------------------------------------------------//

  template<typename ScalarT>
  KOKKOS_FUNCTION
  static ScalarT get_exner(const ScalarT& p_mid);

  template<typename ScalarT, typename InputProviderP>
  KOKKOS_FUNCTION
  static void get_exner(const MemberType& team,
                        const InputProviderP& p_mid, 
                        const view_1d<ScalarT>& exner);

}; // struct PhysicsFunctions

} // namespace physics
} // namespace scream

#endif // COMMON_PHYSICS_FUNCTIONS_HPP
