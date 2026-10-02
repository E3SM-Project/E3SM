#ifndef RAYLEIGH_FRICTION_UNIT_TESTS_COMMON_HPP
#define RAYLEIGH_FRICTION_UNIT_TESTS_COMMON_HPP

#include "rayleigh_friction_functions.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/physics/physics_test_data.hpp"

namespace scream {
namespace rayleigh_friction {
namespace unit_test {

/*
 * Unit test infrastructure for rayleigh_friction unit tests.
 *
 * rayleigh_friction entities can friend scream::rayleigh_friction::unit_test::UnitWrap
 * to give unit tests access to private members.
 *
 * All unit test impls should be within an inner struct of UnitWrap::UnitTest for
 * easy access to useful types.
 */

struct UnitWrap {

  template <typename D=DefaultDevice>
  struct UnitTest : public KokkosTypes<D> {

    using Device      = D;
    using MemberType  = typename KokkosTypes<Device>::MemberType;
    using TeamPolicy  = typename KokkosTypes<Device>::TeamPolicy;
    using RangePolicy = typename KokkosTypes<Device>::RangePolicy;
    using ExeSpace    = typename KokkosTypes<Device>::ExeSpace;

    template <typename S>
    using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
    template <typename S>
    using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
    template <typename S>
    using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

    using Functions = scream::rayleigh_friction::Functions<Real, Device>;
    using Scalar    = typename Functions::Scalar;

    struct Base : public UnitBase {

      Base() :
        UnitBase()
      {
        // no global rayleigh_friction data to initialize
      }

      ~Base() = default;
    };

    // Put struct decls here
    struct TestRayleighFrictionInit;
    struct TestRayleighFrictionTend;
  };

};

} // namespace unit_test
} // namespace rayleigh_friction
} // namespace scream

#endif
