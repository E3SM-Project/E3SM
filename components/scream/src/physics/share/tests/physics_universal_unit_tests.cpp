
#include "catch2/catch.hpp"

#include "physics/share/physics_functions.hpp"
#include "physics/share/physics_universal_impl.hpp"
#include "physics_unit_tests_common.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace scream {
namespace physics {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestUniversal
{

// Test Ideas:
// 1. Test exner:  
// 	Randomize pressures between 200 and 1100
// 	Caclulate Exner directly
// 	Calculate Exner again using universal function
// 	Compare
// 2. Test T_to_th_atm and th_atm_to_T
// 	Make dummy exner value of 2.0
// 	Randomize temperature, T_atm,  between 250 and 350
// 	Caclulate th_atm using dummy exner and T_atm
// 	th_atm should be 2x bigger than T_atm
// 	Calculate a new T_atm using T_atm and dummy exner
// 	the new T_atm should be half the size of the original T_atm
// 	Test error capability:
// 	  pass exner = 0, should throw an error
// 	  pass exner = -1, should throw an error
// 	  pass T_atm = 0, should throw an error
// 	  pass th_atm = 0, should throw an error
// 3. Test T_to_th vs. th_to_T inverse properties
// 	Randomize pressures between 200 and 1100
// 	Randomize temperature between 250 and 350
// 	make a 2d array of all T,P combinations
// 	Caclulate Exner from 2d P array
// 	Calculate th_atm from 2d T array and Exner array
// 	Calculate new T from th_atm and Exner
// 	new T and old T should match.
// 4. Test zi thickness calculation
// 	Create a set of zi with thickness growing by 1 each level, note, typical scream goes from top to bottom:
// 	zi(0) = N*N/2, for k=1...N, zi(k) = zi(k-1)-k
// 	Calculate dz based on zi, the layer thickness should monotonically increase by 1.
// 	Pass zi_top and zi_bot in reverse order, should throw an error, negative pressure thickness.

  static void run()
  {
    int nerr = 0;
    REQUIRE(nerr==0);
  } // run
}; // end of TestUniversal struct

} // namespace unit_test
} // namespace physics
} // namespace scream

namespace{

TEST_CASE("physics_universal_test", "[physics_universal_test]"){
  scream::physics::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUniversal::run();

 } // TEST_CASE

} // namespace
