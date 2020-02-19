#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iomanip>      // std::setprecision

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Unit-tests for p3 ice collection functions.
 */
template <typename D>
struct UnitWrap::UnitTest<D>::TestRainSelfCollection {

  static void run_rain_self_collection_bfb(){

    using KTH = KokkosTypes<HostDevice>;

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);
    
    RainSelfCollectionData dc[max_pack_size] = {
    //  rho, qr_incld, nr_incld, nrslf
    {}

    };

  }


}; //TestRainselfCollection



} // scream
} // p3
} // unit_test