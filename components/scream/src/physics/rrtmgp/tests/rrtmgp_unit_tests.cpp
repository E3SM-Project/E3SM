#include "catch2/catch.hpp"
#include "physics/rrtmgp/rrtmgp_heating_rate.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "YAKL/YAKL.h"
#include "physics/share/physics_constants.hpp"
TEST_CASE("rrtmgp_test_heating") {
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }

    // Test heating rate function by passing simple inputs
    auto dp = real2d("dp", 1, 1);
    auto flux_up = real2d("flux_up", 1, 2);
    auto flux_dn = real2d("flux_dn", 1, 2);
    auto heating = real2d("heating", 1, 1);
    // Simple no-heating test
    // NOTE: parallel_for because we need to do these in a kernel on the device
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        dp(1, 1) = 10;
        flux_up(1, 1) = 1.0;
        flux_up(1, 2) = 1.0;
        flux_dn(1, 1) = 1.0;
        flux_dn(1, 2) = 1.0;
    });
    scream::rrtmgp::compute_heating_rate(flux_up, flux_dn, dp, heating);
    REQUIRE(heating.createHostCopy()(1,1) == 0);

    // Simple net postive heating; net flux into layer should be 1.0
    // NOTE: parallel_for because we need to do these in a kernel on the device
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        flux_up(1, 1) = 1.0;
        flux_up(1, 2) = 1.0;
        flux_dn(1, 1) = 1.5;
        flux_dn(1, 2) = 0.5;
    });
    using physconst = scream::physics::Constants<double>;
    auto g = physconst::gravit; //9.81;
    auto cp_air = physconst::Cpair; //1005.0;
    auto pdel = dp.createHostCopy()(1,1);
    auto heating_ref = 1.0 * g / (cp_air * pdel);
    scream::rrtmgp::compute_heating_rate(flux_up, flux_dn, dp, heating);
    REQUIRE(heating.createHostCopy()(1,1) == heating_ref);

    // Simple net negative heating; net flux into layer should be -1.0
    // NOTE: parallel_for because we need to do these in a kernel on the device
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        flux_up(1,1) = 1.5;
        flux_up(1,2) = 0.5;
        flux_dn(1,1) = 1.0;
        flux_dn(1,2) = 1.0;
    });
    heating_ref = -1.0 * g / (cp_air * pdel);
    scream::rrtmgp::compute_heating_rate(flux_up, flux_dn, dp, heating);
    REQUIRE(heating.createHostCopy()(1,1) == heating_ref);

    // Clean up
    dp.deallocate();
    flux_up.deallocate();
    flux_dn.deallocate();
    heating.deallocate();
    yakl::finalize();
}

TEST_CASE("rrtmgp_test_mixing_ratio_to_cloud_mass") {
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }

    using physconst = scream::physics::Constants<double>;

    // Test mixing ratio to cloud mass function by passing simple inputs
    auto dp = real2d("dp", 1, 1);
    auto mixing_ratio = real2d("mixing_ratio", 1, 1);
    auto cloud_fraction = real2d("cloud_fration", 1, 1);
    auto cloud_mass = real2d("cloud_mass", 1, 1);

    // Test with cell completely filled with cloud
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0001;
        cloud_fraction(1,1) = 1.0;
    });
    auto cloud_mass_ref = mixing_ratio(1,1) * dp(1,1) / physconst::gravit;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass(1,1) == cloud_mass_ref);

    // Test with no cloud
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0;
        cloud_fraction(1,1) = 0.0;
    });
    cloud_mass_ref = 0.0;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass(1,1) == cloud_mass_ref);
 
     // Test with empty clouds (cloud fraction but with no associated mixing ratio)
     // This case could happen if we use a total cloud fraction, but compute layer
     // cloud mass separately for liquid and ice.
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0;
        cloud_fraction(1,1) = 0.1;
    });
    cloud_mass_ref = 0.0;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass(1,1) == cloud_mass_ref);
 
    // Test with cell half filled with cloud
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0001;
        cloud_fraction(1,1) = 0.5;
    });
    cloud_mass_ref = mixing_ratio(1,1) / cloud_fraction(1,1) * dp(1,1) / physconst::gravit;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass(1,1) == cloud_mass_ref);
}
