#include "catch2/catch.hpp"
#include "physics/rrtmgp/rrtmgp_heating_rate.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
//#include "physics/rrtmgp/zenith.hpp"
#include "YAKL/YAKL.h"
#include "physics/share/physics_constants.hpp"
#include "physics/rrtmgp/share/shr_orb_mod_c.hpp"
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
    auto cloud_mass_ref = mixing_ratio.createHostCopy()(1,1) / cloud_fraction.createHostCopy()(1,1) * dp.createHostCopy()(1,1) / physconst::gravit;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass.createHostCopy()(1,1) == cloud_mass_ref);

    // Test with no cloud
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0;
        cloud_fraction(1,1) = 0.0;
    });
    cloud_mass_ref = 0.0;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass.createHostCopy()(1,1) == cloud_mass_ref);
 
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
    REQUIRE(cloud_mass.createHostCopy()(1,1) == cloud_mass_ref);
 
    // Test with cell half filled with cloud
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0001;
        cloud_fraction(1,1) = 0.5;
    });
    cloud_mass_ref = mixing_ratio.createHostCopy()(1,1) / cloud_fraction.createHostCopy()(1,1) * dp.createHostCopy()(1,1) / physconst::gravit;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass.createHostCopy()(1,1) == cloud_mass_ref);

    // Clean up
    dp.deallocate();
    mixing_ratio.deallocate();
    cloud_fraction.deallocate();
    cloud_mass.deallocate();
    yakl::finalize();
}

TEST_CASE("rrtmgp_test_limit_to_bounds") {
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }

    // Test limiter function
    auto arr = real2d("arr", 2, 2);
    auto arr_limited = real2d("arr_limited", 2, 2);

    // Setup dummy array
    parallel_for(1, YAKL_LAMBDA(int dummy) {
        arr(1,1) = 1.0;
        arr(1,2) = 2.0;
        arr(2,1) = 3.0;
        arr(2,2) = 4.0;
    });

    // Limit to bounds that contain the data; should be no change in values
    scream::rrtmgp::limit_to_bounds(arr, 0.0, 5.0, arr_limited);
    REQUIRE(arr.createHostCopy()(1,1) == arr_limited.createHostCopy()(1,1));
    REQUIRE(arr.createHostCopy()(1,2) == arr_limited.createHostCopy()(1,2));
    REQUIRE(arr.createHostCopy()(2,1) == arr_limited.createHostCopy()(2,1));
    REQUIRE(arr.createHostCopy()(2,2) == arr_limited.createHostCopy()(2,2));

    // Limit to bounds that do not completely contain the data; should be a change in values!
    scream::rrtmgp::limit_to_bounds(arr, 1.5, 3.5, arr_limited);
    REQUIRE(arr_limited.createHostCopy()(1,1) == 1.5);
    REQUIRE(arr_limited.createHostCopy()(1,2) == 2.0);
    REQUIRE(arr_limited.createHostCopy()(2,1) == 3.0);
    REQUIRE(arr_limited.createHostCopy()(2,2) == 3.5);
    arr.deallocate();
    arr_limited.deallocate();
    yakl::finalize();
}

TEST_CASE("rrtmgp_test_zenith") {

    // Create some dummy data
    int orbital_year = 1990;
    double calday = 1.0000000000000000;
    double eccen_ref = 1.6707719799280658E-002;
    double mvelpp_ref = 4.9344679089867318;
    double lambm0_ref = -3.2503635878519378E-002;
    double obliqr_ref = 0.40912382465788016;
    double delta_ref = -0.40302893695478670;
    double eccf_ref = 1.0342222039093694;
    double lat = -7.7397590528644963E-002;
    double lon = 2.2584340271163548;
    double coszrs_ref = 0.61243613606766745;

    // Test shr_orb_params()
    // Get orbital parameters based on calendar day
    double eccen;
    double obliq;  // obliquity in degrees
    double mvelp;  // moving vernal equinox long of perihelion; degrees?
    double obliqr;
    double lambm0;
    double mvelpp;
    bool flag_print = false;
    shr_orb_params_c(&orbital_year, &eccen, &obliq, &mvelp, 
                     &obliqr, &lambm0, &mvelpp); //, flag_print); // Note fortran code has optional arg
    REQUIRE(eccen == eccen_ref);
    REQUIRE(obliqr == obliqr_ref);
    REQUIRE(mvelpp == mvelpp_ref);
    REQUIRE(lambm0 == lambm0_ref);
    REQUIRE(mvelpp == mvelpp_ref);

    // Test shr_orb_decl()
    double delta;
    double eccf;
    shr_orb_decl_c(calday, eccen, mvelpp, lambm0,
                   obliqr, &delta, &eccf);
    REQUIRE(delta == delta_ref);
    REQUIRE(eccf  == eccf_ref );

    double dt_avg = 0.; //3600.0000000000000;
    double coszrs = shr_orb_cosz_c(calday, lat, lon, delta, dt_avg);
    REQUIRE(coszrs == coszrs_ref);

    // Another case, this time WITH dt_avg flag:
    //  DEBUG: calday, eccen, mvelpp, lambm0, obliqr, delta, eccf, lat, lon, dt_avg, coszrs =    1.0833333333333333        1.67077
    //  19799280658E-002   4.9344679089867318       -3.2503635878519378E-002  0.40912382465788016      -0.40292121709083456
    //  1.0342248931660425       -1.0724153591027763        4.5284876076962712        3600.0000000000000       0.14559973262047626
}
