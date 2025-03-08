#include "catch2/catch.hpp"
#include "physics/rrtmgp/rrtmgp_utils.hpp"
#include "physics/rrtmgp/eamxx_rrtmgp_interface.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"
#include "mo_load_coefficients.h"

#ifdef RRTMGP_ENABLE_YAKL
#include "YAKL.h"
#endif

namespace {

#ifdef RRTMGP_ENABLE_KOKKOS
template <typename View>
auto chc(const View& view)
{
  return Kokkos::create_mirror_view_and_copy(HostDevice(), view);
}
#endif

// Names of input files we will need.
std::string coefficients_file_sw = SCREAM_DATA_DIR "/init/rrtmgp-data-sw-g112-210809.nc";
std::string coefficients_file_lw = SCREAM_DATA_DIR "/init/rrtmgp-data-lw-g128-210809.nc";
std::string cloud_optics_file_sw = SCREAM_DATA_DIR "/init/rrtmgp-cloud-optics-coeffs-sw.nc";
std::string cloud_optics_file_lw = SCREAM_DATA_DIR "/init/rrtmgp-cloud-optics-coeffs-lw.nc";

#ifdef RRTMGP_ENABLE_YAKL
TEST_CASE("rrtmgp_test_heating") {
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }

    // Test heating rate function by passing simple inputs
    auto dp = real2d("dp", 1, 1);
    auto flux_up = real2d("flux_up", 1, 2);
    auto flux_dn = real2d("flux_dn", 1, 2);
    auto heating = real2d("heating", 1, 1);
    // Simple no-heating test
    // NOTE: yakl::fortran::parallel_for because we need to do these in a kernel on the device
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        dp(1, 1) = 10;
        flux_up(1, 1) = 1.0;
        flux_up(1, 2) = 1.0;
        flux_dn(1, 1) = 1.0;
        flux_dn(1, 2) = 1.0;
    });
    scream::rrtmgp::compute_heating_rate(flux_up, flux_dn, dp, heating);
    REQUIRE(heating.createHostCopy()(1,1) == 0);

    // Simple net postive heating; net flux into layer should be 1.0
    // NOTE: yakl::fortran::parallel_for because we need to do these in a kernel on the device
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
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
    // NOTE: yakl::fortran::parallel_for because we need to do these in a kernel on the device
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
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
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0001;
        cloud_fraction(1,1) = 1.0;
    });
    auto cloud_mass_ref = mixing_ratio.createHostCopy()(1,1) / cloud_fraction.createHostCopy()(1,1) * dp.createHostCopy()(1,1) / physconst::gravit;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass.createHostCopy()(1,1) == cloud_mass_ref);

    // Test with no cloud
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
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
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        dp(1,1) = 10.0;
        mixing_ratio(1,1) = 0.0;
        cloud_fraction(1,1) = 0.1;
    });
    cloud_mass_ref = 0.0;
    scream::rrtmgp::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
    REQUIRE(cloud_mass.createHostCopy()(1,1) == cloud_mass_ref);

    // Test with cell half filled with cloud
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
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
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
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
    // bool flag_print = false;
    shr_orb_params_c2f(&orbital_year, &eccen, &obliq, &mvelp,
                     &obliqr, &lambm0, &mvelpp); //, flag_print); // Note fortran code has optional arg
    REQUIRE(eccen == eccen_ref);
    REQUIRE(obliqr == obliqr_ref);
    REQUIRE(mvelpp == mvelpp_ref);
    REQUIRE(lambm0 == lambm0_ref);
    REQUIRE(mvelpp == mvelpp_ref);

    // Test shr_orb_decl()
    double delta;
    double eccf;
    shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0,
                   obliqr, &delta, &eccf);
    REQUIRE(delta == delta_ref);
    REQUIRE(eccf  == eccf_ref );

    double dt_avg = 0.; //3600.0000000000000;
    double coszrs = shr_orb_cosz_c2f(calday, lat, lon, delta, dt_avg);
    REQUIRE(std::abs(coszrs-coszrs_ref)<1e-14);

    // Another case, this time WITH dt_avg flag:
    calday = 1.0833333333333333;
    eccen = 1.6707719799280658E-002;
    mvelpp = 4.9344679089867318;
    lambm0 = -3.2503635878519378E-002;
    obliqr = 0.40912382465788016;
    delta = -0.40292121709083456;
    eccf = 1.0342248931660425;
    lat = -1.0724153591027763;
    lon = 4.5284876076962712;
    dt_avg = 3600.0000000000000;
    coszrs_ref = 0.14559973262047626;
    coszrs = shr_orb_cosz_c2f(calday, lat, lon, delta, dt_avg);
    REQUIRE(std::abs(coszrs-coszrs_ref)<1e-14);

}

TEST_CASE("rrtmgp_test_compute_broadband_surface_flux") {
    using namespace ekat::logger;
    using logger_t = Logger<LogNoFile,LogRootRank>;

    ekat::Comm comm(MPI_COMM_WORLD);
    auto logger = std::make_shared<logger_t>("",LogLevel::info,comm);

    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }

    // Create arrays
    const int ncol = 1;
    const int nlay = 1;
    const int nbnd = 14;
    const int kbot = nlay + 1;
    auto sfc_flux_dir_nir = real1d("sfc_flux_dir_nir", ncol);
    auto sfc_flux_dir_vis = real1d("sfc_flux_dir_vis", ncol);
    auto sfc_flux_dif_nir = real1d("sfc_flux_dif_nir", ncol);
    auto sfc_flux_dif_vis = real1d("sfc_flux_dif_vis", ncol);

    // Need to initialize RRTMGP with dummy gases
    logger->info("Init gases...\n");
    GasConcs gas_concs;
    string1dv gas_names = {"h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"};
    gas_concs.init(gas_names,ncol,nlay);
    logger->info("Init RRTMGP...\n");
    scream::rrtmgp::rrtmgp_initialize(gas_concs, coefficients_file_sw, coefficients_file_lw, cloud_optics_file_sw, cloud_optics_file_lw, logger);

    // Create simple test cases; We expect, given the input data, that band 10
    // will straddle the NIR and VIS, bands 1-9 will be purely NIR, and bands 11-14
    // will be purely VIS. The implementation in EAMF90 was hard-coded with this
    // band information, but our implementation of compute_broadband_surface_fluxes
    // actually checks the wavenumber limits. These tests will mostly check to make
    // sure our implementation of that is doing what we think it is.

    // ---------------------------------
    // Test case: flux only in straddled band
    auto sw_bnd_flux_dir = real3d("sw_bnd_flux_dir", ncol, nlay+1, nbnd);
    auto sw_bnd_flux_dif = real3d("sw_bnd_flux_dif", ncol, nlay+1, nbnd);
    logger->info("Populate band-resolved 3d fluxes for test case with only transition band flux...\n");
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<3>(nbnd,nlay+1,ncol), YAKL_LAMBDA(int ibnd, int ilay, int icol) {
        if (ibnd < 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
        } else if (ibnd == 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 1;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 1;
        } else {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
        }
    });
    // Compute surface fluxes
    logger->info("Compute broadband surface fluxes...\n");
    scream::rrtmgp::compute_broadband_surface_fluxes(
        ncol, kbot, nbnd,
        sw_bnd_flux_dir, sw_bnd_flux_dif,
        sfc_flux_dir_vis, sfc_flux_dir_nir,
        sfc_flux_dif_vis, sfc_flux_dif_nir
    );
    // Check computed surface fluxes
    logger->info("Check computed fluxes...\n");
    const double tol = 1e-10;  // tolerance on floating point inequality for assertions
    REQUIRE(std::abs(sfc_flux_dir_nir.createHostCopy()(1) - 0.5) < tol);
    REQUIRE(std::abs(sfc_flux_dir_vis.createHostCopy()(1) - 0.5) < tol);
    REQUIRE(std::abs(sfc_flux_dif_nir.createHostCopy()(1) - 0.5) < tol);
    REQUIRE(std::abs(sfc_flux_dif_vis.createHostCopy()(1) - 0.5) < tol);
    // ---------------------------------

    // ---------------------------------
    // Test case, only flux in NIR bands
    logger->info("Populate band-resolved 3d fluxes for test case with only NIR flux...\n");
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<3>(nbnd,nlay+1,ncol), YAKL_LAMBDA(int ibnd, int ilay, int icol) {
        if (ibnd < 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 1;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 1;
        } else if (ibnd == 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
        } else {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
        }
    });
    // Compute surface fluxes
    logger->info("Compute broadband surface fluxes...\n");
    scream::rrtmgp::compute_broadband_surface_fluxes(
        ncol, kbot, nbnd,
        sw_bnd_flux_dir, sw_bnd_flux_dif,
        sfc_flux_dir_vis, sfc_flux_dir_nir,
        sfc_flux_dif_vis, sfc_flux_dif_nir
    );
    // Check computed surface fluxes
    logger->info("Check computed fluxes...\n");
    REQUIRE(std::abs(sfc_flux_dir_nir.createHostCopy()(1) - 9.0) < tol);
    REQUIRE(std::abs(sfc_flux_dir_vis.createHostCopy()(1) - 0.0) < tol);
    REQUIRE(std::abs(sfc_flux_dif_nir.createHostCopy()(1) - 9.0) < tol);
    REQUIRE(std::abs(sfc_flux_dif_vis.createHostCopy()(1) - 0.0) < tol);
    // ---------------------------------

    // ---------------------------------
    // Test case, only flux in VIS bands
    logger->info("Populate band-resolved 3d fluxes for test case with only VIS/UV flux...\n");
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<3>(nbnd,nlay+1,ncol), YAKL_LAMBDA(int ibnd, int ilay, int icol) {
        if (ibnd < 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
        } else if (ibnd == 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
        } else {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 1;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 1;
        }
    });
    // Compute surface fluxes
    logger->info("Compute broadband surface fluxes...\n");
    scream::rrtmgp::compute_broadband_surface_fluxes(
        ncol, kbot, nbnd,
        sw_bnd_flux_dir, sw_bnd_flux_dif,
        sfc_flux_dir_vis, sfc_flux_dir_nir,
        sfc_flux_dif_vis, sfc_flux_dif_nir
    );
    // Check computed surface fluxes
    logger->info("Check computed fluxes...\n");
    REQUIRE(std::abs(sfc_flux_dir_nir.createHostCopy()(1) - 0.0) < tol);
    REQUIRE(std::abs(sfc_flux_dir_vis.createHostCopy()(1) - 4.0) < tol);
    REQUIRE(std::abs(sfc_flux_dif_nir.createHostCopy()(1) - 0.0) < tol);
    REQUIRE(std::abs(sfc_flux_dif_vis.createHostCopy()(1) - 4.0) < tol);
    // ---------------------------------

    // ---------------------------------
    // Test case, only flux in all bands
    logger->info("Populate band-resolved 3d fluxes for test with non-zero flux in all bands...\n");
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<3>(nbnd,nlay+1,ncol), YAKL_LAMBDA(int ibnd, int ilay, int icol) {
        if (ibnd < 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 1.0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 2.0;
        } else if (ibnd == 10) {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 3.0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 4.0;
        } else {
            sw_bnd_flux_dir(icol,ilay,ibnd) = 5.0;
            sw_bnd_flux_dif(icol,ilay,ibnd) = 6.0;
        }
    });
    // Compute surface fluxes
    logger->info("Compute broadband surface fluxes...\n");
    scream::rrtmgp::compute_broadband_surface_fluxes(
        ncol, kbot, nbnd,
        sw_bnd_flux_dir, sw_bnd_flux_dif,
        sfc_flux_dir_vis, sfc_flux_dir_nir,
        sfc_flux_dif_vis, sfc_flux_dif_nir
    );
    // Check computed surface fluxes
    logger->info("Check computed fluxes...\n");
    REQUIRE(std::abs(sfc_flux_dir_nir.createHostCopy()(1) - 10.5) < tol);
    REQUIRE(std::abs(sfc_flux_dir_vis.createHostCopy()(1) - 21.5) < tol);
    REQUIRE(std::abs(sfc_flux_dif_nir.createHostCopy()(1) - 20.0) < tol);
    REQUIRE(std::abs(sfc_flux_dif_vis.createHostCopy()(1) - 26.0) < tol);
    // ---------------------------------

    // Finalize YAKL
    logger->info("Free memory...\n");
    scream::rrtmgp::rrtmgp_finalize();
    gas_concs.reset();
    sw_bnd_flux_dir.deallocate();
    sw_bnd_flux_dif.deallocate();
    sfc_flux_dir_nir.deallocate();
    sfc_flux_dir_vis.deallocate();
    sfc_flux_dif_nir.deallocate();
    sfc_flux_dif_vis.deallocate();
    if (yakl::isInitialized()) { yakl::finalize(); }
}

TEST_CASE("rrtmgp_test_radiation_do") {
    // If we specify rad every step, radiation_do should always be true
    REQUIRE(scream::rrtmgp::radiation_do(1, 0) == true);
    REQUIRE(scream::rrtmgp::radiation_do(1, 1) == true);
    REQUIRE(scream::rrtmgp::radiation_do(1, 2) == true);

    // Test cases where we want rad called every other step
    REQUIRE(scream::rrtmgp::radiation_do(2, 0) == true);
    REQUIRE(scream::rrtmgp::radiation_do(2, 1) == false);
    REQUIRE(scream::rrtmgp::radiation_do(2, 2) == true);
    REQUIRE(scream::rrtmgp::radiation_do(2, 3) == false);

    // Test cases where we want rad every third step
    REQUIRE(scream::rrtmgp::radiation_do(3, 0) == true);
    REQUIRE(scream::rrtmgp::radiation_do(3, 1) == false);
    REQUIRE(scream::rrtmgp::radiation_do(3, 2) == false);
    REQUIRE(scream::rrtmgp::radiation_do(3, 3) == true);
    REQUIRE(scream::rrtmgp::radiation_do(3, 4) == false);
    REQUIRE(scream::rrtmgp::radiation_do(3, 5) == false);
    REQUIRE(scream::rrtmgp::radiation_do(3, 6) == true);
}

TEST_CASE("rrtmgp_test_check_range") {
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }
    // Create some dummy data and test with both values inside valid range and outside
    auto dummy = real2d("dummy", 2, 1);
    // All values within range
    memset(dummy, 0.1);
    REQUIRE(scream::rrtmgp::check_range(dummy, 0.0, 1.0, "dummy") == true);
    // At least one value below lower bound
    yakl::fortran::parallel_for(1, YAKL_LAMBDA (int i) {dummy(i, 1) = -0.1;});
    REQUIRE(scream::rrtmgp::check_range(dummy, 0.0, 1.0, "dummy") == false);
    // At least one value above upper bound
    yakl::fortran::parallel_for(1, YAKL_LAMBDA (int i) {dummy(i, 1) = 1.1;});
    REQUIRE(scream::rrtmgp::check_range(dummy, 0.0, 1.0, "dummy") == false);
    dummy.deallocate();
    if (yakl::isInitialized()) { yakl::finalize(); }
}

TEST_CASE("rrtmgp_test_subcol_gen") {
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }
    // Create dummy data
    const int ncol = 1;
    const int nlay = 4;
    const int ngpt = 10;
    auto cldfrac = real2d("cldfrac", ncol, nlay);
    // Set cldfrac values
    memset(cldfrac, 0.0);
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac(1,1) = 1;
        cldfrac(1,2) = 0.5;
        cldfrac(1,3) = 0;
        cldfrac(1,4) = 1;
    });
    auto cldmask = int3d("cldmask", ncol, nlay, ngpt);
    auto cldfrac_from_mask = real2d("cldfrac_from_mask", ncol, nlay);
    // Run subcol gen, make sure we get what we expect; do this for some different seed values
    for (unsigned seed = 1; seed <= 10; seed++) {
        auto seeds = int1d("seeds", ncol);
        memset(seeds, seed);
        cldmask = scream::rrtmgp::get_subcolumn_mask(ncol, nlay, ngpt, cldfrac, 1, seeds);
        // Check answers by computing new cldfrac from mask
        memset(cldfrac_from_mask, 0.0);
        yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
            for (int igpt = 1; igpt <= ngpt; ++igpt) {
                real cldmask_real = cldmask(icol,ilay,igpt);
                cldfrac_from_mask(icol,ilay) += cldmask_real;
            }
        });
        yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
            cldfrac_from_mask(icol,ilay) = cldfrac_from_mask(icol,ilay) / ngpt;
        });
        // For cldfrac 1 we should get 1, for cldfrac 0 we should get 0, but in between we cannot be sure
        // deterministically, since the computed cloud mask depends on pseudo-random numbers
        REQUIRE(cldfrac_from_mask.createHostCopy()(1,1) == 1);
        REQUIRE(cldfrac_from_mask.createHostCopy()(1,2) <= 1);
        REQUIRE(cldfrac_from_mask.createHostCopy()(1,3) == 0);
        REQUIRE(cldfrac_from_mask.createHostCopy()(1,4) == 1);
    }

    // For maximum-random overlap, vertically-contiguous layers maximimally overlap,
    // thus if we have non-zero cloud fraction in two adjacent layers, then every subcolumn
    // that has cloud in the layer above must also have cloud in the layer below; test
    // this property by creating two layers with non-zero cloud fraction, creating subcolums,
    // and verifying that every subcolumn with cloud in layer 1 has cloud in layer 2
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac(1,1) = 0.5;
        cldfrac(1,2) = 0.5;
        cldfrac(1,3) = 0;
        cldfrac(1,4) = 0;
    });
    for (unsigned seed = 1; seed <= 10; seed++) {
        auto seeds = int1d("seeds", ncol);
        memset(seeds, seed);
        cldmask = scream::rrtmgp::get_subcolumn_mask(ncol, nlay, ngpt, cldfrac, 1, seeds);
        auto cldmask_h = cldmask.createHostCopy();
        for (int igpt = 1; igpt <= ngpt; igpt++) {
            if (cldmask_h(1,1,igpt) == 1) {
                REQUIRE(cldmask_h(1,2,igpt) == 1);
            }
        }
    }
    // Clean up after test
    cldfrac.deallocate();
    cldmask.deallocate();
    cldfrac_from_mask.deallocate();
    yakl::finalize();
}


TEST_CASE("rrtmgp_cloud_area") {
    // Initialize YAKL
    if (!yakl::isInitialized()) { yakl::init(); }
    // Create dummy data
    const int ncol = 1;
    const int nlay = 2;
    const int ngpt = 3;
    auto cldtau = real3d("cldtau", ncol, nlay, ngpt);
    auto cldtot = real1d("cldtot", ncol);
    auto pmid = real2d("pmid", ncol, nlay);

    // Set up pressure levels for test problem
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        pmid(1,1) = 100;
        pmid(1,2) = 200;
    });

    // Case:
    //
    // 0 0 0
    // 0 0 0
    //
    // should give cldtot = 0.0
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        cldtau(1,1,1) = 0;
        cldtau(1,1,2) = 0;
        cldtau(1,1,3) = 0;
        cldtau(1,2,1) = 0;
        cldtau(1,2,2) = 0;
        cldtau(1,2,3) = 0;
    });
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 0.0);

    // Case:
    //
    // 1 1 1
    // 1 1 1
    //
    // should give cldtot = 1.0
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        cldtau(1,1,1) = 1;
        cldtau(1,1,2) = 1;
        cldtau(1,1,3) = 1;
        cldtau(1,2,1) = 1;
        cldtau(1,2,2) = 1;
        cldtau(1,2,3) = 1;
    });
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 1.0);

    // Case:
    //
    // 1 1 0  100
    // 0 0 1  200
    //
    // should give cldtot = 1.0
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        cldtau(1,1,1) = 0.1;
        cldtau(1,1,2) = 1.5;
        cldtau(1,1,3) = 0;
        cldtau(1,2,1) = 0;
        cldtau(1,2,2) = 0;
        cldtau(1,2,3) = 1.0;
    });
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 1.0);
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 0, 150, pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 2.0 / 3.0);
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 110, 250, pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 1.0 / 3.0);

    // Case:
    //
    // 1 0 0
    // 1 0 1
    //
    // should give cldtot = 2/3
    yakl::fortran::parallel_for(1, YAKL_LAMBDA(int /* dummy */) {
        cldtau(1,1,1) = 1;
        cldtau(1,1,2) = 0;
        cldtau(1,1,3) = 0;
        cldtau(1,2,1) = 1;
        cldtau(1,2,2) = 0;
        cldtau(1,2,3) = 1;
    });
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 2.0 / 3.0);
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 0, 100, pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 0.0);
    scream::rrtmgp::compute_cloud_area(ncol, nlay, ngpt, 100, 300, pmid, cldtau, cldtot);
    REQUIRE(cldtot.createHostCopy()(1) == 2.0 / 3.0);
    pmid.deallocate();
    cldtau.deallocate();
    cldtot.deallocate();
    yakl::finalize();
}

TEST_CASE("rrtmgp_aerocom_cloudtop") {
  // Initialize YAKL
  if(!yakl::isInitialized()) {
    yakl::init();
  }
  // Create dummy data
  const int ncol = 1;
  const int nlay = 9;
  // Set up input fields
  auto tmid        = real2d("tmid", ncol, nlay);
  auto pmid        = real2d("pmid", ncol, nlay);
  auto p_del       = real2d("p_del", ncol, nlay);
  auto z_del       = real2d("z_del", ncol, nlay);
  auto qc          = real2d("qc", ncol, nlay);
  auto qi          = real2d("qi", ncol, nlay);
  auto rel         = real2d("rel", ncol, nlay);
  auto rei         = real2d("rei", ncol, nlay);
  auto cldfrac_tot = real2d("cldfrac_tot", ncol, nlay);
  auto nc          = real2d("nc", ncol, nlay);
  // Set up output fields
  auto tmid_at_cldtop          = real1d("tmid_at_cldtop", ncol);
  auto pmid_at_cldtop          = real1d("pmid_at_cldtop", ncol);
  auto cldfrac_ice_at_cldtop   = real1d("cldfrac_ice_at_cldtop", ncol);
  auto cldfrac_liq_at_cldtop   = real1d("cldfrac_liq_at_cldtop", ncol);
  auto cldfrac_tot_at_cldtop   = real1d("cldfrac_tot_at_cldtop", ncol);
  auto cdnc_at_cldtop          = real1d("cdnc_at_cldtop", ncol);
  auto eff_radius_qc_at_cldtop = real1d("eff_radius_qc_at_cldtop", ncol);
  auto eff_radius_qi_at_cldtop = real1d("eff_radius_qi_at_cldtop", ncol);

  // Case 1: if no clouds, everything goes to zero
  memset(tmid, 300.0);
  memset(pmid, 100.0);
  memset(p_del, 10.0);
  memset(z_del, 100.0);
  memset(qc, 1.0);
  memset(qi, 1.0);
  memset(cldfrac_tot, 0.0);
  memset(nc, 5.0);
  memset(rel, 10.0);
  memset(rei, 10.0);
  // Call the function
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  // Check the results
  REQUIRE(tmid_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(pmid_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(cldfrac_liq_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(cldfrac_ice_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(cdnc_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(eff_radius_qc_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(eff_radius_qi_at_cldtop.createHostCopy()(1) == 0.0);

  // Case 2: if all clouds, everything goes to 1 * its value
  memset(cldfrac_tot, 1.0);
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(tmid_at_cldtop.createHostCopy()(1) == 300.0);
  REQUIRE(pmid_at_cldtop.createHostCopy()(1) == 100.0);
  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) == 1.0);
  REQUIRE(cldfrac_liq_at_cldtop.createHostCopy()(1) == 0.5);
  REQUIRE(cldfrac_ice_at_cldtop.createHostCopy()(1) == 0.5);
  REQUIRE(cdnc_at_cldtop.createHostCopy()(1) > 0.0);
  REQUIRE(eff_radius_qc_at_cldtop.createHostCopy()(1) > 0.0);
  REQUIRE(eff_radius_qi_at_cldtop.createHostCopy()(1) > 0.0);

  // Case 3: test max overlap (if contiguous cloudy layers, then max)
  memset(cldfrac_tot, 0.0);
  yakl::fortran::parallel_for(
      1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac_tot(1, 2) = 0.5;
        cldfrac_tot(1, 3) = 0.7;
        cldfrac_tot(1, 4) = 0.3;
        cldfrac_tot(1, 5) = 0.2;
      });
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) == .7);

  // Case 3xtra: test max overlap
  // This case produces >0.7 due to slight enhancement in the presence of a
  // local minimum (0.1 is the local minimum between 0.2 and 0.4)
  yakl::fortran::parallel_for(
      1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac_tot(1, 5) = 0.1;
        cldfrac_tot(1, 6) = 0.4;
        cldfrac_tot(1, 7) = 0.2;
      });
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) > .7);

  // Case 4: test random overlap (if non-contiguous cloudy layers, then
  // random)
  yakl::fortran::parallel_for(
      1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac_tot(1, 5) = 0.0;
        cldfrac_tot(1, 6) = 0.1;
      });
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) >
          .7);  // larger than the max

  // Case 5a: test independence of ice and liquid fractions
  yakl::fortran::parallel_for(
      1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac_tot(1, 2) = 1.0;
        cldfrac_tot(1, 7) = 1.0;
        cldfrac_tot(1, 8) = 0.2;
      });
  memset(qc, 1.0);
  memset(qi, 0.0);
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) == 1.0);
  REQUIRE(cldfrac_liq_at_cldtop.createHostCopy()(1) == 1.0);
  REQUIRE(cldfrac_ice_at_cldtop.createHostCopy()(1) == 0.0);

  // Case 5b: test independence of ice and liquid fractions
  yakl::fortran::parallel_for(
      1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac_tot(1, 2) = 1.0;
        cldfrac_tot(1, 7) = 1.0;
        cldfrac_tot(1, 8) = 0.2;
      });
  memset(qc, 0.0);
  memset(qi, 1.0);
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) == 1.0);
  REQUIRE(cldfrac_liq_at_cldtop.createHostCopy()(1) == 0.0);
  REQUIRE(cldfrac_ice_at_cldtop.createHostCopy()(1) == 1.0);

  // Case 6: test independence of ice and liquid fractions
  // There is NOT complete independence...
  // Essentially, higher ice clouds mask lower liquid clouds
  // This can be problematic if the ice clouds are thin...
  // We will revisit and validate this assumption later
  memset(cldfrac_tot, 0.0);
  memset(qc, 0.0);
  memset(qi, 0.0);
  yakl::fortran::parallel_for(
      1, YAKL_LAMBDA(int /* dummy */) {
        cldfrac_tot(1, 2) = 0.5;  // ice
        cldfrac_tot(1, 3) = 0.7;  // ice ------> max
        cldfrac_tot(1, 4) = 0.3;  // ice
        // note cldfrac_tot(1, 5) is 0
        cldfrac_tot(1, 6) = 0.2;  // liq
        cldfrac_tot(1, 7) = 0.5;  // liq ------> not max
        cldfrac_tot(1, 8) = 0.1;  // liq
        // note cldfrac_tot(1, 9) is 0
        qi(1, 2) = 100;
        qi(1, 3) = 200;
        qi(1, 4) = 50;
        // note qc(1, 5) is 0
        // note qi(1, 5) is 0
        qc(1, 6) = 20;
        qc(1, 7) = 50;
        qc(1, 8) = 10;
      });
  scream::rrtmgp::compute_aerocom_cloudtop(
      ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
      tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
      cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
      eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);
  REQUIRE(cldfrac_tot_at_cldtop.createHostCopy()(1) > 0.70);  // unaffected
  REQUIRE(cldfrac_liq_at_cldtop.createHostCopy()(1) < 0.50);  // not max
  REQUIRE(cldfrac_ice_at_cldtop.createHostCopy()(1) == 0.7);  // max

  // cleanup
  tmid.deallocate();
  pmid.deallocate();
  p_del.deallocate();
  z_del.deallocate();
  qc.deallocate();
  qi.deallocate();
  rel.deallocate();
  rei.deallocate();
  cldfrac_tot.deallocate();
  nc.deallocate();

  tmid_at_cldtop.deallocate();
  pmid_at_cldtop.deallocate();
  cldfrac_ice_at_cldtop.deallocate();
  cldfrac_liq_at_cldtop.deallocate();
  cldfrac_tot_at_cldtop.deallocate();
  cdnc_at_cldtop.deallocate();
  eff_radius_qc_at_cldtop.deallocate();
  eff_radius_qi_at_cldtop.deallocate();

  yakl::finalize();
}
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
using interface_t = scream::rrtmgp::rrtmgp_interface<>;
using pool_t = interface_t::pool_t;
using real1dk = interface_t::view_t<scream::Real*>;
using real2dk = interface_t::view_t<scream::Real**>;
using real3dk = interface_t::view_t<scream::Real***>;
using int1dk = interface_t::view_t<int*>;
using int2dk = interface_t::view_t<int**>;
using int3dk = interface_t::view_t<int***>;
using MDRP = interface_t::MDRP;

TEST_CASE("rrtmgp_test_heating_k") {
  // Initialize Kokkos
  scream::init_kls();
  pool_t::init(10000);

  // Test heating rate function by passing simple inputs
  auto dp = real2dk("dp", 1, 1);
  auto flux_up = real2dk("flux_up", 1, 2);
  auto flux_dn = real2dk("flux_dn", 1, 2);
  auto heating = real2dk("heating", 1, 1);
  // Simple no-heating test
  // NOTE: Kokkos::parallel_for because we need to do these in a kernel on the device
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    dp(0, 0) = 10;
    flux_up(0, 0) = 1.0;
    flux_up(0, 1) = 1.0;
    flux_dn(0, 0) = 1.0;
    flux_dn(0, 1) = 1.0;
  });
  scream::rrtmgp::compute_heating_rate(flux_up, flux_dn, dp, heating);
  REQUIRE(chc(heating)(0,0) == 0);

  // Simple net postive heating; net flux into layer should be 1.0
  // NOTE: Kokkos::parallel_for because we need to do these in a kernel on the device
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    flux_up(0, 0) = 1.0;
    flux_up(0, 1) = 1.0;
    flux_dn(0, 0) = 1.5;
    flux_dn(0, 1) = 0.5;
  });
  using physconst = scream::physics::Constants<double>;
  auto g = physconst::gravit; //9.81;
  auto cp_air = physconst::Cpair; //1005.0;
  auto pdel = chc(dp)(0,0);
  auto heating_ref = 1.0 * g / (cp_air * pdel);
  scream::rrtmgp::compute_heating_rate(flux_up, flux_dn, dp, heating);
  REQUIRE(chc(heating)(0,0) == heating_ref);

  // Simple net negative heating; net flux into layer should be -1.0
  // NOTE: Kokkos::parallel_for because we need to do these in a kernel on the device
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    flux_up(0,0) = 1.5;
    flux_up(0,1) = 0.5;
    flux_dn(0,0) = 1.0;
    flux_dn(0,1) = 1.0;
  });
  heating_ref = -1.0 * g / (cp_air * pdel);
  scream::rrtmgp::compute_heating_rate(flux_up, flux_dn, dp, heating);
  REQUIRE(chc(heating)(0,0) == heating_ref);

  // Clean up
  pool_t::finalize();
  scream::finalize_kls();
}

TEST_CASE("rrtmgp_test_mixing_ratio_to_cloud_mass_k") {
  // Initialize YAKL
  scream::init_kls();
  pool_t::init(10000);

  using physconst = scream::physics::Constants<double>;

  // Test mixing ratio to cloud mass function by passing simple inputs
  auto dp = real2dk("dp", 1, 1);
  auto mixing_ratio = real2dk("mixing_ratio", 1, 1);
  auto cloud_fraction = real2dk("cloud_fration", 1, 1);
  auto cloud_mass = real2dk("cloud_mass", 1, 1);

  // Test with cell completely filled with cloud
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    dp(0,0) = 10.0;
    mixing_ratio(0,0) = 0.0001;
    cloud_fraction(0,0) = 1.0;
  });
  auto cloud_mass_ref = chc(mixing_ratio)(0,0) / chc(cloud_fraction)(0,0) * chc(dp)(0,0) / physconst::gravit;
  interface_t::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
  REQUIRE(chc(cloud_mass)(0,0) == cloud_mass_ref);

  // Test with no cloud
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    dp(0,0) = 10.0;
    mixing_ratio(0,0) = 0.0;
    cloud_fraction(0,0) = 0.0;
  });
  cloud_mass_ref = 0.0;
  interface_t::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
  REQUIRE(chc(cloud_mass)(0,0) == cloud_mass_ref);

  // Test with empty clouds (cloud fraction but with no associated mixing ratio)
  // This case could happen if we use a total cloud fraction, but compute layer
  // cloud mass separately for liquid and ice.
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    dp(0,0) = 10.0;
    mixing_ratio(0,0) = 0.0;
    cloud_fraction(0,0) = 0.1;
  });
  cloud_mass_ref = 0.0;
  interface_t::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
  REQUIRE(chc(cloud_mass)(0,0) == cloud_mass_ref);

  // Test with cell half filled with cloud
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    dp(0,0) = 10.0;
    mixing_ratio(0,0) = 0.0001;
    cloud_fraction(0,0) = 0.5;
  });
  cloud_mass_ref = chc(mixing_ratio)(0,0) / chc(cloud_fraction)(0,0) * chc(dp)(0,0) / physconst::gravit;
  interface_t::mixing_ratio_to_cloud_mass(mixing_ratio, cloud_fraction, dp, cloud_mass);
  REQUIRE(chc(cloud_mass)(0,0) == cloud_mass_ref);

  // Clean up
  pool_t::finalize();
  scream::finalize_kls();
}

TEST_CASE("rrtmgp_test_limit_to_bounds_k") {
  // Initialize YAKL
  scream::init_kls();
  pool_t::init(10000);

  // Test limiter function
  auto arr = real2dk("arr", 2, 2);
  auto arr_limited = real2dk("arr_limited", 2, 2);

  // Setup dummy array
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    arr(0,0) = 1.0;
    arr(0,1) = 2.0;
    arr(1,0) = 3.0;
    arr(1,1) = 4.0;
  });

  // Limit to bounds that contain the data; should be no change in values
  interface_t::limit_to_bounds_k(arr, 0.0, 5.0, arr_limited);
  REQUIRE(chc(arr)(0,0) == chc(arr_limited)(0,0));
  REQUIRE(chc(arr)(0,1) == chc(arr_limited)(0,1));
  REQUIRE(chc(arr)(1,0) == chc(arr_limited)(1,0));
  REQUIRE(chc(arr)(1,1) == chc(arr_limited)(1,1));

  // Limit to bounds that do not completely contain the data; should be a change in values!
  interface_t::limit_to_bounds_k(arr, 1.5, 3.5, arr_limited);
  REQUIRE(chc(arr_limited)(0,0) == 1.5);
  REQUIRE(chc(arr_limited)(0,1) == 2.0);
  REQUIRE(chc(arr_limited)(1,0) == 3.0);
  REQUIRE(chc(arr_limited)(1,1) == 3.5);

  pool_t::finalize();
  scream::finalize_kls();
}

TEST_CASE("rrtmgp_test_zenith_k") {

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
  // bool flag_print = false;
  shr_orb_params_c2f(&orbital_year, &eccen, &obliq, &mvelp,
                     &obliqr, &lambm0, &mvelpp); //, flag_print); // Note fortran code has optional arg
  REQUIRE(eccen == eccen_ref);
  REQUIRE(obliqr == obliqr_ref);
  REQUIRE(mvelpp == mvelpp_ref);
  REQUIRE(lambm0 == lambm0_ref);
  REQUIRE(mvelpp == mvelpp_ref);

  // Test shr_orb_decl()
  double delta;
  double eccf;
  shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0,
                   obliqr, &delta, &eccf);
  REQUIRE(delta == delta_ref);
  REQUIRE(eccf  == eccf_ref );

  double dt_avg = 0.; //3600.0000000000000;
  double coszrs = shr_orb_cosz_c2f(calday, lat, lon, delta, dt_avg);
  REQUIRE(std::abs(coszrs-coszrs_ref)<1e-14);

  // Another case, this time WITH dt_avg flag:
  calday = 1.0833333333333333;
  eccen = 1.6707719799280658E-002;
  mvelpp = 4.9344679089867318;
  lambm0 = -3.2503635878519378E-002;
  obliqr = 0.40912382465788016;
  delta = -0.40292121709083456;
  eccf = 1.0342248931660425;
  lat = -1.0724153591027763;
  lon = 4.5284876076962712;
  dt_avg = 3600.0000000000000;
  coszrs_ref = 0.14559973262047626;
  coszrs = shr_orb_cosz_c2f(calday, lat, lon, delta, dt_avg);
  REQUIRE(std::abs(coszrs-coszrs_ref)<1e-14);
}

TEST_CASE("rrtmgp_test_compute_broadband_surface_flux_k") {
  using namespace ekat::logger;
  using logger_t = Logger<LogNoFile,LogRootRank>;

  ekat::Comm comm(MPI_COMM_WORLD);
  auto logger = std::make_shared<logger_t>("",LogLevel::info,comm);

  // Initialize YAKL
  scream::init_kls();
  pool_t::init(10000);

  // Create arrays
  const int ncol = 1;
  const int nlay = 1;
  const int nbnd = 14;
  const int kbot = nlay;
  auto sfc_flux_dir_nir = real1dk("sfc_flux_dir_nir", ncol);
  auto sfc_flux_dir_vis = real1dk("sfc_flux_dir_vis", ncol);
  auto sfc_flux_dif_nir = real1dk("sfc_flux_dif_nir", ncol);
  auto sfc_flux_dif_vis = real1dk("sfc_flux_dif_vis", ncol);

  // Need to initialize RRTMGP with dummy gases
  logger->info("Init gases...\n");
  GasConcsK<scream::Real, Kokkos::LayoutRight, DefaultDevice> gas_concs;
  string1dv gas_names = {"h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"};
  gas_concs.init(gas_names,ncol,nlay);
  logger->info("Init RRTMGP...\n");
  interface_t::rrtmgp_initialize(gas_concs, coefficients_file_sw, coefficients_file_lw, cloud_optics_file_sw, cloud_optics_file_lw, logger);

  // Create simple test cases; We expect, given the input data, that band 10
  // will straddle the NIR and VIS, bands 1-9 will be purely NIR, and bands 11-14
  // will be purely VIS. The implementation in EAMF90 was hard-coded with this
  // band information, but our implementation of compute_broadband_surface_fluxes
  // actually checks the wavenumber limits. These tests will mostly check to make
  // sure our implementation of that is doing what we think it is.

  // ---------------------------------
  // Test case: flux only in straddled band
  auto sw_bnd_flux_dir = real3dk("sw_bnd_flux_dir", ncol, nlay+1, nbnd);
  auto sw_bnd_flux_dif = real3dk("sw_bnd_flux_dif", ncol, nlay+1, nbnd);
  logger->info("Populate band-resolved 3d fluxes for test case with only transition band flux...\n");
  Kokkos::parallel_for(MDRP::template get<3>({ncol, nlay+1, nbnd}), KOKKOS_LAMBDA(int icol, int ilay, int ibnd) {
    if (ibnd < 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
    } else if (ibnd == 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 1;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 1;
    } else {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
    }
  });
  // Compute surface fluxes
  logger->info("Compute broadband surface fluxes...\n");
  interface_t::compute_broadband_surface_fluxes(
    ncol, kbot, nbnd,
    sw_bnd_flux_dir, sw_bnd_flux_dif,
    sfc_flux_dir_vis, sfc_flux_dir_nir,
    sfc_flux_dif_vis, sfc_flux_dif_nir
  );
  // Check computed surface fluxes
  logger->info("Check computed fluxes...\n");
  const double tol = 1e-10;  // tolerance on floating point inequality for assertions
  REQUIRE(std::abs(chc(sfc_flux_dir_nir)(0) - 0.5) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dir_vis)(0) - 0.5) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_nir)(0) - 0.5) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_vis)(0) - 0.5) < tol);
  // ---------------------------------

  // ---------------------------------
  // Test case, only flux in NIR bands
  logger->info("Populate band-resolved 3d fluxes for test case with only NIR flux...\n");
  Kokkos::parallel_for(MDRP::template get<3>({ncol, nlay+1, nbnd}), KOKKOS_LAMBDA(int icol, int ilay, int ibnd) {
    if (ibnd < 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 1;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 1;
    } else if (ibnd == 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
    } else {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
    }
  });
  // Compute surface fluxes
  logger->info("Compute broadband surface fluxes...\n");
  interface_t::compute_broadband_surface_fluxes(
    ncol, kbot, nbnd,
    sw_bnd_flux_dir, sw_bnd_flux_dif,
    sfc_flux_dir_vis, sfc_flux_dir_nir,
    sfc_flux_dif_vis, sfc_flux_dif_nir
  );
  // Check computed surface fluxes
  logger->info("Check computed fluxes...\n");
  REQUIRE(std::abs(chc(sfc_flux_dir_nir)(0) - 9.0) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dir_vis)(0) - 0.0) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_nir)(0) - 9.0) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_vis)(0) - 0.0) < tol);
  // ---------------------------------

  // ---------------------------------
  // Test case, only flux in VIS bands
  logger->info("Populate band-resolved 3d fluxes for test case with only VIS/UV flux...\n");
  Kokkos::parallel_for(MDRP::template get<3>({ncol, nlay+1, nbnd}), KOKKOS_LAMBDA(int icol, int ilay, int ibnd) {
    if (ibnd < 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
    } else if (ibnd == 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 0;
    } else {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 1;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 1;
    }
  });
  // Compute surface fluxes
  logger->info("Compute broadband surface fluxes...\n");
  interface_t::compute_broadband_surface_fluxes(
    ncol, kbot, nbnd,
    sw_bnd_flux_dir, sw_bnd_flux_dif,
    sfc_flux_dir_vis, sfc_flux_dir_nir,
    sfc_flux_dif_vis, sfc_flux_dif_nir
  );
  // Check computed surface fluxes
  logger->info("Check computed fluxes...\n");
  REQUIRE(std::abs(chc(sfc_flux_dir_nir)(0) - 0.0) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dir_vis)(0) - 4.0) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_nir)(0) - 0.0) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_vis)(0) - 4.0) < tol);
  // ---------------------------------

  // ---------------------------------
  // Test case, only flux in all bands
  logger->info("Populate band-resolved 3d fluxes for test with non-zero flux in all bands...\n");
  Kokkos::parallel_for(MDRP::template get<3>({ncol, nlay+1, nbnd}), KOKKOS_LAMBDA(int icol, int ilay, int ibnd) {
    if (ibnd < 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 1.0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 2.0;
    } else if (ibnd == 9) {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 3.0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 4.0;
    } else {
      sw_bnd_flux_dir(icol,ilay,ibnd) = 5.0;
      sw_bnd_flux_dif(icol,ilay,ibnd) = 6.0;
    }
  });
  // Compute surface fluxes
  logger->info("Compute broadband surface fluxes...\n");
  interface_t::compute_broadband_surface_fluxes(
    ncol, kbot, nbnd,
    sw_bnd_flux_dir, sw_bnd_flux_dif,
    sfc_flux_dir_vis, sfc_flux_dir_nir,
    sfc_flux_dif_vis, sfc_flux_dif_nir
  );
  // Check computed surface fluxes
  logger->info("Check computed fluxes...\n");
  REQUIRE(std::abs(chc(sfc_flux_dir_nir)(0) - 10.5) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dir_vis)(0) - 21.5) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_nir)(0) - 20.0) < tol);
  REQUIRE(std::abs(chc(sfc_flux_dif_vis)(0) - 26.0) < tol);
  // ---------------------------------

  // Finalize YAKL
  logger->info("Free memory...\n");
  interface_t::rrtmgp_finalize();
  gas_concs.reset();
  pool_t::finalize();
  scream::finalize_kls();
}

TEST_CASE("rrtmgp_test_radiation_do_k") {
  // If we specify rad every step, radiation_do should always be true
  REQUIRE(scream::rrtmgp::radiation_do(1, 0) == true);
  REQUIRE(scream::rrtmgp::radiation_do(1, 1) == true);
  REQUIRE(scream::rrtmgp::radiation_do(1, 2) == true);

    // Test cases where we want rad called every other step
  REQUIRE(scream::rrtmgp::radiation_do(2, 0) == true);
  REQUIRE(scream::rrtmgp::radiation_do(2, 1) == false);
  REQUIRE(scream::rrtmgp::radiation_do(2, 2) == true);
  REQUIRE(scream::rrtmgp::radiation_do(2, 3) == false);

  // Test cases where we want rad every third step
  REQUIRE(scream::rrtmgp::radiation_do(3, 0) == true);
  REQUIRE(scream::rrtmgp::radiation_do(3, 1) == false);
  REQUIRE(scream::rrtmgp::radiation_do(3, 2) == false);
  REQUIRE(scream::rrtmgp::radiation_do(3, 3) == true);
  REQUIRE(scream::rrtmgp::radiation_do(3, 4) == false);
  REQUIRE(scream::rrtmgp::radiation_do(3, 5) == false);
  REQUIRE(scream::rrtmgp::radiation_do(3, 6) == true);
}

TEST_CASE("rrtmgp_test_check_range_k") {
  // Initialize YAKL
  scream::init_kls();
  pool_t::init(10000);
  // Create some dummy data and test with both values inside valid range and outside
  auto dummy = real2dk("dummy", 2, 1);
  // All values within range
  Kokkos::deep_copy(dummy, 0.1);
  REQUIRE(scream::rrtmgp::check_range_k(dummy, 0.0, 1.0, "dummy") == true);
  // At least one value below lower bound
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int i) {dummy(i, 0) = -0.1;});
  REQUIRE(scream::rrtmgp::check_range_k(dummy, 0.0, 1.0, "dummy") == false);
  // At least one value above upper bound
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int i) {dummy(i, 0) = 1.1;});
  REQUIRE(scream::rrtmgp::check_range_k(dummy, 0.0, 1.0, "dummy") == false);
  pool_t::finalize();
  scream::finalize_kls();
}

TEST_CASE("rrtmgp_test_subcol_gen_k") {
  // Initialize YAKL
  scream::init_kls();
  pool_t::init(10000);
  // Create dummy data
  const int ncol = 1;
  const int nlay = 4;
  const int ngpt = 10;
  auto cldfrac = real2dk("cldfrac", ncol, nlay);
  // Set cldfrac values
  Kokkos::deep_copy(cldfrac, 0.0);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
      cldfrac(0,0) = 1;
      cldfrac(0,1) = 0.5;
      cldfrac(0,2) = 0;
      cldfrac(0,3) = 1;
  });
  auto cldmask = int3dk("cldmask", ncol, nlay, ngpt);
  auto cldfrac_from_mask = real2dk("cldfrac_from_mask", ncol, nlay);
  // Run subcol gen, make sure we get what we expect; do this for some different seed values
  for (unsigned seed = 0; seed < 10; seed++) {
    auto seeds = int1dk("seeds", ncol);
    Kokkos::deep_copy(seeds, seed);
    interface_t::get_subcolumn_mask(ncol, nlay, ngpt, cldfrac, 1, seeds, cldmask);
    // Check answers by computing new cldfrac from mask
    Kokkos::deep_copy(cldfrac_from_mask, 0.0);
    Kokkos::parallel_for(MDRP::template get<2>({ncol, nlay}), KOKKOS_LAMBDA(int icol, int ilay) {
      for (int igpt = 0; igpt < ngpt; ++igpt) {
        real cldmask_real = cldmask(icol,ilay,igpt);
        cldfrac_from_mask(icol,ilay) += cldmask_real;
      }
    });
    Kokkos::parallel_for(MDRP::template get<2>({ncol, nlay}), KOKKOS_LAMBDA(int icol, int ilay) {
      cldfrac_from_mask(icol,ilay) = cldfrac_from_mask(icol,ilay) / ngpt;
    });
    // For cldfrac 1 we should get 1, for cldfrac 0 we should get 0, but in between we cannot be sure
    // deterministically, since the computed cloud mask depends on pseudo-random numbers
    REQUIRE(chc(cldfrac_from_mask)(0,0) == 1);
    REQUIRE(chc(cldfrac_from_mask)(0,1) <= 1);
    REQUIRE(chc(cldfrac_from_mask)(0,2) == 0);
    REQUIRE(chc(cldfrac_from_mask)(0,3) == 1);
  }

  // For maximum-random overlap, vertically-contiguous layers maximimally overlap,
  // thus if we have non-zero cloud fraction in two adjacent layers, then every subcolumn
  // that has cloud in the layer above must also have cloud in the layer below; test
  // this property by creating two layers with non-zero cloud fraction, creating subcolums,
  // and verifying that every subcolumn with cloud in layer 1 has cloud in layer 2
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldfrac(0,0) = 0.5;
    cldfrac(0,1) = 0.5;
    cldfrac(0,2) = 0;
    cldfrac(0,3) = 0;
  });
  for (unsigned seed = 0; seed < 10; seed++) {
    auto seeds = int1dk("seeds", ncol);
    Kokkos::deep_copy(seeds, seed);
    interface_t::get_subcolumn_mask(ncol, nlay, ngpt, cldfrac, 1, seeds, cldmask);
    auto cldmask_h = chc(cldmask);
    for (int igpt = 0; igpt < ngpt; igpt++) {
      if (cldmask_h(0,0,igpt) == 1) {
        REQUIRE(cldmask_h(0,1,igpt) == 1);
      }
    }
  }
  // Clean up after test
  pool_t::finalize();
  scream::finalize_kls();
}

TEST_CASE("rrtmgp_cloud_area_k") {
  // Initialize YAKL
  scream::init_kls();
  pool_t::init(10000);
  // Create dummy data
  const int ncol = 1;
  const int nlay = 2;
  const int ngpt = 3;
  auto cldtau = real3dk("cldtau", ncol, nlay, ngpt);
  auto cldtot = real1dk("cldtot", ncol);
  auto pmid = real2dk("pmid", ncol, nlay);

  // Set up pressure levels for test problem
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    pmid(0,0) = 100;
    pmid(0,1) = 200;
  });

  // Case:
  //
  // 0 0 0
  // 0 0 0
  //
  // should give cldtot = 0.0
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldtau(0,0,0) = 0;
    cldtau(0,0,1) = 0;
    cldtau(0,0,2) = 0;
    cldtau(0,1,0) = 0;
    cldtau(0,1,1) = 0;
    cldtau(0,1,2) = 0;
  });
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 0.0);

  // Case:
  //
  // 1 1 1
  // 1 1 1
  //
  // should give cldtot = 1.0
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldtau(0,0,0) = 1;
    cldtau(0,0,1) = 1;
    cldtau(0,0,2) = 1;
    cldtau(0,1,0) = 1;
    cldtau(0,1,1) = 1;
    cldtau(0,1,2) = 1;
  });
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 1.0);

  // Case:
  //
  // 1 1 0  100
  // 0 0 1  200
  //
  // should give cldtot = 1.0
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
      cldtau(0,0,0) = 0.1;
      cldtau(0,0,1) = 1.5;
      cldtau(0,0,2) = 0;
      cldtau(0,1,0) = 0;
      cldtau(0,1,1) = 0;
      cldtau(0,1,2) = 1.0;
  });
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 1.0);
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 0, 150, pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 2.0 / 3.0);
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 110, 250, pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 1.0 / 3.0);

  // Case:
  //
  // 1 0 0
  // 1 0 1
  //
  // should give cldtot = 2/3
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldtau(0,0,0) = 1;
    cldtau(0,0,1) = 0;
    cldtau(0,0,2) = 0;
    cldtau(0,1,0) = 1;
    cldtau(0,1,1) = 0;
    cldtau(0,1,2) = 1;
  });
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 0, std::numeric_limits<scream::Real>::max(), pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 2.0 / 3.0);
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 0, 100, pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 0.0);
  interface_t::compute_cloud_area(ncol, nlay, ngpt, 100, 300, pmid, cldtau, cldtot);
  REQUIRE(chc(cldtot)(0) == 2.0 / 3.0);
  pool_t::finalize();
  scream::finalize_kls();
}

TEST_CASE("rrtmgp_aerocom_cloudtop_k") {
  // Initialize YAKL
  scream::init_kls();
  pool_t::init(10000);

  // Create dummy data
  const int ncol = 1;
  const int nlay = 9;
  // Set up input fields
  auto tmid        = real2dk("tmid", ncol, nlay);
  auto pmid        = real2dk("pmid", ncol, nlay);
  auto p_del       = real2dk("p_del", ncol, nlay);
  auto z_del       = real2dk("z_del", ncol, nlay);
  auto qc          = real2dk("qc", ncol, nlay);
  auto qi          = real2dk("qi", ncol, nlay);
  auto rel         = real2dk("rel", ncol, nlay);
  auto rei         = real2dk("rei", ncol, nlay);
  auto cldfrac_tot = real2dk("cldfrac_tot", ncol, nlay);
  auto nc          = real2dk("nc", ncol, nlay);
  // Set up output fields
  auto tmid_at_cldtop          = real1dk("tmid_at_cldtop", ncol);
  auto pmid_at_cldtop          = real1dk("pmid_at_cldtop", ncol);
  auto cldfrac_ice_at_cldtop   = real1dk("cldfrac_ice_at_cldtop", ncol);
  auto cldfrac_liq_at_cldtop   = real1dk("cldfrac_liq_at_cldtop", ncol);
  auto cldfrac_tot_at_cldtop   = real1dk("cldfrac_tot_at_cldtop", ncol);
  auto cdnc_at_cldtop          = real1dk("cdnc_at_cldtop", ncol);
  auto eff_radius_qc_at_cldtop = real1dk("eff_radius_qc_at_cldtop", ncol);
  auto eff_radius_qi_at_cldtop = real1dk("eff_radius_qi_at_cldtop", ncol);

  // Case 1: if no clouds, everything goes to zero
  Kokkos::deep_copy(tmid, 300.0);
  Kokkos::deep_copy(pmid, 100.0);
  Kokkos::deep_copy(p_del, 10.0);
  Kokkos::deep_copy(z_del, 100.0);
  Kokkos::deep_copy(qc, 1.0);
  Kokkos::deep_copy(qi, 1.0);
  Kokkos::deep_copy(cldfrac_tot, 0.0);
  Kokkos::deep_copy(nc, 5.0);
  Kokkos::deep_copy(rel, 10.0);
  Kokkos::deep_copy(rei, 10.0);
  // Call the function
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  // Check the results
  REQUIRE(chc(tmid_at_cldtop)(0) == 0.0);
  REQUIRE(chc(pmid_at_cldtop)(0) == 0.0);
  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) == 0.0);
  REQUIRE(chc(cldfrac_liq_at_cldtop)(0) == 0.0);
  REQUIRE(chc(cldfrac_ice_at_cldtop)(0) == 0.0);
  REQUIRE(chc(cdnc_at_cldtop)(0) == 0.0);
  REQUIRE(chc(eff_radius_qc_at_cldtop)(0) == 0.0);
  REQUIRE(chc(eff_radius_qi_at_cldtop)(0) == 0.0);

  // Case 2: if all clouds, everything goes to 1 * its value
  Kokkos::deep_copy(cldfrac_tot, 1.0);
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(chc(tmid_at_cldtop)(0) == 300.0);
  REQUIRE(chc(pmid_at_cldtop)(0) == 100.0);
  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) == 1.0);
  REQUIRE(chc(cldfrac_liq_at_cldtop)(0) == 0.5);
  REQUIRE(chc(cldfrac_ice_at_cldtop)(0) == 0.5);
  REQUIRE(chc(cdnc_at_cldtop)(0) > 0.0);
  REQUIRE(chc(eff_radius_qc_at_cldtop)(0) > 0.0);
  REQUIRE(chc(eff_radius_qi_at_cldtop)(0) > 0.0);

  // Case 3: test max overlap (if contiguous cloudy layers, then max)
  Kokkos::deep_copy(cldfrac_tot, 0.0);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldfrac_tot(0, 1) = 0.5;
    cldfrac_tot(0, 2) = 0.7;
    cldfrac_tot(0, 3) = 0.3;
    cldfrac_tot(0, 4) = 0.2;
  });
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) == .7);

  // Case 3xtra: test max overlap
  // This case produces >0.7 due to slight enhancement in the presence of a
  // local minimum (0.1 is the local minimum between 0.2 and 0.4)
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldfrac_tot(0, 4) = 0.1;
    cldfrac_tot(0, 5) = 0.4;
    cldfrac_tot(0, 6) = 0.2;
  });
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) > .7);

  // Case 4: test random overlap (if non-contiguous cloudy layers, then
  // random)
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldfrac_tot(0, 4) = 0.0;
    cldfrac_tot(0, 5) = 0.1;
  });
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) > .7);  // larger than the max

  // Case 5a: test independence of ice and liquid fractions
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldfrac_tot(0, 1) = 1.0;
    cldfrac_tot(0, 6) = 1.0;
    cldfrac_tot(0, 7) = 0.2;
  });
  Kokkos::deep_copy(qc, 1.0);
  Kokkos::deep_copy(qi, 0.0);
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) == 1.0);
  REQUIRE(chc(cldfrac_liq_at_cldtop)(0) == 1.0);
  REQUIRE(chc(cldfrac_ice_at_cldtop)(0) == 0.0);

  // Case 5b: test independence of ice and liquid fractions
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldfrac_tot(0, 1) = 1.0;
    cldfrac_tot(0, 6) = 1.0;
    cldfrac_tot(0, 7) = 0.2;
  });
  Kokkos::deep_copy(qc, 0.0);
  Kokkos::deep_copy(qi, 1.0);
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);

  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) == 1.0);
  REQUIRE(chc(cldfrac_liq_at_cldtop)(0) == 0.0);
  REQUIRE(chc(cldfrac_ice_at_cldtop)(0) == 1.0);

  // Case 6: test independence of ice and liquid fractions
  // There is NOT complete independence...
  // Essentially, higher ice clouds mask lower liquid clouds
  // This can be problematic if the ice clouds are thin...
  // We will revisit and validate this assumption later
  Kokkos::deep_copy(cldfrac_tot, 0.0);
  Kokkos::deep_copy(qc, 0.0);
  Kokkos::deep_copy(qi, 0.0);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int /* dummy */) {
    cldfrac_tot(0, 1) = 0.5;  // ice
    cldfrac_tot(0, 2) = 0.7;  // ice ------> max
    cldfrac_tot(0, 3) = 0.3;  // ice
    // note cldfrac_tot(1, 5) is 0
    cldfrac_tot(0, 5) = 0.2;  // liq
    cldfrac_tot(0, 6) = 0.5;  // liq ------> not max
    cldfrac_tot(0, 7) = 0.1;  // liq
    // note cldfrac_tot(1, 9) is 0
    qi(0, 1) = 100;
    qi(0, 2) = 200;
    qi(0, 3) = 50;
    // note qc(1, 5) is 0
    // note qi(1, 5) is 0
    qc(0, 5) = 20;
    qc(0, 6) = 50;
    qc(0, 7) = 10;
  });
  interface_t::compute_aerocom_cloudtop(
    ncol, nlay, tmid, pmid, p_del, z_del, qc, qi, rel, rei, cldfrac_tot, nc,
    tmid_at_cldtop, pmid_at_cldtop, cldfrac_ice_at_cldtop,
    cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
    eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);
  REQUIRE(chc(cldfrac_tot_at_cldtop)(0) > 0.70);  // unaffected
  REQUIRE(chc(cldfrac_liq_at_cldtop)(0) < 0.50);  // not max
  REQUIRE(chc(cldfrac_ice_at_cldtop)(0) == 0.7);  // max

  // cleanup
  pool_t::finalize();
  scream::finalize_kls();
}
#endif

}
