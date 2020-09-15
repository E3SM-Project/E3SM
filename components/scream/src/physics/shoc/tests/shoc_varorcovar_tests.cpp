#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocVarorCovar {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 4;
    static constexpr auto nlevi   = nlev + 1;

    //NOTE: This routine does not compute the (co)variance
    // for boundary points.  Input Grid values that are never used
    // are set to ZERO so as not to confuse the test data
    // that is actually used (and is also an implicit test to be
    // sure that operations are not done at the boundary).
    // Output values at boundary set to something other than
    // zero to check that they are not modified.

    //IN: tunefac, isotropy_zi, tkh_zi, dz_zi,
    //IN: invar1, invar2
    //OUT: varorcovar

    //********************************************
    // TEST ONE
    // Verify the variance calculation is valid

    // Define data.  Some of this is recycled for tests two and three.

    // Define delta z on the nlevi grid [m]
    static constexpr Real dz_zi[nlevi] = {0.0, 100.0, 50.0, 20.0, 0.0};
    // Eddy coefficients on the nlevi grid [m2s-1]
    static constexpr Real tkh_zi[nlevi] = {0.0, 10.0, 10.0, 10.0, 0.0};
    // Return to isotropy timescale [s]
    static constexpr Real isotropy_zi[nlevi] = {0.0, 100.0, 100.0, 100.0, 0.0};
    // Invar (in this example we use potential temperature) [K]
    //   this variable should be on the nlev grid
    static constexpr Real invar_theta[nlev] = {320.0, 315.0, 315.0, 316.0};
    // Initialize (co)variance on the nlevi grid, set boundaries
    //   to large values to make sure they are not modified
    Real varorcovar[nlevi] = {100., 0., 0., 0., 100.};
    // Tuning factor
    static constexpr Real tunefac=1.0;

    // Initialzie data structure for bridgeing to F90
    SHOCVarorcovarData SDS(shcol, nlev, nlevi, tunefac);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev && SDS.nlevi() && SDS.tunefac == tunefac) );
    REQUIRE(nlevi - nlev == 1);
    REQUIRE(shcol > 0);

    // Fill in test data
    for(Int s = 0; s < shcol; ++s) {
      // First on the nlev grid
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	// Fill invar1 and invar2 with the SAME
	//  variable for this case to test the
	//  variance calculation of this function
	SDS.invar1[offset] = invar_theta[n];
	SDS.invar2[offset] = invar_theta[n];
      }

      // Now for data on the nlevi grid
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset = n + s * nlevi;

	SDS.tkh_zi[offset] = tkh_zi[n];
	SDS.isotropy_zi[offset] = isotropy_zi[n];
	SDS.dz_zi[offset] = dz_zi[n];
	SDS.varorcovar[offset] = varorcovar[n];
      }
    }

    // Check that the inputs make sense

    REQUIRE(SDS.tunefac > 0.0);
    // Check to make sure that dz_zi, tkh_zi, and isotropy_zi
    //  (outside of the boundaries) are greater than zero
    for(Int s = 0; s < shcol; ++s) {
      // do NOT check boundaries!
      for(Int n = 1; n < nlevi-1; ++n) {
	const auto offset = n + s * nlevi;
	REQUIRE(SDS.dz_zi[offset] > 0.0);
	REQUIRE(SDS.tkh_zi[offset] > 0.0);
	REQUIRE(SDS.isotropy_zi[offset] > 0.0);
      }
      // For this test make sure that invar1 = invar2
      for(Int n = 0; n < nlev; ++n){
	const auto offset = n + s * nlev;
	REQUIRE(SDS.invar1[offset] == SDS.invar2[offset]);
      }
    }

    // Call the fortran implementation for variance
    calc_shoc_varorcovar(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset = n + s * nlevi;

	// validate that the boundary points have NOT been modified
	if (n == 0 || n == nlevi){
          REQUIRE(SDS.varorcovar[offset] == 100.0);
	}
	else{

	// Validate that all values are greater to
	//   or equal to zero

          REQUIRE(SDS.varorcovar[offset] >= 0.0);

          // well mixed layer test
          if ((invar_theta[n-1] - invar_theta[n]) == 0.0){
            REQUIRE(SDS.varorcovar[offset] == 0.0);
          }
	}
      }
    }

    //********************************************
    // TEST TWO
    // Verify the Covariance calculation is valid

    // Now we check the covariance calculation
    // Define an array to be total water [g/kg]
    Real invar_qw[nlev] = {5.0, 10.0, 15.0, 20.0};

    // NOTE: all other inputs are reused from test one

    // convert total water to [kg/kg]
    for (Int n = 0; n < nlev; ++n){
      invar_qw[n] = invar_qw[n]/1000.;
    }

    // Update invar2 to be total water
    for(Int s = 0; s < shcol; ++s) {
      // First on the nlev grid
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.invar2[offset] = invar_qw[n];
      }
    }

    // Check inputs
    for(Int s = 0; s < shcol; ++s) {
      // For this test make sure that invar1 is NOT
      //  equal to invar2
      for(Int n = 0; n < nlev; ++n){
	const auto offset = n + s * nlev;
	REQUIRE(SDS.invar1[offset] != SDS.invar2[offset]);
      }
    }

    // Call the fortran implementation for covariance
    calc_shoc_varorcovar(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset = n + s * nlevi;

	// validate that the boundary points
	//   have NOT been modified
	if (n == 0 || n == nlevi){
          REQUIRE(SDS.varorcovar[offset] == 100.);
	}
	else{

          // well mixed layer test
          if ((invar_theta[n-1] - invar_theta[n]) == 0.0 ||
	      (invar_qw[n-1] - invar_qw[n]) == 0.0){
            REQUIRE(SDS.varorcovar[offset] == 0.0);
          }

	  // validate values are NEGATIVE if potential
	  //  temperature INCREASES with height and total water
	  //  DECREASES with height
          if ((invar_theta[n-1] - invar_theta[n]) > 0.0 &&
	      (invar_qw[n-1] - invar_qw[n]) < 0.0){
            REQUIRE(SDS.varorcovar[offset] < 0.0);
          }

	  // validate values are POSITIVE if both
	  //   potential temperature and total water
	  //   DECREASE with height
          if ((invar_theta[n-1] - invar_theta[n]) < 0.0 &&
	      (invar_qw[n-1] - invar_qw[n]) < 0.0){
            REQUIRE(SDS.varorcovar[offset] > 0.0);
          }

	}
      }
    }

    //********************************************
    // TEST THREE
    // Verify that vertical differencing is done correctly

    //  Assume that isotropy and tkh are constant with height
    //  Assign values of potential temperature that are increasing
    //  at a constant rate per GRID box.  Assign dz that INCREASES
    //  with height and check that varorcovar DECREASES with height

    static constexpr Real invar_th2[nlev] = {320.0, 315.0, 310.0, 305.0};

    // Update invar1 and invar2 to be identical
    for(Int s = 0; s < shcol; ++s) {
      // First on the nlev grid
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	// Set both inputs to the same value
	SDS.invar1[offset] = invar_th2[n];
	SDS.invar2[offset] = invar_th2[n];
      }
    }

    // Check the inputs
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 1; n < nlevi-1; ++n) {
	const auto offset = n + s * nlevi;
	// Validate that values of dz_zi are INCREASING with height
	REQUIRE(SDS.dz_zi[offset]-SDS.dz_zi[offset+1] > 0.0);
      }
      // For this test make sure that invar1 = invar2
      for(Int n = 0; n < nlev; ++n){
	const auto offset = n + s * nlev;
	REQUIRE(SDS.invar1[offset] == SDS.invar2[offset]);
      }
    }

    // Call the fortran implementation for variance
    calc_shoc_varorcovar(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 1; n < nlev-1; ++n) {
	const auto offset = n + s * nlevi;

	// Validate that values of varorcovar
	//  are decreasing with height
	REQUIRE(SDS.varorcovar[offset]-SDS.varorcovar[offset+1] < 0.0);

      }
    }
  }

static void run_bfb()
  {
    SHOCVarorcovarData SDS_f90[] = {
      //               shcol, nlev, nlevi, tunefac
      SHOCVarorcovarData(10, 71, 72, 1),
      SHOCVarorcovarData(10, 12, 13, 1),
      SHOCVarorcovarData(7,  16, 17, 1),
      SHOCVarorcovarData(2, 7, 8, 0.005),
    };

    static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(SHOCVarorcovarData);

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    SHOCVarorcovarData SDS_cxx[] = {
      SHOCVarorcovarData(SDS_f90[0]),
      SHOCVarorcovarData(SDS_f90[1]),
      SHOCVarorcovarData(SDS_f90[2]),
      SHOCVarorcovarData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      calc_shoc_varorcovar(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::util::TransposeDirection::c2f>();
      // expects data in fortran layout
      calc_shoc_varorcovar_f(d.shcol(), d.nlev(), d.nlevi(),
                             d.tunefac, d.isotropy_zi,
                             d.tkh_zi, d.dz_zi,
                             d.invar1, d.invar2, d.varorcovar);
      d.transpose<ekat::util::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      SHOCVarorcovarData& d_f90 = SDS_f90[i];
      SHOCVarorcovarData& d_cxx = SDS_cxx[i];
      for (Int k = 0; k < d_f90.total1x3(); ++k) {
        REQUIRE(d_f90.varorcovar[k] == d_cxx.varorcovar[k]);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_varorcovar_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocVarorCovar;

  TestStruct::run_property();
}

TEST_CASE("shoc_varorcovar_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocVarorCovar;

  TestStruct::run_bfb();
}

} // namespace
