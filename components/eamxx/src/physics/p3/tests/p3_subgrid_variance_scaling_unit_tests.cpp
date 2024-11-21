#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iomanip>      // std::setprecision

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3SubgridVarianceScaling : public UnitWrap::UnitTest<D>::Base
{

  //-----------------------------------------------------------------
  void run_bfb_tests() {
    //test that C++ and F90 implementations are BFB

    //Set of relvar values to loop over
    Scalar relvars[16] = {
      0.1,0.5,1.0,2.0,
      3.0,4.0,5.0,6.0,
      6.5,7.0,8.0,9.0,
      9.1,9.5,9.8,10.};

    //Set of exponents to loop over
    Scalar expons[3] = {1.0,2.47,0.1};

    Scalar baseline_scaling;

    //Make C++ output available on host and device
    view_1d<Scalar> scaling_device("c scaling",1);
    auto scaling_host = Kokkos::create_mirror_view(scaling_device);

    for (Int i = 0; i < 3; ++i) {  // loop over exponents
      for (Int j = 0; j < 16; ++j) { // loop over relvars

	// Get baseline solution
	// ----------------------------------
        if (this->m_baseline_action == COMPARE) {
          ekat::read(&baseline_scaling, 1, Base::m_fid);
        }

	// Get C++ solution
	// ----------------------------------

	//Make scalar copies so available on device
	Scalar expon = expons[i];
	Scalar relvar = relvars[j];

	RangePolicy my_policy(0,1);
	Kokkos::parallel_for(my_policy,KOKKOS_LAMBDA(int /* i */){
	    Spack scalings = Functions::subgrid_variance_scaling(Spack(relvar),expon );

	    //all elements of scalings are identical. just copy 1 back to host.
	    scaling_device(0) = scalings[0];
	  });

	// Copy results back to host
	Kokkos::deep_copy(scaling_host, scaling_device);

	// Validate results
        if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
          REQUIRE(baseline_scaling == scaling_host(0) );
        }
        else if (this->m_baseline_action == GENERATE) {
          ekat::write(&scaling_host(0), 1, Base::m_fid);
        }
      } //end loop over relvar[j]
    } //end loop over expons[i]
  } //end function run_bfb_tests

  //-----------------------------------------------------------------
  KOKKOS_FUNCTION static void subgrid_variance_scaling_linearity_test(const Scalar& relvar,
    int& errors) {
    //If expon=1, subgrid_variance_scaling should be 1

    Scalar tol = C::macheps * 1e3; //1e3 is scale factor to make pass, essentially an estimate of numerical error

    //Get value from C++ code
    const Spack relvars(relvar);
    Spack c_scaling = Functions::subgrid_variance_scaling(relvars,1.0);

    if ( std::abs(c_scaling[0] -  1) > tol ){
      printf("subgrid_variance_scaling should be 1 for expon=1, but is %e. "
	     "Diff = %e, Tol = %e\n",c_scaling[0],c_scaling[0]-1, tol);
	errors++;}
  }

  //-----------------------------------------------------------------
  KOKKOS_FUNCTION static void subgrid_variance_scaling_relvar1_test(int& errors) {
    //If relvar=1, subgrid_variance_scaling should be factorial(expon)

    Scalar tol = C::macheps * 1e3; //1e3 is scale factor to make pass, essentially an estimate of numerical error

    //Get value from C++ code
    const Spack ones(1);
    Spack c_scaling = Functions::subgrid_variance_scaling(ones,4.0);

    Real fact = std::tgamma(5.0); //factorial(n) = gamma(n+1)

    if ( std::abs(c_scaling[0] -  fact) > tol ){
      printf("subgrid_variance_scaling should be factorial(expon) when relvar=1. "
	     "For expon=4, should be %f but is=%f\n Diff = %e, Tol = %e\n",
	     fact,c_scaling[0], c_scaling[0] -  fact, tol);
      errors++;}
  }

  //-----------------------------------------------------------------
  KOKKOS_FUNCTION static void subgrid_variance_scaling_relvar3_test(int& errors) {
  //If expon=3, subgrid variance scaling should be relvar^3+3*relvar^2+2*relvar/relvar^3

  Scalar tol = C::macheps * 100; //100 is a fudge factor to make sure tests pass. 10 was too small for gnu on CPU.

  Real relvar_info[max_pack_size] = {0.1,0.5,1.0,2.0,
                                     3.0,4.0,5.0,6.0,
                                     6.5,7.0,8.0,9.0,
                                     9.1,9.5,9.8,10.};

  for (Int s = 0; s < 16; ++s) {
    Spack relvars=Spack(relvar_info[s]);

    //Get value from C++ code
    Spack c_scaling = Functions::subgrid_variance_scaling(relvars,3.0);

    //Get analytic expected value
    Real targ=1+3/relvar_info[s] + 2/std::pow(relvar_info[s],2.0);

    //Expected relative discrepancy is relative condition # * tolerance
    //For expon=3, expected val is 1+3/relvar + 2/relvar**2.
    //Condition number is x*f'(x)/f(x) = (3*relvar + 4)/(relvar**2. + 3*relvar+2)
    const Real cond_num = (3.*relvar_info[s] + 4.)/(std::pow(relvar_info[s],2.0) +3*relvar_info[s]+2.0);
    const Real max_tol = tol*cond_num;

    if ( std::abs(targ - c_scaling[0]) > max_tol * targ ){
      printf("When expon=3, subgrid_variance_scaling doesn't match analytic expectation. "
	     "Val = %e, expected = %e, rel diff = %e, tol = %e\n",
	     c_scaling[0],targ, (targ-c_scaling[0]), max_tol*targ );
      errors++;
    } // end if
  }   //end for
  }   //end relvar3_test

  //-----------------------------------------------------------------
  void run_property_tests() {
    /*This function executes all the SGS variance scaling tests by looping
     *over a bunch of test and summing their return statuses.
     *If that sum is zero, no errors have occurred. Otherwise you have errors.
     *We do that loop in the parallel reduce below.
     */
    int nerr = 0;

    //functions below use Spack size <16 but can't deal w/ exceptions on GPU, so do it here.
    TeamPolicy policy(ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("SGSvarScaling::run", policy,
      KOKKOS_LAMBDA(const MemberType& /* team */, int& errors) {
        errors = 0;

        //If expon=1, subgrid_variance_scaling should be 1
        //                            args = relvar,return error count
        subgrid_variance_scaling_linearity_test(10.,errors);
        subgrid_variance_scaling_linearity_test(0.1,errors);

        //If relvar=1, subgrid_variance_scaling should be factorial(expon)
        //                            args = return error count
        subgrid_variance_scaling_relvar1_test(errors);

        //If expon=3, subgrid variance scaling should be relvar^3+3*relvar^2+2*relvar/relvar^3
        //                          args = return error count
        subgrid_variance_scaling_relvar3_test(errors);
      }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);
  } //end of TestP3SubgridVarianceScaling struct

}; // UnitWrap

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_subgrid_variance_scaling_test", "[p3_subgrid_variance_scaling_test]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3SubgridVarianceScaling;

  T t;
  t.run_bfb_tests();
  t.run_property_tests();
}

} // namespace
