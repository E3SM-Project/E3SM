#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestIceDepositionSublimation {

  static void run_property(){
    //Note that a lot of property tests are included in run_bfb for simplicity

    //Choose default values (just grabbed a row from the array in run_bfb)
    Spack qi_incld =  5.1000E-03;
    Spack ni_incld=4.7790E+05;
    Spack T_atm=2.7292E+02;
    Spack qv_sat_l=4.5879E-03;
    Spack qv_sat_i=4.5766E-03;
    Spack abi=2.0649E+00;
    Spack qv=5.0000E-03;
    Scalar inv_dt=1.666667e-02;

    //init output vars
    Spack qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend;
    
    //CHECK THAT UNREASONABLY LARGE VAPOR DEPOSITION DOESN'T LEAVE QV SUBSATURATED WRT QI
    Spack epsi_tmp=1e6; //make 1/(sat removal timescale) huge so vapdep rate removes all supersat in 1 dt.
    Functions::ice_deposition_sublimation(qi_incld, ni_incld, T_atm, qv_sat_l, qv_sat_i,
	        epsi_tmp, abi, qv, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend);
    REQUIRE( (qv2qi_vapdep_tend[0]==0 || std::abs( qv2qi_vapdep_tend[0] - (qv[0] - qv_sat_i[0])*inv_dt) <1e-8) );
   
    //CHECK THAT HUGE SUBLIMATION DOESN'T LEAVE QV SUPERSATURATED WRT QI
    Spack qv_sat_i_tmp=1e-2;
    Functions::ice_deposition_sublimation(qi_incld, ni_incld, T_atm, qv_sat_l, qv_sat_i_tmp,
	        epsi_tmp, abi, qv, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend);
    REQUIRE( (qi2qv_sublim_tend[0]==0 || std::abs( qi2qv_sublim_tend[0] - (qv_sat_i_tmp[0] - qv[0])*inv_dt) <1e-8) );

    //CHECK BEHAVIOR AS DT->0?

  }
  
  static void run_bfb()
  {
    IceDepositionSublimationData f90_data[max_pack_size] = {
      {1.0000E-04,4.5010E+05,2.8750E+02,1.1279E-02,1.1279E-02,0.0000E+00,3.3648E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.5370E+05,2.8542E+02,9.9759E-03,9.9759E-03,0.0000E+00,3.1223E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.5742E+05,2.8334E+02,8.8076E-03,8.8076E-03,0.0000E+00,2.9014E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.6125E+05,2.8125E+02,7.7615E-03,7.7615E-03,0.0000E+00,2.7005E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.6521E+05,2.7917E+02,6.8265E-03,6.8265E-03,0.0000E+00,2.5180E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.6930E+05,2.7709E+02,5.9921E-03,5.9921E-03,0.0000E+00,2.3526E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.7353E+05,2.7501E+02,5.2488E-03,5.2488E-03,0.0000E+00,2.2028E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.7790E+05,2.7292E+02,4.5879E-03,4.5766E-03,6.2108E-02,2.0649E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.8241E+05,2.7084E+02,4.0015E-03,3.9112E-03,6.1911E-02,1.9241E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.8709E+05,2.6876E+02,3.4821E-03,3.3349E-03,6.1708E-02,1.8002E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.9193E+05,2.6667E+02,3.0231E-03,2.8368E-03,6.1502E-02,1.6914E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,4.9695E+05,2.6459E+02,2.6183E-03,2.4074E-03,6.1290E-02,1.5960E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,5.0216E+05,2.6251E+02,2.2621E-03,2.0379E-03,6.1073E-02,1.5125E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,5.0756E+05,2.6042E+02,1.9495E-03,1.7207E-03,6.0850E-02,1.4397E+00,5.0000E-03,1.666667e-02},
      {5.1000E-03,5.1317E+05,2.5834E+02,1.6757E-03,1.4491E-03,6.0620E-02,1.3763E+00,5.0000E-03,1.666667e-02},
      {5.0000E-08,5.4479E+05,2.4793E+02,7.5430E-04,5.8895E-04,4.6769E-04,1.1661E+00,1.5278E-04,1.666667e-02},
    };
    
    static constexpr Int num_runs = sizeof(f90_data) / sizeof(IceDepositionSublimationData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    //for (auto& d : f90_data) {
    //  d.randomize();
    //}

    // Create copies of data for use by cxx and sync it to device. Needs to happen before fortran calls so that
    // inout data is in original state
    view_1d<IceDepositionSublimationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Get data from fortran
    for (auto& d : f90_data) {
      ice_deposition_sublimation(d);
    }

    // Get data from cxx. Run ice_deposition_sublimation from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Scalar inv_dt;
      Spack abi, epsi, ni_incld, qi_incld, qv, qv_sat_i, qv_sat_l, T_atm;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        abi[s] = cxx_device(vs).abi;
        epsi[s] = cxx_device(vs).epsi;
        ni_incld[s] = cxx_device(vs).ni_incld;
        qi_incld[s] = cxx_device(vs).qi_incld;
        qv[s] = cxx_device(vs).qv;
        qv_sat_i[s] = cxx_device(vs).qv_sat_i;
        qv_sat_l[s] = cxx_device(vs).qv_sat_l;
        T_atm[s] = cxx_device(vs).T_atm;
	inv_dt = cxx_device(vs).inv_dt;
      }

      // Init outputs
      Spack ni_sublim_tend(0), qi2qv_sublim_tend(0), qc2qi_berg_tend(0), qv2qi_vapdep_tend(0);


      Functions::ice_deposition_sublimation(qi_incld, ni_incld, T_atm, qv_sat_l, qv_sat_i, epsi, abi, qv, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend);

      // Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_device(vs).ni_sublim_tend = ni_sublim_tend[s];
        cxx_device(vs).qi2qv_sublim_tend = qi2qv_sublim_tend[s];
        cxx_device(vs).qc2qi_berg_tend = qc2qi_berg_tend[s];
        cxx_device(vs).qv2qi_vapdep_tend = qv2qi_vapdep_tend[s];
      }

    });

    Kokkos::deep_copy(cxx_host, cxx_device);

    for (Int i = 0; i < num_runs; ++i) {
      // Verify BFB results
      IceDepositionSublimationData& d_f90 = f90_data[i];
      IceDepositionSublimationData& d_cxx = cxx_host[i];
      REQUIRE(d_f90.qv2qi_vapdep_tend == d_cxx.qv2qi_vapdep_tend);
      REQUIRE(d_f90.qi2qv_sublim_tend == d_cxx.qi2qv_sublim_tend);
      REQUIRE(d_f90.ni_sublim_tend == d_cxx.ni_sublim_tend);
      REQUIRE(d_f90.qc2qi_berg_tend == d_cxx.qc2qi_berg_tend);

      //MAKE SURE OUTPUT IS WITHIN EXPECTED BOUNDS:
      REQUIRE(d_cxx.qv2qi_vapdep_tend >=0);
      REQUIRE(d_cxx.qi2qv_sublim_tend >=0);
      REQUIRE(d_cxx.ni_sublim_tend >=0);
      REQUIRE(d_cxx.qc2qi_berg_tend >=0);

      //vapdep should only occur when qv>qv_sat_i
      REQUIRE( (d_cxx.qv2qi_vapdep_tend==0 || d_cxx.qv + d_cxx.qv2qi_vapdep_tend*d_cxx.inv_dt >= d_cxx.qv_sat_i) );
      //sublim should only occur when qv<qv_sat_i
      REQUIRE( (d_cxx.qi2qv_sublim_tend==0 || d_cxx.qv + d_cxx.qi2qv_sublim_tend*d_cxx.inv_dt <= d_cxx.qv_sat_i) );

      //if T>frz, berg and vapdep should be 0:
      REQUIRE( (d_cxx.T_atm<C::T_zerodegc || d_cxx.qc2qi_berg_tend==0) );
      REQUIRE( (d_cxx.T_atm<C::T_zerodegc || d_cxx.qv2qi_vapdep_tend==0) );

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("ice_deposition_sublimation_property", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceDepositionSublimation;

  TestStruct::run_property();
}
  
TEST_CASE("ice_deposition_sublimation_bfb", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceDepositionSublimation;

  TestStruct::run_bfb();
}

} // empty namespace
