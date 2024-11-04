#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"

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
struct UnitWrap::UnitTest<D>::TestIceNucleation {

  static void run_ice_nucleation_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    constexpr Scalar inv_dt           = 9.558E-04;

    //Loop over logicals being both true and false. Use boolean:integer equivalence to loop.
    for (bool do_predict_nc : {false, true}) {
      for (bool do_prescribed_CCN : {false, true}) {

	IceNucleationData self[max_pack_size] = {
	  // temp,    inv_rho,   ni,     ni_activated,      qv_supersat_i,      inv_dt, do_predict_nc, do_prescribed_CCN
	  {2.106E+02, 8.852E-01, 0.974E+04, 9.221E+03, 5.100E-01, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.052E+02, 8.852E-01, 0.874E+04, 8.221E+03, 4.100E-01, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.352E+02, 8.900E-01, 0.723E+04, 7.221E+03, 3.100E-01, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.212E+02, 9.900E-01, 0.623E+04, 6.221E+03, 2.100E-01, inv_dt, do_predict_nc, do_prescribed_CCN },

	  {2.251E+02, 0.100E+01, 0.574E+04, 5.221E+03, 1.100E-01, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.710E+02, 0.100E+01, 0.474E+04, 4.221E+03, 8.100E-02, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.052E+02, 0.100E+01, 0.323E+04, 3.221E+03, 4.100E-02, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.502E+02, 0.100E+01, 0.223E+04, 2.221E+03, 2.100E-02, inv_dt, do_predict_nc, do_prescribed_CCN },

	  {2.552E+02, 0.950E+00, 0.150E+04, 9.221E+02, 9.952E-02, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.452E+02, 0.950E+00, 0.974E+03, 8.221E+02, 4.952E-02, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.352E+02, 0.950E+00, 0.823E+03, 7.221E+02, 1.952E-02, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.252E+02, 0.950E+00, 0.723E+03, 6.221E+02, 9.952E-02, inv_dt, do_predict_nc, do_prescribed_CCN },

	  {1.990E+02, 1.069E+00, 0.674E+03, 5.221E+01, 6.952E-01, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.952E+02, 1.069E+00, 0.574E+03, 4.221E+01, 3.952E-01, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.852E+02, 1.069E+00, 0.423E+03, 3.221E+01, 1.952E-01, inv_dt, do_predict_nc, do_prescribed_CCN },
	  {2.702E+02, 1.069E+00, 0.323E+03, 2.221E+01, 9.952E-01, inv_dt, do_predict_nc, do_prescribed_CCN }
	};

	// Run the fortran code
	for (Int i = 0; i < max_pack_size; ++i) {
	  ice_nucleation(self[i]);
	}

	// Sync to device
	KTH::view_1d<IceNucleationData> self_host("self_host", max_pack_size);
	view_1d<IceNucleationData> self_device("self_host", max_pack_size);
	std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
	Kokkos::deep_copy(self_device, self_host);

	// Run the lookup from a kernel and copy results back to host
	Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
	    const Int offset = i * Spack::n;

	    // Init pack inputs
	    Spack temp, inv_rho, ni, ni_activated, qv_supersat_i;
	    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
	      temp[s]            = self_device(vs).temp;
	      inv_rho[s]         = self_device(vs).inv_rho;
	      ni[s]              = self_device(vs).ni;
	      ni_activated[s]    = self_device(vs).ni_activated;
	      qv_supersat_i[s]   = self_device(vs).qv_supersat_i;
	    }
	    // outputs
	    Spack qv2qi_nucleat_tend{0.0};
	    Spack ni_nucleat_tend{0.0};
	    Functions::ice_nucleation(temp, inv_rho, ni, ni_activated, qv_supersat_i, self_device(0).inv_dt, do_predict_nc,
				      do_prescribed_CCN, qv2qi_nucleat_tend, ni_nucleat_tend, physics::P3_Constants<Real>());

	    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
	      self_device(vs).qv2qi_nucleat_tend = qv2qi_nucleat_tend[s];
	      self_device(vs).ni_nucleat_tend = ni_nucleat_tend[s];
	    }
	  });

	Kokkos::deep_copy(self_host, self_device);

        if (SCREAM_BFB_TESTING) {
          for (Int s = 0; s < max_pack_size; ++s) {
            REQUIRE(self[s].qv2qi_nucleat_tend == self_host(s).qv2qi_nucleat_tend);
            REQUIRE(self[s].ni_nucleat_tend    == self_host(s).ni_nucleat_tend);
          }
        }
      } //end for do_predict_nc
    } //end for do_prescribed_CCN
  }

  static void run_ice_nucleation_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_ice_nucleation", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceNucleation;

  TD::run_ice_nucleation_phys();
  TD::run_ice_nucleation_bfb();
}

}
