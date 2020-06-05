#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

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
struct UnitWrap::UnitTest<D>::TestDropletActivation {

  static void run_droplet_activation_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    constexpr bool   log_predictNc = false;
    constexpr Scalar odt = 9.558E-04;

    DropletActivationData self[max_pack_size] = {
      // temp,    pres,      qv,        qc,        inv_rho,   sup,       xxlv,      npccn,     odt, log,           qcnuc,     ncnuc
      {2.106E+02, 1.020E+01, 1.002E-02, 2.034E-02, 8.852E-01, 0.974E-03, 9.221E+03, 5.100E-01, odt, log_predictNc, 9.012E+01, 1.023E+02},
      {2.052E+02, 2.010E+01, 2.982E-02, 2.134E-02, 8.852E-01, 0.874E-03, 8.221E+03, 4.100E-01, odt, log_predictNc, 8.014E+01, 2.310E+02},
      {2.352E+02, 3.010E+01, 3.982E-02, 3.420E-02, 8.900E-01, 0.723E-03, 7.221E+03, 3.100E-01, odt, log_predictNc, 7.345E+01, 3.210E+02},
      {2.212E+02, 4.023E+01, 4.982E-02, 4.423E-02, 9.900E-01, 0.623E-03, 6.221E+03, 2.100E-01, odt, log_predictNc, 6.421E+01, 4.231E+02},

      {2.251E+02, 5.901E+01, 5.981E-02, 5.592E-02, 0.100E+01, 0.574E-03, 5.221E+03, 1.100E-01, odt, log_predictNc, 5.423E+01, 5.671E+02},
      {2.710E+02, 6.310E+01, 6.902E-02, 6.213E-02, 0.100E+01, 0.474E-03, 4.221E+03, 8.100E-02, odt, log_predictNc, 4.012E+01, 6.045E+02},
      {2.052E+02, 7.120E+01, 7.902E-02, 7.123E-02, 0.100E+01, 0.323E-07, 3.221E+03, 4.100E-02, odt, log_predictNc, 3.214E+01, 7.231E+02},
      {2.502E+02, 8.900E+01, 8.901E-02, 8.014E-02, 0.100E+01, 0.223E-03, 2.221E+03, 2.100E-02, odt, log_predictNc, 2.190E+01, 8.923E+02},

      {2.552E+02, 9.230E+01, 9.120E-02, 9.234E-02, 0.950E+00, 0.150E-06, 9.221E+02, 9.952E-02, odt, log_predictNc, 1.320E+01, 9.821E+02},
      {2.452E+02, 1.220E+02, 9.320E-02, 1.902E-01, 0.950E+00, 1.974E-06, 8.221E+02, 4.952E-02, odt, log_predictNc, 9.024E+00, 1.092E+03},
      {2.352E+02, 1.320E+02, 1.023E-01, 2.983E-01, 0.950E+00, 0.823E-06, 7.221E+02, 1.952E-02, odt, log_predictNc, 8.723E+00, 1.231E+03},
      {2.252E+02, 1.456E+02, 1.243E-01, 3.234E-01, 0.950E+00, 0.723E-06, 6.221E+02, 9.952E-02, odt, log_predictNc, 7.324E+00, 1.346E+03},

      {1.990E+02, 1.623E+02, 1.334E-01, 4.231E-01, 1.069E+00, 0.674E-06, 5.221E+01, 6.952E-01, odt, log_predictNc, 6.832E+00, 1.532E+03},
      {2.952E+02, 1.670E+02, 1.445E-01, 5.782E-01, 1.069E+00, 1.574E-06, 4.221E+01, 3.952E-01, odt, log_predictNc, 5.346E+00, 1.753E+03},
      {2.852E+02, 1.980E+02, 1.650E-01, 6.743E-01, 1.069E+00, 0.423E-06, 3.221E+01, 1.952E-01, odt, log_predictNc, 4.312E+00, 1.982E+03},
      {2.702E+02, 2.091E+02, 1.982E-01, 9.621E-01, 1.069E+00, 0.323E-06, 2.221E+01, 9.952E-01, odt, log_predictNc, 3.245E+00, 2.130E+03}
    };

    // Sync to device
    KTH::view_1d<DropletActivationData> self_host("self_host", max_pack_size);
    view_1d<DropletActivationData> self_device("self_host", max_pack_size);
    std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      droplet_activation(self[i]);
     }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack temp, pres, qv, qc, inv_rho, sup, xxlv, npccn, qcnuc, ncnuc;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        temp[s]          = self_device(vs).temp;
        pres[s]          = self_device(vs).pres;
        qv[s]            = self_device(vs).qv;
        qc[s]            = self_device(vs).qc;
        inv_rho[s]       = self_device(vs).inv_rho;
        sup[s]           = self_device(vs).sup;
        xxlv[s]          = self_device(vs).xxlv;
        npccn[s]         = self_device(vs).npccn;
        qcnuc[s]         = self_device(vs).qcnuc;
        ncnuc[s]         = self_device(vs).ncnuc;
      }
      Functions::droplet_activation(temp, pres, qv, qc, inv_rho, sup, xxlv, npccn, log_predictNc, self_device(0).odt,
                                    qcnuc, ncnuc);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        self_device(vs).qcnuc = qcnuc[s];
        self_device(vs).ncnuc = ncnuc[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(self[s].qcnuc == self_host(s).qcnuc);
      REQUIRE(self[s].ncnuc == self_host(s).ncnuc);
    }
  }

  static void run_droplet_activation_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_droplet_activation", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDropletActivation;

  TD::run_droplet_activation_phys();
  TD::run_droplet_activation_bfb();
}

}
