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

template <typename D>
struct UnitWrap::UnitTest<D>::TestIceNucleation {

  static void run_ice_nucleation_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    constexpr bool   log_predictNc = false;
    constexpr Scalar odt           = 9.558E-04;

    IceNucleationData self[max_pack_size] = {
      // temp,    inv_rho,   nitot,     naai,      supi,      odt, log_predictNc
      {2.106E+02, 8.852E-01, 0.974E+04, 9.221E+03, 5.100E-01, odt, log_predictNc },
      {2.052E+02, 8.852E-01, 0.874E+04, 8.221E+03, 4.100E-01, odt, log_predictNc },
      {2.352E+02, 8.900E-01, 0.723E+04, 7.221E+03, 3.100E-01, odt, log_predictNc },
      {2.212E+02, 9.900E-01, 0.623E+04, 6.221E+03, 2.100E-01, odt, log_predictNc },

      {2.251E+02, 0.100E+01, 0.574E+04, 5.221E+03, 1.100E-01, odt, log_predictNc},
      {2.710E+02, 0.100E+01, 0.474E+04, 4.221E+03, 8.100E-02, odt, log_predictNc},
      {2.052E+02, 0.100E+01, 0.323E+04, 3.221E+03, 4.100E-02, odt, log_predictNc},
      {2.502E+02, 0.100E+01, 0.223E+04, 2.221E+03, 2.100E-02, odt, log_predictNc},

      {2.552E+02, 0.950E+00, 0.150E+04, 9.221E+02, 9.952E-02, odt, log_predictNc},
      {2.452E+02, 0.950E+00, 0.974E+03, 8.221E+02, 4.952E-02, odt, log_predictNc},
      {2.352E+02, 0.950E+00, 0.823E+03, 7.221E+02, 1.952E-02, odt, log_predictNc},
      {2.252E+02, 0.950E+00, 0.723E+03, 6.221E+02, 9.952E-02, odt, log_predictNc},

      {1.990E+02, 1.069E+00, 0.674E+03, 5.221E+01, 6.952E-01, odt, log_predictNc },
      {2.952E+02, 1.069E+00, 0.574E+03, 4.221E+01, 3.952E-01, odt, log_predictNc },
      {2.852E+02, 1.069E+00, 0.423E+03, 3.221E+01, 1.952E-01, odt, log_predictNc },
      {2.702E+02, 1.069E+00, 0.323E+03, 2.221E+01, 9.952E-01, odt, log_predictNc }
    };

    // Get data from fortran
    for (Int i = 0; i < Spack::n; ++i) {
      ice_nucleation(self[i]);
    }

    // Sync to device
    KTH::view_1d<IceNucleationData> self_host("self_host", Spack::n);
    view_1d<IceNucleationData> self_device("self_host", Spack::n);
    std::copy(&self[0], &self[0] + Spack::n, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
     // Init pack inputs
      Spack temp, inv_rho, nitot, naai, supi;
      for (Int s = 0; s < Spack::n; ++s) {
        temp[s]          = self_device(s).temp;
        inv_rho[s]       = self_device(s).inv_rho;
        nitot[s]         = self_device(s).nitot;
        naai[s]          = self_device(s).naai;
        supi[s]          = self_device(s).supi;
      }
      // outputs
      Spack qinuc{0.0};
      Spack ninuc{0.0};
      Functions::ice_nucleation(temp, inv_rho, nitot, naai, supi, self_device(0).odt, log_predictNc, qinuc, ninuc);

      for (Int s = 0; s < Spack::n; ++s) {
        self_device(s).qinuc = qinuc[s];
        self_device(s).ninuc = ninuc[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(self[s].qinuc == self_host(s).qinuc);
      REQUIRE(self[s].ninuc == self_host(s).ninuc);
    }
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
