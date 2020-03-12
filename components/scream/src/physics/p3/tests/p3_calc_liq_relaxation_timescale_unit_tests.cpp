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
struct UnitWrap::UnitTest<D>::TestCalcLiqRelaxationTimescale {

  static void run_phys()
  {
    // TODO
  }

  static void run_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    CalcLiqRelaxationData self[max_pack_size];
    for (Int i = 0; i < Spack::n; ++i) {
      self[i].randomize();
    }

    // Get data from fortran
    for (Int i = 0; i < Spack::n; ++i) {
      calc_liq_relaxation_timescale(self[i]);
    }

    // Sync to device
    KTH::view_1d<CalcLiqRelaxationData> self_host("self_host", Spack::n);
    view_1d<CalcLiqRelaxationData> self_device("self_host", Spack::n);
    std::copy(&self[0], &self[0] + Spack::n, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
    // Init pack inputs
    Spack rho, f1r, f2r, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld,
      qc_incld;

      for (Int s = 0; s < Spack::n; ++s) {
        rho[s]      = self_device(s).rho;
        f1r[s]      = self_device(s).f1r;
        f2r[s]      = self_device(s).f2r;
        dv[s]       = self_device(s).dv;
        mu[s]       = self_device(s).mu;
        sc[s]       = self_device(s).sc;
        mu_r[s]     = self_device(s).mu_r;
        lamr[s]     = self_device(s).lamr;
        cdistr[s]   = self_device(s).cdistr;
        cdist[s]    = self_device(s).cdist;
        qr_incld[s] = self_device(s).qr_incld;
        qc_incld[s] = self_device(s).qc_incld;
      }

      Spack epsr{0.0}, epsc{0.0};
      Functions::calc_liq_relaxation_timescale(rho, f1r, f2r, dv, mu, sc, mu_r,
        lamr, cdistr, cdist, qr_incld, qc_incld, epsr, epsc);

      for (Int s = 0; s < Spack::n; ++s) {
        self_device(s).epsr = epsr[s];
        self_device(s).epsc = epsc[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(self[s].epsr == self_host(s).epsr);
      REQUIRE(self[s].epsc == self_host(s).epsc);
    }
  }

};

}
}
}

namespace {

TEST_CASE("p3_calc_liq_relaxation_timescale", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcLiqRelaxationTimescale;

  TD::run_phys();
  TD::run_bfb();
}

}
