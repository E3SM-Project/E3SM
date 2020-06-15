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
struct UnitWrap::UnitTest<D>::TestIncloudMixing {

  static void run_incloud_mixing_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    constexpr Scalar qsmall = C::QSMALL;
    constexpr Scalar qc0 = 0.1*qsmall, qc1 = 0.9*qsmall, qc2 = 2.*qsmall, qc3 = 100.*qsmall;
    constexpr Scalar qr0 = 0.3*qsmall, qr1 = 0.6*qsmall, qr2 = 5.*qsmall, qr3 = 200.*qsmall;
    constexpr Scalar qitot0 = 0.1*qsmall, qitot1 = 0.5*qsmall, qitot2 = 12.*qsmall, qitot3 = 300.*qsmall;
    constexpr Scalar qirim0 = 0.1*qsmall, qirim1 = 0.5*qsmall, qirim2 = 12.*qsmall, qirim3 = 300.*qsmall;

    constexpr Scalar nc0 = 1.052E+05, nc1 = 4.952E+06, nc2 = 1.340E+07, nc3 = 9.652E+08;
    constexpr Scalar nr0 = 2.052E+05, nr1 = 5.952E+06, nr2 = 2.340E+07, nr3 = 8.652E+08;
    constexpr Scalar nitot0 = 3.052E+05, nitot1 = 6.952E+06, nitot2 = 3.340E+07, nitot3 = 7.652E+08;
    constexpr Scalar birim0 = 4.052E+05, birim1 = 7.952E+06, birim2 = 4.340E+07, birim3 = 6.652E+08;

    constexpr Scalar inv_lcldm0 = 1.052E+01, inv_lcldm1 = 2.952E+02, inv_lcldm2 = 1.340E+03, inv_lcldm3 = 1.652E+04;
    constexpr Scalar inv_icldm0 = 2.052E+01, inv_icldm1 = 4.952E+02, inv_icldm2 = 3.340E+03, inv_icldm3 = 2.652E+04;
    constexpr Scalar inv_rcldm0 = 3.052E+01, inv_rcldm1 = 5.952E+02, inv_rcldm2 = 3.340E+03, inv_rcldm3 = 3.652E+04;

    IncloudMixingData self[max_pack_size] = {
      //qc, qr, qitot, qirim, nc, nr, nitot, birim, inv_lcldm, inv_icldm, inv_rcldm
      {qc0, qr3, qitot0, qirim3, nc0, nr3, nitot0, birim3, inv_lcldm0, inv_icldm3, inv_rcldm0},
      {qc1, qr2, qitot1, qirim2, nc1, nr2, nitot1, birim2, inv_lcldm1, inv_icldm2, inv_rcldm1},
      {qc2, qr1, qitot2, qirim1, nc2, nr1, nitot2, birim1, inv_lcldm2, inv_icldm1, inv_rcldm2},
      {qc3, qr0, qitot3, qirim0, nc3, nr0, nitot3, birim0, inv_lcldm3, inv_icldm0, inv_rcldm3},

      {qc0, qr0, qitot0, qirim0, nc0, nr0, nitot0, birim0, inv_lcldm0, inv_icldm0, inv_rcldm0},
      {qc1, qr1, qitot1, qirim1, nc1, nr1, nitot1, birim1, inv_lcldm1, inv_icldm1, inv_rcldm1},
      {qc2, qr2, qitot2, qirim2, nc2, nr2, nitot2, birim2, inv_lcldm2, inv_icldm2, inv_rcldm2},
      {qc3, qr3, qitot3, qirim3, nc3, nr3, nitot3, birim3, inv_lcldm3, inv_icldm3, inv_rcldm3},

      {qc3, qr0, qitot3, qirim0, nc3, nr0, nitot3, birim0, inv_lcldm3, inv_icldm0, inv_rcldm3},
      {qc2, qr1, qitot2, qirim1, nc2, nr1, nitot2, birim1, inv_lcldm2, inv_icldm1, inv_rcldm2},
      {qc1, qr2, qitot1, qirim2, nc1, nr2, nitot1, birim2, inv_lcldm1, inv_icldm2, inv_rcldm1},
      {qc0, qr3, qitot0, qirim3, nc0, nr3, nitot0, birim3, inv_lcldm0, inv_icldm3, inv_rcldm0},

      {qc3, qr2, qitot1, qirim3, nc2, nr1, nitot2, birim1, inv_lcldm1, inv_icldm3, inv_rcldm1},
      {qc2, qr3, qitot2, qirim2, nc3, nr2, nitot1, birim2, inv_lcldm3, inv_icldm1, inv_rcldm2},
      {qc0, qr0, qitot3, qirim1, nc0, nr3, nitot3, birim3, inv_lcldm2, inv_icldm0, inv_rcldm3},
      {qc1, qr1, qitot0, qirim0, nc1, nr0, nitot0, birim0, inv_lcldm0, inv_icldm2, inv_rcldm0}
    };

    // Sync to device
    KTH::view_1d<IncloudMixingData> self_host("self_host", max_pack_size);
    view_1d<IncloudMixingData> self_device("self_host", max_pack_size);
    std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
       calculate_incloud_mixingratios(self[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

     // Init pack inputs
      Spack qc, qr, qitot, qirim, nc, nr, nitot, birim, inv_lcldm, inv_icldm, inv_rcldm;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qc[s]            = self_device(vs).qc;
        qr[s]            = self_device(vs).qr;
        qitot[s]         = self_device(vs).qitot;
        qirim[s]         = self_device(vs).qirim;
        nc[s]            = self_device(vs).nc;
        nr[s]            = self_device(vs).nr;
        nitot[s]         = self_device(vs).nitot;
        birim[s]         = self_device(vs).birim;
        inv_lcldm[s]     = self_device(vs).inv_lcldm;
        inv_icldm[s]     = self_device(vs).inv_icldm;
        inv_rcldm[s]     = self_device(vs).inv_rcldm;
      }
      // outputs
      Spack qc_incld{0.}, qr_incld{0.}, qitot_incld{0.}, qirim_incld{0.};
      Spack nc_incld{0.}, nr_incld{0.}, nitot_incld{0.}, birim_incld{0.};

      Functions::calculate_incloud_mixingratios(qc, qr, qitot, qirim, nc, nr, nitot, birim, inv_lcldm, inv_icldm, inv_rcldm,
                                                qc_incld, qr_incld, qitot_incld, qirim_incld,
                                                nc_incld, nr_incld, nitot_incld, birim_incld);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        self_device(vs).qc_incld    = qc_incld[s];
        self_device(vs).qr_incld    = qr_incld[s];
        self_device(vs).qitot_incld = qitot_incld[s];
        self_device(vs).qirim_incld = qirim_incld[s];
        self_device(vs).nc_incld    = nc_incld[s];
        self_device(vs).nr_incld    = nr_incld[s];
        self_device(vs).nitot_incld = nitot_incld[s];
        self_device(vs).birim_incld = birim_incld[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(self[s].qc_incld    == self_host(s).qc_incld);
      REQUIRE(self[s].qr_incld    == self_host(s).qr_incld);
      REQUIRE(self[s].qitot_incld == self_host(s).qitot_incld);
      REQUIRE(self[s].qirim_incld == self_host(s).qirim_incld);
      REQUIRE(self[s].nc_incld    == self_host(s).nc_incld);
      REQUIRE(self[s].nr_incld    == self_host(s).nr_incld);
      REQUIRE(self[s].nitot_incld == self_host(s).nitot_incld);
      REQUIRE(self[s].birim_incld == self_host(s).birim_incld);
    }
  }

  static void run_incloud_mixing_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_incloud_mixingratios", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIncloudMixing;

  TD::run_incloud_mixing_phys();
  TD::run_incloud_mixing_bfb();
}

}
