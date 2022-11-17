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
struct UnitWrap::UnitTest<D>::TestIncloudMixing {

  static void run_incloud_mixing_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    constexpr Scalar qsmall = C::QSMALL;
    constexpr Scalar qc0 = 0.1*qsmall, qc1 = 0.9*qsmall, qc2 = 2.*qsmall, qc3 = 100.*qsmall;
    constexpr Scalar qr0 = 0.3*qsmall, qr1 = 0.6*qsmall, qr2 = 5.*qsmall, qr3 = 200.*qsmall;
    constexpr Scalar qi0 = 0.1*qsmall, qi1 = 0.5*qsmall, qi2 = 12.*qsmall, qi3 = 300.*qsmall;
    constexpr Scalar qm0 = 0.1*qsmall, qm1 = 0.5*qsmall, qm2 = 12.*qsmall, qm3 = 300.*qsmall;

    constexpr Scalar nc0 = 1.052E+05, nc1 = 4.952E+06, nc2 = 1.340E+07, nc3 = 9.652E+08;
    constexpr Scalar nr0 = 2.052E+05, nr1 = 5.952E+06, nr2 = 2.340E+07, nr3 = 8.652E+08;
    constexpr Scalar ni0 = 3.052E+05, ni1 = 6.952E+06, ni2 = 3.340E+07, ni3 = 7.652E+08;
    constexpr Scalar bm0 = 4.052E+05, bm1 = 7.952E+06, bm2 = 4.340E+07, bm3 = 6.652E+08;

    constexpr Scalar inv_cld_frac_l0 = 1.052E+01, inv_cld_frac_l1 = 2.952E+02, inv_cld_frac_l2 = 1.340E+03, inv_cld_frac_l3 = 1.652E+04;
    constexpr Scalar inv_cld_frac_i0 = 2.052E+01, inv_cld_frac_i1 = 4.952E+02, inv_cld_frac_i2 = 3.340E+03, inv_cld_frac_i3 = 2.652E+04;
    constexpr Scalar inv_cld_frac_r0 = 3.052E+01, inv_cld_frac_r1 = 5.952E+02, inv_cld_frac_r2 = 3.340E+03, inv_cld_frac_r3 = 3.652E+04;

    IncloudMixingData self[max_pack_size] = {
      //qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r
      {qc0, qr3, qi0, qm3, nc0, nr3, ni0, bm3, inv_cld_frac_l0, inv_cld_frac_i3, inv_cld_frac_r0},
      {qc1, qr2, qi1, qm2, nc1, nr2, ni1, bm2, inv_cld_frac_l1, inv_cld_frac_i2, inv_cld_frac_r1},
      {qc2, qr1, qi2, qm1, nc2, nr1, ni2, bm1, inv_cld_frac_l2, inv_cld_frac_i1, inv_cld_frac_r2},
      {qc3, qr0, qi3, qm0, nc3, nr0, ni3, bm0, inv_cld_frac_l3, inv_cld_frac_i0, inv_cld_frac_r3},

      {qc0, qr0, qi0, qm0, nc0, nr0, ni0, bm0, inv_cld_frac_l0, inv_cld_frac_i0, inv_cld_frac_r0},
      {qc1, qr1, qi1, qm1, nc1, nr1, ni1, bm1, inv_cld_frac_l1, inv_cld_frac_i1, inv_cld_frac_r1},
      {qc2, qr2, qi2, qm2, nc2, nr2, ni2, bm2, inv_cld_frac_l2, inv_cld_frac_i2, inv_cld_frac_r2},
      {qc3, qr3, qi3, qm3, nc3, nr3, ni3, bm3, inv_cld_frac_l3, inv_cld_frac_i3, inv_cld_frac_r3},

      {qc3, qr0, qi3, qm0, nc3, nr0, ni3, bm0, inv_cld_frac_l3, inv_cld_frac_i0, inv_cld_frac_r3},
      {qc2, qr1, qi2, qm1, nc2, nr1, ni2, bm1, inv_cld_frac_l2, inv_cld_frac_i1, inv_cld_frac_r2},
      {qc1, qr2, qi1, qm2, nc1, nr2, ni1, bm2, inv_cld_frac_l1, inv_cld_frac_i2, inv_cld_frac_r1},
      {qc0, qr3, qi0, qm3, nc0, nr3, ni0, bm3, inv_cld_frac_l0, inv_cld_frac_i3, inv_cld_frac_r0},

      {qc3, qr2, qi1, qm3, nc2, nr1, ni2, bm1, inv_cld_frac_l1, inv_cld_frac_i3, inv_cld_frac_r1},
      {qc2, qr3, qi2, qm2, nc3, nr2, ni1, bm2, inv_cld_frac_l3, inv_cld_frac_i1, inv_cld_frac_r2},
      {qc0, qr0, qi3, qm1, nc0, nr3, ni3, bm3, inv_cld_frac_l2, inv_cld_frac_i0, inv_cld_frac_r3},
      {qc1, qr1, qi0, qm0, nc1, nr0, ni0, bm0, inv_cld_frac_l0, inv_cld_frac_i2, inv_cld_frac_r0}
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
      Spack qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qc[s]             = self_device(vs).qc;
        qr[s]             = self_device(vs).qr;
        qi[s]             = self_device(vs).qi;
        qm[s]             = self_device(vs).qm;
        nc[s]             = self_device(vs).nc;
        nr[s]             = self_device(vs).nr;
        ni[s]             = self_device(vs).ni;
        bm[s]             = self_device(vs).bm;
        inv_cld_frac_l[s] = self_device(vs).inv_cld_frac_l;
        inv_cld_frac_i[s] = self_device(vs).inv_cld_frac_i;
        inv_cld_frac_r[s] = self_device(vs).inv_cld_frac_r;
      }
      // outputs
      Spack qc_incld{0.}, qr_incld{0.}, qi_incld{0.}, qm_incld{0.};
      Spack nc_incld{0.}, nr_incld{0.}, ni_incld{0.}, bm_incld{0.};

      Functions::calculate_incloud_mixingratios(qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
                                                qc_incld, qr_incld, qi_incld, qm_incld,
                                                nc_incld, nr_incld, ni_incld, bm_incld);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        self_device(vs).qc_incld = qc_incld[s];
        self_device(vs).qr_incld = qr_incld[s];
        self_device(vs).qi_incld = qi_incld[s];
        self_device(vs).qm_incld = qm_incld[s];
        self_device(vs).nc_incld = nc_incld[s];
        self_device(vs).nr_incld = nr_incld[s];
        self_device(vs).ni_incld = ni_incld[s];
        self_device(vs).bm_incld = bm_incld[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    if (SCREAM_BFB_TESTING) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(self[s].qc_incld == self_host(s).qc_incld);
        REQUIRE(self[s].qr_incld == self_host(s).qr_incld);
        REQUIRE(self[s].qi_incld == self_host(s).qi_incld);
        REQUIRE(self[s].qm_incld == self_host(s).qm_incld);
        REQUIRE(self[s].nc_incld == self_host(s).nc_incld);
        REQUIRE(self[s].nr_incld == self_host(s).nr_incld);
        REQUIRE(self[s].ni_incld == self_host(s).ni_incld);
        REQUIRE(self[s].bm_incld == self_host(s).bm_incld);
      }
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
