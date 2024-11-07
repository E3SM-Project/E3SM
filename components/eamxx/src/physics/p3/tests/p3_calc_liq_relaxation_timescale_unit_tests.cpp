#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

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
struct UnitWrap::UnitTest<D>::TestCalcLiqRelaxationTimescale : public UnitWrap::UnitTest<D>::Base {

  void run_phys()
  {
    // TODO
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Read in tables
    view_2d_table vn_table_vals, vm_table_vals, revap_table_vals;
    view_1d_table mu_r_table_vals;
    view_dnu_table dnu;
    Functions::init_kokkos_tables(vn_table_vals, vm_table_vals, revap_table_vals, mu_r_table_vals, dnu);

    using KTH = KokkosTypes<HostDevice>;

    // Set up input data.
    constexpr Scalar qsmall = C::QSMALL;
    constexpr Scalar qr_small = 0.9 * qsmall;
    constexpr Scalar qr_not_small = 2.0 * qsmall;
    constexpr Scalar qc_small = 0.9 * qsmall;
    constexpr Scalar qc_not_small = 2.0 * qsmall;
    CalcLiqRelaxationData self[max_pack_size];
    for (Int i = 0; i < max_pack_size; ++i) {
      self[i].randomize(engine);
      self[i].qr_incld = (i % 2) ? qr_small : qr_not_small;
      self[i].qc_incld = ((i/2) % 2) ? qc_small : qc_not_small;
      self[i].f1r = C::f1r;
      self[i].f2r = C::f2r;
    }

  // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        self[i].read(Base::m_fid);
      }
    }

    // Sync to device
    KTH::view_1d<CalcLiqRelaxationData> self_host("self_host", max_pack_size);
    view_1d<CalcLiqRelaxationData> self_device("self_host", max_pack_size);
    std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack rho, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rho[s]      = self_device(vs).rho;
        dv[s]       = self_device(vs).dv;
        mu[s]       = self_device(vs).mu;
        sc[s]       = self_device(vs).sc;
        mu_r[s]     = self_device(vs).mu_r;
        lamr[s]     = self_device(vs).lamr;
        cdistr[s]   = self_device(vs).cdistr;
        cdist[s]    = self_device(vs).cdist;
        qr_incld[s] = self_device(vs).qr_incld;
        qc_incld[s] = self_device(vs).qc_incld;
      }

      Spack epsr{0.0}, epsc{0.0};
      Functions::calc_liq_relaxation_timescale(revap_table_vals, rho, self_device(0).f1r, self_device(0).f2r, dv,
        mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld, epsr, epsc);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        self_device(vs).epsr = epsr[s];
        self_device(vs).epsc = epsc[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(self[s].epsr == self_host(s).epsr);
        REQUIRE(self[s].epsc == self_host(s).epsc);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        self_host(s).write(Base::m_fid);
      }
    }
  }

};

}
}
}

namespace {

TEST_CASE("p3_calc_liq_relaxation_timescale", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcLiqRelaxationTimescale;

  T t;
  t.run_phys();
  t.run_bfb();
}

}
