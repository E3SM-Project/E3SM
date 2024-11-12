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
struct UnitWrap::UnitTest<D>::TestIceRelaxationTimescale : public UnitWrap::UnitTest<D>::Base {

  void run_ice_relaxation_timescale_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    IceRelaxationData self[max_pack_size] = {

     // rho,     temp,    rhofaci,   table_val_qi2qr_melting,    table_val_qi2qr_vent_melt,     dv,       mu,          sc,      qi_incld, ni_incld, epsi,epsi_tot
      {4.056E-03, 1.021E+01, 8.852E-01, 0.174E+00, 0.021E+00, 1.221E-14, 5.100E-03, 9.558E-04, 1.234E-03, 9.952E+03},
      {6.852E-02, 2.022E+01, 8.852E-01, 0.374E+00, 0.042E+00, 1.221E-13, 4.100E-03, 9.558E-04, 2.670E-03, 9.952E+03},
      {8.852E-02, 3.086E+01, 8.900E-01, 0.123E+00, 0.081E+00, 1.221E-12, 3.100E-03, 9.558E-04, 3.451E-03, 9.952E+03},
      {1.902E-01, 5.078E+01, 9.900E-01, 0.123E+00, 0.101E+00, 1.221E-11, 2.100E-03, 9.558E-04, 4.135E-03, 9.952E+03},

      {2.201E-01, 1.300E+02, 0.100E+01, 0.174E+00, 0.112E+00, 1.221E-10, 1.100E-03, 2.558E-05, 5.672E-03, 9.952E+04},
      {3.502E-01, 2.409E+02, 0.100E+01, 0.374E+00, 0.140E+00, 1.221E-09, 8.100E-04, 2.558E-05, 6.432E-03, 9.952E+04},
      {4.852E-01, 3.490E+02, 0.100E+01, 0.123E+00, 0.210E+00, 1.221E-08, 4.100E-04, 2.558E-05, 7.412E-03, 9.952E+04},
      {5.852E-01, 4.690E+02, 0.100E+01, 0.123E+00, 0.321E+00, 1.221E-07, 2.100E-04, 2.558E-05, 8.021E-03, 9.952E+04},

      {6.852E-01, 5.021E+02, 0.950E+00, 0.150E+00, 0.432E+00, 1.221E-06, 9.952E-05, 4.596E-05, 9.834E-03, 1.734E+04},
      {7.852E-01, 6.213E+02, 0.950E+00, 0.374E+00, 0.543E+00, 1.221E-05, 4.952E-05, 4.596E-05, 1.213E-02, 1.734E+04},
      {8.852E-01, 7.012E+02, 0.950E+00, 0.123E+00, 0.671E+00, 1.221E-04, 1.952E-05, 4.596E-05, 1.346E-02, 1.734E+04},
      {9.852E-01, 8.123E+02, 0.950E+00, 0.123E+00, 0.982E+00, 1.221E-03, 9.952E-06, 4.596E-05, 3.589E-02, 1.734E+04},

      {1.002E+01, 9.321E+02, 1.069E+00, 0.174E+00, 1.201E+00, 1.221E-02, 6.952E-06, 6.596E-05, 6.982E-02, 1.734E+04},
      {1.152E+01, 1.023E+03, 1.069E+00, 0.374E+00, 1.678E+00, 1.221E-02, 3.952E-06, 6.596E-05, 9.234E-02, 1.734E+04},
      {1.252E+01, 2.012E+03, 1.069E+00, 0.123E+00, 2.312E+00, 1.221E-02, 1.952E-06, 6.596E-05, 2.345E-01, 1.734E+04},
      {1.352E+01, 3.210E+03, 1.069E+00, 0.123E+00, 3.456E+00, 1.221E-02, 9.952E-07, 6.596E-05, 4.532E-01, 1.734E+04}
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        self[i].read(Base::m_fid);
      }
    }

    // Sync to device
    KTH::view_1d<IceRelaxationData> self_host("self_host", max_pack_size);
    view_1d<IceRelaxationData> self_device("self_host", max_pack_size);
    std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, mu, sc, qi_incld, ni_incld;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rho[s]                        = self_device(vs).rho;
        temp[s]                       = self_device(vs).temp;
        rhofaci[s]                    = self_device(vs).rhofaci;
        table_val_qi2qr_melting[s]    = self_device(vs).table_val_qi2qr_melting;
        table_val_qi2qr_vent_melt[s]  = self_device(vs).table_val_qi2qr_vent_melt;
        dv[s]                         = self_device(vs).dv;
        mu[s]                         = self_device(vs).mu;
        sc[s]                         = self_device(vs).sc;
        qi_incld[s]                   = self_device(vs).qi_incld;
        ni_incld[s]                   = self_device(vs).ni_incld;
      }

      Spack epsi{0.0};
      Spack epsi_tot{0.0};
      Functions::ice_relaxation_timescale(rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, mu, sc, qi_incld, ni_incld,
                                          epsi, epsi_tot);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        self_device(vs).epsi     = epsi[s];
        self_device(vs).epsi_tot = epsi_tot[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(self[s].epsi     == self_host(s).epsi);
        REQUIRE(self[s].epsi_tot == self_host(s).epsi_tot);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        self_host(s).write(Base::m_fid);
      }
    }
  }

  void run_ice_relaxation_timescale_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_ice_relaxation_timescale", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceRelaxationTimescale;

  T t;
  t.run_ice_relaxation_timescale_phys();
  t.run_ice_relaxation_timescale_bfb();
}

}
