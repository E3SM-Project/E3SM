#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
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
struct UnitWrap::UnitTest<D>::TestIceCldliqWetGrowth : public UnitWrap::UnitTest<D>::Base {

  void run_ice_cldliq_wet_growth_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    constexpr Scalar latvap = C::LatVap;
    constexpr Scalar latice = C::LatIce;

    IceWetGrowthData self[max_pack_size] = {
      // rho,temp,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,kap,mu,sc,qv,qc_incld,qi_incld,ni_incld,qr_incld,log_wetgrowth,qr2qi_collect_tend,qc2qi_collect_tend,qc_growth_rate,nr_ice_shed_tend,qc2qr_ice_shed_tend
      {4.056000E-03, 1.023000E+02, 1.201000E+02, 9.002000E-04, 8.215000E-04, 8.852000E-01, latvap, latice, 5.100000E-03, 9.558000E-04, 1.213000E-03, 9.653000E-04, 1.023000E-01, 4.098000E-02, 2.098000E-02, 9.952000E+03, 1.023000E-05, false, 1.241456E-04, 9.021345E-02, 1.043000E-01, 1.921000E-02, 0.242000E-02},
      {6.852000E-02, 1.120000E+02, 2.450000E+02, 9.321000E-04, 9.124000E-04, 8.852000E-01, latvap, latice, 4.100000E-03, 9.558000E-04, 2.560000E-03, 1.764000E-03, 2.346000E-01, 5.632000E-02, 3.024000E-02, 9.952000E+03, 2.093000E-05, false, 2.341678E-04, 1.092432E-02, 2.903000E-01, 2.125000E-02, 0.342000E-02},
      {8.852000E-02, 1.210000E+02, 3.420000E+02, 9.623000E-04, 9.432000E-04, 8.900000E-01, latvap, latice, 3.100000E-03, 9.558000E-04, 3.211000E-03, 3.421000E-03, 3.421000E-01, 6.542000E-02, 4.567000E-02, 9.952000E+03, 3.091000E-05, false, 3.215234E-04, 2.098987E-02, 3.450000E-01, 3.490000E-02, 0.932000E-02},
      {1.902000E-01, 1.326000E+02, 4.321000E+02, 9.982000E-04, 9.623000E-04, 9.900000E-01, latvap, latice, 2.100000E-03, 9.558000E-04, 4.121000E-03, 4.569000E-03, 4.673000E-01, 7.902000E-02, 5.321000E-02, 9.952000E+03, 4.521000E-05, false, 4.675567E-04, 3.214982E-02, 4.290000E-01, 4.590000E-02, 1.025000E-01},

      {2.201000E-01, 1.456000E+02, 5.670000E+02, 1.234000E-03, 9.723000E-04, 0.100000E+01, latvap, latice, 1.100000E-03, 2.550008E-05, 4.980000E-03, 5.621000E-03, 5.420000E-01, 8.021000E-02, 6.902000E-02, 9.952000E+04, 5.678000E-05, false, 5.389236E-04, 4.125969E-02, 5.098000E-01, 5.921000E-02, 2.031000E-01},
      {3.502000E-01, 1.780009E+02, 6.832000E+02, 1.562000E-03, 1.024000E-03, 0.100000E+01, latvap, latice, 8.100000E-04, 2.558000E-05, 5.643000E-03, 7.367000E-03, 6.782000E-01, 9.253000E-02, 8.045000E-02, 9.952000E+04, 6.902000E-05, false, 6.432654E-04, 5.389457E-02, 6.723000E-01, 6.093000E-02, 4.098000E-01},
      {4.852000E-01, 2.100009E+02, 7.090000E+02, 2.101000E-03, 1.235000E-03, 0.100000E+01, latvap, latice, 4.100000E-04, 2.558000E-05, 7.892000E-03, 9.087000E-03, 8.213000E-01, 1.256000E-01, 9.134000E-02, 9.952000E+04, 8.367000E-05, false, 7.210983E-04, 6.476985E-02, 8.902000E-01, 8.345000E-02, 8.023000E-01},
      {5.852000E-01, 2.310000E+02, 9.215000E+02, 2.312000E-03, 1.456000E-03, 0.100000E+01, latvap, latice, 2.100000E-04, 2.558000E-05, 9.321000E-03, 1.245000E-02, 1.067000E-00, 2.347000E-01, 1.092000E-01, 9.952000E+04, 9.098000E-05, false, 8.543367E-04, 8.213186E-02, 9.021000E-01, 9.321000E-02, 9.098000E-01},

      {6.852000E-01, 2.563000E+02, 1.089000E+03, 3.601000E-03, 1.864000E-03, 0.950000E+00, latvap, latice, 9.952000E-05, 4.596000E-05, 1.453000E-02, 2.543000E-02, 2.345000E-00, 3.578000E-01, 2.873000E-01, 1.734000E+04, 1.023000E-04, false, 9.021215E-04, 9.023367E-02, 1.023000E-00, 1.056000E-01, 1.256000E-00},
      {7.852000E-01, 2.789000E+02, 3.754000E+03, 3.891000E-03, 2.093000E-03, 0.950000E+00, latvap, latice, 4.952000E-05, 4.596000E-05, 2.789000E-02, 4.367000E-02, 3.890000E-00, 4.980000E-01, 3.468000E-01, 1.734000E+04, 2.146000E-04, false, 1.043468E-05, 1.094854E-02, 2.012000E-00, 2.893000E-01, 2.903000E-00},
      {8.852000E-01, 3.123000E+02, 8.902000E+03, 4.872000E-03, 2.345000E-03, 0.950000E+00, latvap, latice, 1.952000E-05, 4.596000E-05, 4.256000E-02, 6.324000E-02, 4.120000E-00, 6.321000E-01, 4.890000E-01, 1.734000E+04, 4.321000E-04, false, 2.341763E-05, 2.126247E-03, 3.120000E-00, 3.456000E-01, 3.912000E-00},
      {9.852000E-01, 4.981000E+02, 1.092000E+04, 5.210000E-03, 3.210000E-03, 0.950000E+00, latvap, latice, 9.952000E-06, 4.596000E-05, 6.821000E-02, 8.789000E-02, 5.320000E-00, 7.982000E-01, 6.921000E-01, 1.734000E+04, 5.821000E-04, false, 3.901479E-05, 3.874763E-03, 5.902000E-00, 5.092000E-01, 4.821000E-00},

      {1.002000E+01, 1.234000E+03, 2.125000E+04, 6.012000E-03, 5.902000E-03, 1.069000E+00, latvap, latice, 6.952000E-06, 6.596000E-05, 8.472000E-02, 1.543000E-01, 6.012000E-00, 8.902000E-01, 9.210000E-01, 1.734000E+04, 6.921000E-04, false, 4.521923E-05, 4.592698E-03, 6.091000E-00, 6.743000E-01, 5.602000E-00},
      {1.152000E+01, 2.120000E+03, 4.568000E+04, 6.342000E-03, 9.210000E-03, 1.069000E+00, latvap, latice, 3.952000E-06, 6.596000E-05, 1.098000E-01, 3.456000E-01, 7.241000E-00, 9.102000E-01, 1.002000E-00, 1.734000E+04, 7.901000E-04, false, 5.236542E-05, 5.678873E-03, 7.231000E-00, 8.321000E-01, 6.092000E-00},
      {1.252000E+01, 3.145000E+03, 8.213000E+04, 9.290000E-03, 1.034000E-02, 1.069000E+00, latvap, latice, 1.952000E-06, 6.596000E-05, 2.340006E-01, 5.632000E-01, 8.452000E-00, 1.003000E-01, 2.145000E-00, 1.734000E+04, 9.212000E-04, false, 6.732276E-05, 7.321873E-03, 8.234000E-00, 9.023000E-01, 7.201000E-00},
      {1.352000E+01, 4.742000E+03, 1.014000E+05, 1.234000E-02, 1.456000E-02, 1.069000E+00, latvap, latice, 9.952000E-07, 6.596000E-05, 4.123000E-01, 6.128000E-01, 9.076000E-00, 2.831000E-01, 3.902000E-00, 1.734000E+04, 1.023000E-03, false, 7.902887E-05, 9.032908E-03, 9.021000E-00, 1.092000E-01, 8.096000E-00}
    };

    // Sync to device
    KTH::view_1d<IceWetGrowthData> self_host("self_host", max_pack_size);
    view_1d<IceWetGrowthData> self_device("self_device", max_pack_size);
    std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        self[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack rho,temp, pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,dv,kap,mu,sc,
        qv,qc_incld,qi_incld,ni_incld,qr_incld;

      Smask log_wetgrowth;

      Spack qr2qi_collect_tend,qc2qi_collect_tend,qc_growth_rate,nr_ice_shed_tend,qc2qr_ice_shed_tend;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rho[s]                        = self_device(vs).rho;
        temp[s]                       = self_device(vs).temp;
        pres[s]                       = self_device(vs).pres;
        rhofaci[s]                    = self_device(vs).rhofaci;
        table_val_qi2qr_melting[s]    = self_device(vs).table_val_qi2qr_melting;
        table_val_qi2qr_vent_melt[s]  = self_device(vs).table_val_qi2qr_vent_melt;
        dv[s]                         = self_device(vs).dv;
        kap[s]                        = self_device(vs).kap;
        mu[s]                         = self_device(vs).mu;
        sc[s]                         = self_device(vs).sc;
        qv[s]                         = self_device(vs).qv;
        qc_incld[s]                   = self_device(vs).qc_incld;
        qi_incld[s]                   = self_device(vs).qi_incld;
        ni_incld[s]                   = self_device(vs).ni_incld;
        qr_incld[s]                   = self_device(vs).qr_incld;
        qr2qi_collect_tend[s]         = self_device(vs).qr2qi_collect_tend;
        qc2qi_collect_tend[s]         = self_device(vs).qc2qi_collect_tend;
        qc_growth_rate[s]             = self_device(vs).qc_growth_rate;
        nr_ice_shed_tend[s]           = self_device(vs).nr_ice_shed_tend;
        qc2qr_ice_shed_tend[s]        = self_device(vs).qc2qr_ice_shed_tend;
        log_wetgrowth.set(s, self_device(vs).log_wetgrowth);
      }

      Functions::ice_cldliq_wet_growth(rho, temp, pres, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, kap, mu, sc,
                                       qv, qc_incld, qi_incld, ni_incld, qr_incld,
                                       log_wetgrowth, qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        self_device(vs).log_wetgrowth       = log_wetgrowth[s];
        self_device(vs).qr2qi_collect_tend  = qr2qi_collect_tend[s];
        self_device(vs).qc2qi_collect_tend  = qc2qi_collect_tend[s];
        self_device(vs).qc_growth_rate      = qc_growth_rate[s];
        self_device(vs).nr_ice_shed_tend    = nr_ice_shed_tend[s];
        self_device(vs).qc2qr_ice_shed_tend = qc2qr_ice_shed_tend[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(static_cast<bool>(self[s].log_wetgrowth) == static_cast<bool>(self_host(s).log_wetgrowth));

        REQUIRE(self[s].qr2qi_collect_tend    == self_host(s).qr2qi_collect_tend);
        REQUIRE(self[s].qc2qi_collect_tend    == self_host(s).qc2qi_collect_tend);
        REQUIRE(self[s].qc_growth_rate        == self_host(s).qc_growth_rate);
        REQUIRE(self[s].nr_ice_shed_tend      == self_host(s).nr_ice_shed_tend);
        REQUIRE(self[s].qc2qr_ice_shed_tend   == self_host(s).qc2qr_ice_shed_tend);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        self_host(s).write(Base::m_fid);
      }
    }
  }

  void run_ice_cldliq_wet_growth_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_ice_cldliq_wet_growth", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceCldliqWetGrowth;

  T t;
  t.run_ice_cldliq_wet_growth_phys();
  t.run_ice_cldliq_wet_growth_bfb();
}

}
