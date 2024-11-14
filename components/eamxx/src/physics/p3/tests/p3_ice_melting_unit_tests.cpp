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
struct UnitWrap::UnitTest<D>::TestP3IceMelting : public UnitWrap::UnitTest<D>::Base
{

void ice_melting_bfb() {
  constexpr Scalar latvap = C::LatVap;
  constexpr Scalar latice = C::LatIce;

  // make array of input data (why not pass actual variables?). Copied 1st 4 rows 4x to fill pack size.
  IceMeltingData IceMelt[max_pack_size] = {
    //rho,     T_atm,        pres,     rhofaci,  table_val_qi2qr_melting,   table_val_qi2qr_vent_melt,   latent_heat_vapor,     latent_heat_fusion,      dv,       sc,       mu,       kap,      qv,       qi_incld,ni_incld
    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,latvap,latice,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,latvap,latice,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,latvap,latice,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,latvap,latice,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05},

    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,latvap,latice,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,latvap,latice,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,latvap,latice,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,latvap,latice,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05},

    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,latvap,latice,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,latvap,latice,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,latvap,latice,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,latvap,latice,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05},

    {0.117E+01,0.299E+03,0.101E+06,0.829E+00,0.122E+01,0.562E-01,latvap,latice,0.263E-04,0.601E+00,0.185E-04,0.261E-01,0.160E-01,0.510E-02,  0.195E-12},
    {0.114E+01,0.296E+03,0.973E+05,0.842E+00,0.122E+01,0.562E-01,latvap,latice,0.268E-04,0.601E+00,0.183E-04,0.259E-01,0.149E-01,0.510E-02,  0.195E-12},
    {0.977E+00,0.287E+03,0.809E+05,0.913E+00,0.122E+01,0.562E-01,latvap,latice,0.306E-04,0.599E+00,0.179E-04,0.253E-01,0.827E-02,0.000E+00,  0.000E+00},
    {0.103E+01,0.289E+03,0.862E+05,0.887E+00,0.636E-03,0.281E-04,latvap,latice,0.291E-04,0.600E+00,0.180E-04,0.254E-01,0.107E-01,0.510E-02,  0.336E+05}
  };

  // Sync to device
  view_1d<IceMeltingData> IceMelt_device("IceMelt", max_pack_size);
  auto IceMelt_host = Kokkos::create_mirror_view(IceMelt_device);
  // This copy only copies the input variables.
  std::copy(&IceMelt[0], &IceMelt[0] + max_pack_size, IceMelt_host.data());
  Kokkos::deep_copy(IceMelt_device, IceMelt_host);

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < max_pack_size; ++i) {
      IceMelt[i].read(Base::m_fid);
    }
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,dv,sc,mu,kap,qv,qi_incld,ni_incld,qi2qr_melt_tend,ni2nr_melt_tend;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      rho[s]                       = IceMelt_device(vs).rho;
      T_atm[s]                     = IceMelt_device(vs).T_atm;
      pres[s]                      = IceMelt_device(vs).pres;
      rhofaci[s]                   = IceMelt_device(vs).rhofaci;
      table_val_qi2qr_melting[s]   = IceMelt_device(vs).table_val_qi2qr_melting;
      table_val_qi2qr_vent_melt[s] = IceMelt_device(vs).table_val_qi2qr_vent_melt;
      dv[s]                        = IceMelt_device(vs).dv;
      sc[s]                        = IceMelt_device(vs).sc;
      mu[s]                        = IceMelt_device(vs).mu;
      kap[s]                       = IceMelt_device(vs).kap;
      qv[s]                        = IceMelt_device(vs).qv;
      qi_incld[s]                  = IceMelt_device(vs).qi_incld;
      ni_incld[s]                  = IceMelt_device(vs).ni_incld;
      qi2qr_melt_tend[s]           = IceMelt_device(vs).qi2qr_melt_tend;
      ni2nr_melt_tend[s]           = IceMelt_device(vs).ni2nr_melt_tend;
    }

    Functions::ice_melting(rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,dv,sc,mu,kap,qv,qi_incld,ni_incld,qi2qr_melt_tend,ni2nr_melt_tend);
    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      IceMelt_device(vs).qi2qr_melt_tend = qi2qr_melt_tend[s];
      IceMelt_device(vs).ni2nr_melt_tend = ni2nr_melt_tend[s];
    }

  });

  // Sync back to host
  Kokkos::deep_copy(IceMelt_host, IceMelt_device);

  // Validate results
  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(IceMelt[s].qi2qr_melt_tend == IceMelt_host(s).qi2qr_melt_tend);
      REQUIRE(IceMelt[s].ni2nr_melt_tend == IceMelt_host(s).ni2nr_melt_tend);
    }
  }
  else if (this->m_baseline_action == GENERATE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      IceMelt_host(s).write(Base::m_fid);
    }
  }

}

}; // UnitWrap

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_ice_melting_test", "[p3_ice_melting_test]") {
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3IceMelting;

  T t;
  t.ice_melting_bfb();
}

} // namespace
