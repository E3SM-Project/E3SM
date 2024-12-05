#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {
/*
 * Unit-tests for p3_functions.
 */
template <typename D>
struct UnitWrap::UnitTest<D>::TestP3Conservation : public UnitWrap::UnitTest<D>::Base
{

  void cloud_water_conservation_tests_device() {

    using KTH = KokkosTypes<HostDevice>;

    CloudWaterConservationData cwdc[1] = {{sp(1e-5), sp(1.1), sp(1e-4), 0.0, 0.0, 0.0, 0.0, 0.0, sp(1.0), sp(1.0)}};

    // Sync to device
    KTH::view_1d<CloudWaterConservationData> cwdc_host("cwdc_host", 1);
    view_1d<CloudWaterConservationData> cwdc_device("cwdc_host", 1);

    // This copy only copies the input variables.
    std::copy(&cwdc[0], &cwdc[0] + 1, cwdc_host.data());
    Kokkos::deep_copy(cwdc_device, cwdc_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      Spack qc(cwdc_device(0).qc);
      Spack qc2qr_autoconv_tend(cwdc_device(0).qc2qr_autoconv_tend);
      Spack qc2qr_accret_tend(cwdc_device(0).qc2qr_accret_tend);
      Spack qc2qi_collect_tend(cwdc_device(0).qc2qi_collect_tend);
      Spack qc2qi_hetero_freeze_tend(cwdc_device(0).qc2qi_hetero_freeze_tend);
      Spack qc2qr_ice_shed_tend(cwdc_device(0).qc2qr_ice_shed_tend);
      Spack qc2qi_berg_tend(cwdc_device(0).qc2qi_berg_tend);
      Spack qi2qv_sublim_tend(cwdc_device(0).qi2qv_sublim_tend);
      Spack qv2qi_vapdep_tend(cwdc_device(0).qv2qi_vapdep_tend);

      Functions::cloud_water_conservation(qc, cwdc_device(0).dt, qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend);

      cwdc_device(0).qc = qc[0];
      cwdc_device(0).qc2qr_autoconv_tend = qc2qr_autoconv_tend[0];
      cwdc_device(0).qc2qr_accret_tend = qc2qr_accret_tend[0];
      cwdc_device(0).qc2qi_collect_tend = qc2qi_collect_tend[0];
      cwdc_device(0).qc2qi_hetero_freeze_tend = qc2qi_hetero_freeze_tend[0];
      cwdc_device(0).qc2qr_ice_shed_tend = qc2qr_ice_shed_tend[0];
      cwdc_device(0).qc2qi_berg_tend = qc2qi_berg_tend[0];
      cwdc_device(0).qi2qv_sublim_tend = qi2qv_sublim_tend[0];
      cwdc_device(0).qv2qi_vapdep_tend = qv2qi_vapdep_tend[0];
    });

    // Sync back to host
    Kokkos::deep_copy(cwdc_host, cwdc_device);

    const auto ratio = cwdc[0].qc/(cwdc[0].qc2qr_autoconv_tend * cwdc[0].dt);
    REQUIRE(std::abs(cwdc_host(0).qc2qr_autoconv_tend - cwdc[0].qc2qr_autoconv_tend*ratio) <= C::macheps);
    REQUIRE(cwdc_host(0).qc2qr_accret_tend == 0.0);
    REQUIRE(cwdc_host(0).qc2qi_collect_tend == 0.0);
    REQUIRE(cwdc_host(0).qc2qi_hetero_freeze_tend == 0.0);
    REQUIRE(cwdc_host(0).qc2qr_ice_shed_tend == 0.0);
    REQUIRE(cwdc_host(0).qc2qi_berg_tend == 0.0);
    REQUIRE(std::abs(cwdc_host(0).qi2qv_sublim_tend -(1.0 - ratio)) <= C::macheps);
    REQUIRE(std::abs(cwdc_host(0).qv2qi_vapdep_tend - (1.0 - ratio)) <= C::macheps);
    REQUIRE(cwdc_host[0].qc2qr_autoconv_tend * cwdc[0].dt <= cwdc_host[0].qc);
  }

  void rain_water_conservation_tests_device() {
    using KTH = KokkosTypes<HostDevice>;

    RainWaterConservationData rwdc[1] = {{sp(1e-5), 0.0, 0.0, 0.0, 0.0, sp(1.1), sp(1e-4), 0.0, 0.0 }};

    // Sync to device
    KTH::view_1d<RainWaterConservationData> rwdc_host("rwdc_host", 1);
    view_1d<RainWaterConservationData> rwdc_device("rwdc_host", 1);

    // This copy only copies the input variables.
    std::copy(&rwdc[0], &rwdc[0] + 1, rwdc_host.data());
    Kokkos::deep_copy(rwdc_device, rwdc_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      Spack qr(rwdc_device(0).qr);
      Spack qc2qr_autoconv_tend(rwdc_device(0).qc2qr_autoconv_tend);
      Spack qc2qr_accret_tend(rwdc_device(0).qc2qr_accret_tend);
      Spack qi2qr_melt_tend(rwdc_device(0).qi2qr_melt_tend);
      Spack qc2qr_ice_shed_tend(rwdc_device(0).qc2qr_ice_shed_tend);
      Spack qr2qv_evap_tend(rwdc_device(0).qr2qv_evap_tend);
      Spack qr2qi_collect_tend(rwdc_device(0).qr2qi_collect_tend);
      Spack qr2qi_immers_freeze_tend(rwdc_device(0).qr2qi_immers_freeze_tend);

      Functions::rain_water_conservation(qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, rwdc_device(0).dt, qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend);

      rwdc_device(0).qr = qr[0];
      rwdc_device(0).qc2qr_autoconv_tend = qc2qr_autoconv_tend[0];
      rwdc_device(0).qc2qr_accret_tend = qc2qr_accret_tend[0];
      rwdc_device(0).qi2qr_melt_tend = qi2qr_melt_tend[0];
      rwdc_device(0).qc2qr_ice_shed_tend = qc2qr_ice_shed_tend[0];
      rwdc_device(0).qr2qv_evap_tend = qr2qv_evap_tend[0];
      rwdc_device(0).qr2qi_collect_tend = qr2qi_collect_tend[0];
      rwdc_device(0).qr2qi_immers_freeze_tend = qr2qi_immers_freeze_tend[0];
    });

    // Sync back to host
    Kokkos::deep_copy(rwdc_host, rwdc_device);
    const auto ratio = rwdc[0].qr/(rwdc[0].qr2qv_evap_tend * rwdc[0].dt);

    //Here we check cases where source > sinks and sinks > 1e-20
    REQUIRE(rwdc_host(0).qc2qr_autoconv_tend == 0.0);
    REQUIRE(rwdc_host(0).qc2qr_accret_tend == 0.0);
    REQUIRE(rwdc_host(0).qi2qr_melt_tend == 0.0);
    REQUIRE(rwdc_host(0).qc2qr_ice_shed_tend == 0.0);

    //Check the value of qr2qv_evap_tend
    REQUIRE(std::abs(rwdc_host(0).qr2qv_evap_tend- rwdc[0].qr2qv_evap_tend*ratio)<= C::macheps);

    //Now test that conservation has actually been enforced
    REQUIRE( rwdc_host(0).qr2qv_evap_tend * rwdc_host(0).dt  <= rwdc_host(0).qr);
  }

  void ice_water_conservation_tests_device() {
    using KTH = KokkosTypes<HostDevice>;

    IceWaterConservationData iwdc[1] = {{sp(1e-5), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sp(1.1), sp(1e-4), 0.0}};

    // Sync to device
    KTH::view_1d<IceWaterConservationData> iwdc_host("iwdc_host", 1);
    view_1d<IceWaterConservationData> iwdc_device("iwdc_host", 1);

    // This copy only copies the input variables.
    std::copy(&iwdc[0], &iwdc[0] + 1, iwdc_host.data());
    Kokkos::deep_copy(iwdc_device, iwdc_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      Spack qi(iwdc_device(0).qi);
      Spack qv2qi_vapdep_tend(iwdc_device(0).qv2qi_vapdep_tend);
      Spack qv2qi_nucleat_tend(iwdc_device(0).qv2qi_nucleat_tend);
      Spack qr2qi_collect_tend(iwdc_device(0).qr2qi_collect_tend);
      Spack qc2qi_collect_tend(iwdc_device(0).qc2qi_collect_tend);
      Spack qr2qi_immers_freeze_tend(iwdc_device(0).qr2qi_immers_freeze_tend);
      Spack qc2qi_hetero_freeze_tend(iwdc_device(0).qc2qi_hetero_freeze_tend);
      Spack qc2qi_berg_tend(iwdc_device(0).qc2qi_berg_tend);
      Spack qi2qv_sublim_tend(iwdc_device(0).qi2qv_sublim_tend);
      Spack qi2qr_melt_tend(iwdc_device(0).qi2qr_melt_tend);

      Functions::ice_water_conservation(qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, qc2qi_berg_tend, iwdc_device(0).dt, qi2qv_sublim_tend, qi2qr_melt_tend);

      iwdc_device(0).qi = qi[0];
      iwdc_device(0).qv2qi_vapdep_tend = qv2qi_vapdep_tend[0];
      iwdc_device(0).qv2qi_nucleat_tend = qv2qi_nucleat_tend[0];
      iwdc_device(0).qr2qi_collect_tend = qr2qi_collect_tend[0];
      iwdc_device(0).qc2qi_collect_tend = qc2qi_collect_tend[0];
      iwdc_device(0).qr2qi_immers_freeze_tend = qr2qi_immers_freeze_tend[0];
      iwdc_device(0).qc2qi_hetero_freeze_tend = qc2qi_hetero_freeze_tend[0];
      iwdc_device(0).qc2qi_berg_tend = qc2qi_berg_tend[0];
      iwdc_device(0).qi2qv_sublim_tend = qi2qv_sublim_tend[0];
      iwdc_device(0).qi2qr_melt_tend = qi2qr_melt_tend[0];
    });

  }

  void run()
  {
    cloud_water_conservation_tests_device();

    rain_water_conservation_tests_device();

    ice_water_conservation_tests_device();
  }

  void cloud_water_conservation_unit_bfb_tests() {

    using KTH = KokkosTypes<HostDevice>;

    // These static asserts are important for many tests. If this test gets
    // removed, please put these lines in another test.
    static_assert(Spack::n <= max_pack_size,     "Unit testing infrastructure does not support this pack size (too big)");
    static_assert(max_pack_size % Spack::n == 0, "Unit testing infrastructure does not support this pack size (does not evenly divide 16)");

    CloudWaterConservationData cwdc[max_pack_size] = {
      //qc, cwdc_device(0).dt, qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend
      {9.9999999999999995e-7, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7},

      {9.9999999999999995e-7, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7},

      {9.9999999999999995e-7, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7},

      {9.9999999999999995e-7, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7}
    };

    // Sync to device
    KTH::view_1d<CloudWaterConservationData> cwdc_host("cwdc_host", max_pack_size);
    view_1d<CloudWaterConservationData> cwdc_device("cwdc_host", max_pack_size);

    // This copy only copies the input variables.
    std::copy(&cwdc[0], &cwdc[0] + max_pack_size, cwdc_host.data());
    Kokkos::deep_copy(cwdc_device, cwdc_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        cwdc[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qc, qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qc[s]     = cwdc_device(vs).qc;
        qc2qr_autoconv_tend[s]  = cwdc_device(vs).qc2qr_autoconv_tend;
        qc2qr_accret_tend[s]  = cwdc_device(vs).qc2qr_accret_tend;
        qc2qi_collect_tend[s]  = cwdc_device(vs).qc2qi_collect_tend;
        qc2qi_hetero_freeze_tend[s] = cwdc_device(vs).qc2qi_hetero_freeze_tend;
        qc2qr_ice_shed_tend[s]  = cwdc_device(vs).qc2qr_ice_shed_tend;
        qc2qi_berg_tend[s] = cwdc_device(vs).qc2qi_berg_tend;
        qi2qv_sublim_tend[s]  = cwdc_device(vs).qi2qv_sublim_tend;
        qv2qi_vapdep_tend[s]  = cwdc_device(vs).qv2qi_vapdep_tend;
      }

      Functions::cloud_water_conservation(qc, cwdc_device(0).dt, qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend);
      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cwdc_device(vs).qc     = qc[s];
        cwdc_device(vs).qc2qr_autoconv_tend  = qc2qr_autoconv_tend[s];
        cwdc_device(vs).qc2qr_accret_tend  = qc2qr_accret_tend[s];
        cwdc_device(vs).qc2qi_collect_tend  = qc2qi_collect_tend[s];
        cwdc_device(vs).qc2qi_hetero_freeze_tend = qc2qi_hetero_freeze_tend[s];
        cwdc_device(vs).qc2qi_berg_tend = qc2qi_berg_tend[s];
        cwdc_device(vs).qi2qv_sublim_tend  = qi2qv_sublim_tend[s];
        cwdc_device(vs).qv2qi_vapdep_tend  = qv2qi_vapdep_tend[s];
      }

    });
    // Sync back to host
    Kokkos::deep_copy(cwdc_host, cwdc_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(cwdc[s].qc     == cwdc_host(s).qc);
        REQUIRE(cwdc[s].qc2qr_autoconv_tend  == cwdc_host(s).qc2qr_autoconv_tend);
        REQUIRE(cwdc[s].qc2qr_accret_tend  == cwdc_host(s).qc2qr_accret_tend);
        REQUIRE(cwdc[s].qc2qi_collect_tend  == cwdc_host(s).qc2qi_collect_tend);
        REQUIRE(cwdc[s].qc2qi_hetero_freeze_tend == cwdc_host(s).qc2qi_hetero_freeze_tend);
        REQUIRE(cwdc[s].qc2qi_berg_tend == cwdc_host(s).qc2qi_berg_tend);
        REQUIRE(cwdc[s].qi2qv_sublim_tend  == cwdc_host(s).qi2qv_sublim_tend);
        REQUIRE(cwdc[s].qv2qi_vapdep_tend  == cwdc_host(s).qv2qi_vapdep_tend);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        cwdc_host(s).write(Base::m_fid);
      }
    }
  }

  void ice_water_conservation_unit_bfb_tests()
  {
    using KTH = KokkosTypes<HostDevice>;

    IceWaterConservationData iwdc[max_pack_size] = {
      // qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, iwdc_device(0).dt, qi2qv_sublim_tend, qi2qr_melt_tend
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0},

      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0},

      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0},

      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0}
    };

    // Sync to device
    KTH::view_1d<IceWaterConservationData> iwdc_host("iwdc_host", max_pack_size);
    view_1d<IceWaterConservationData> iwdc_device("iwdc_host", max_pack_size);

    // This copy only copies the input variables.
    std::copy(&iwdc[0], &iwdc[0] + max_pack_size, iwdc_host.data());
    Kokkos::deep_copy(iwdc_device, iwdc_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        iwdc[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qi,qv2qi_vapdep_tend,qv2qi_nucleat_tend,qc2qi_berg_tend,qr2qi_collect_tend,qc2qi_collect_tend,qr2qi_immers_freeze_tend,qc2qi_hetero_freeze_tend,qi2qv_sublim_tend,qi2qr_melt_tend;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qi[s]  = iwdc_device(vs).qi;
        qv2qi_vapdep_tend[s]  = iwdc_device(vs).qv2qi_vapdep_tend;
        qv2qi_nucleat_tend[s]  = iwdc_device(vs).qv2qi_nucleat_tend;
        qc2qi_berg_tend[s] = iwdc_device(vs).qc2qi_berg_tend;
        qr2qi_collect_tend[s]  = iwdc_device(vs).qr2qi_collect_tend;
        qc2qi_collect_tend[s]  = iwdc_device(vs).qc2qi_collect_tend;
        qr2qi_immers_freeze_tend[s] = iwdc_device(vs).qr2qi_immers_freeze_tend;
        qc2qi_hetero_freeze_tend[s] = iwdc_device(vs).qc2qi_hetero_freeze_tend;
        qi2qv_sublim_tend[s] = iwdc_device(vs).qi2qv_sublim_tend;
        qi2qr_melt_tend[s] = iwdc_device(vs).qi2qr_melt_tend;
      }

      Functions::ice_water_conservation(qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, iwdc_device(0).dt, qi2qv_sublim_tend, qi2qr_melt_tend);
      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        iwdc_device(vs).qi = qi[s];
        iwdc_device(vs).qv2qi_vapdep_tend = qv2qi_vapdep_tend[s];
        iwdc_device(vs).qv2qi_nucleat_tend = qv2qi_nucleat_tend[s];
        iwdc_device(vs).qc2qi_berg_tend = qc2qi_berg_tend[s];
        iwdc_device(vs).qr2qi_collect_tend = qr2qi_collect_tend[s];
        iwdc_device(vs).qc2qi_collect_tend = qc2qi_collect_tend[s];
        iwdc_device(vs).qr2qi_immers_freeze_tend = qr2qi_immers_freeze_tend[s];
        iwdc_device(vs).qc2qi_hetero_freeze_tend = qc2qi_hetero_freeze_tend[s];
        iwdc_device(vs).qi2qv_sublim_tend = qi2qv_sublim_tend[s];
        iwdc_device(vs).qi2qr_melt_tend = qi2qr_melt_tend[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(iwdc_host, iwdc_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(iwdc[s].qi == iwdc_host(s).qi);
        REQUIRE(iwdc[s].qv2qi_vapdep_tend == iwdc_host(s).qv2qi_vapdep_tend );
        REQUIRE(iwdc[s].qv2qi_nucleat_tend == iwdc_host(s).qv2qi_nucleat_tend);
        REQUIRE(iwdc[s].qc2qi_berg_tend == iwdc_host(s).qc2qi_berg_tend);
        REQUIRE(iwdc[s].qr2qi_collect_tend  == iwdc_host(s).qr2qi_collect_tend);
        REQUIRE(iwdc[s].qc2qi_collect_tend == iwdc_host(s).qc2qi_collect_tend);
        REQUIRE(iwdc[s].qr2qi_immers_freeze_tend == iwdc_host(s).qr2qi_immers_freeze_tend);
        REQUIRE(iwdc[s].qc2qi_hetero_freeze_tend == iwdc_host(s).qc2qi_hetero_freeze_tend);
        REQUIRE(iwdc[s].qi2qv_sublim_tend == iwdc_host(s).qi2qv_sublim_tend);
        REQUIRE(iwdc[s].qi2qr_melt_tend == iwdc_host(s).qi2qr_melt_tend);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        iwdc_host(s).write(Base::m_fid);
      }
    }
  }

  void rain_water_conservation_unit_bfb_tests() {

    using KTH = KokkosTypes<HostDevice>;

    RainWaterConservationData rwdc[max_pack_size] = {
      // qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, rwdc_device(0).dt, qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend
      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10},

      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10},

      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10},

      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10}
    };

    // Sync to device
    KTH::view_1d<RainWaterConservationData> rwdc_host("rwdc_host", max_pack_size);
    view_1d<RainWaterConservationData> rwdc_device("rwdc_host", max_pack_size);

    // This copy only copies the input variables.
    std::copy(&rwdc[0], &rwdc[0] + max_pack_size, rwdc_host.data());
    Kokkos::deep_copy(rwdc_device, rwdc_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        rwdc[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qr[s]     = rwdc_device(vs).qr;
        qc2qr_autoconv_tend[s]  = rwdc_device(vs).qc2qr_autoconv_tend;
        qc2qr_accret_tend[s]  = rwdc_device(vs).qc2qr_accret_tend;
        qi2qr_melt_tend[s]  = rwdc_device(vs).qi2qr_melt_tend;
        qc2qr_ice_shed_tend[s]  = rwdc_device(vs).qc2qr_ice_shed_tend;
        qr2qv_evap_tend[s]  = rwdc_device(vs).qr2qv_evap_tend;
        qr2qi_collect_tend[s]  = rwdc_device(vs).qr2qi_collect_tend;
        qr2qi_immers_freeze_tend[s] = rwdc_device(vs).qr2qi_immers_freeze_tend;
      }

      Functions::rain_water_conservation(qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, rwdc_device(0).dt, qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend);
      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rwdc_device(vs).qr     = qr[s];
        rwdc_device(vs).qc2qr_autoconv_tend  = qc2qr_autoconv_tend[s];
        rwdc_device(vs).qc2qr_accret_tend  = qc2qr_accret_tend[s];
        rwdc_device(vs).qi2qr_melt_tend  = qi2qr_melt_tend[s];
        rwdc_device(vs).qc2qr_ice_shed_tend  = qc2qr_ice_shed_tend[s];
        rwdc_device(vs).qr2qv_evap_tend  = qr2qv_evap_tend[s];
        rwdc_device(vs).qr2qi_collect_tend  = qr2qi_collect_tend[s];
        rwdc_device(vs).qr2qi_immers_freeze_tend = qr2qi_immers_freeze_tend[s];
      }

    });

    // Sync back to host
    Kokkos::deep_copy(rwdc_host, rwdc_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(rwdc[s].qr     == rwdc_host(s).qr);
        REQUIRE(rwdc[s].qc2qr_autoconv_tend  == rwdc_host(s).qc2qr_autoconv_tend);
        REQUIRE(rwdc[s].qc2qr_accret_tend  == rwdc_host(s).qc2qr_accret_tend);
        REQUIRE(rwdc[s].qi2qr_melt_tend  == rwdc_host(s).qi2qr_melt_tend);
        REQUIRE(rwdc[s].qc2qr_ice_shed_tend  == rwdc_host(s).qc2qr_ice_shed_tend);
        REQUIRE(rwdc[s].qr2qv_evap_tend  == rwdc_host(s).qr2qv_evap_tend);
        REQUIRE(rwdc[s].qr2qi_collect_tend  == rwdc_host(s).qr2qi_collect_tend);
        REQUIRE(rwdc[s].qr2qi_immers_freeze_tend == rwdc_host(s).qr2qi_immers_freeze_tend);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        rwdc_host(s).write(Base::m_fid);
      }
    }
  }

  void run_bfb() {
    cloud_water_conservation_unit_bfb_tests();

    rain_water_conservation_unit_bfb_tests();

    ice_water_conservation_unit_bfb_tests();
  }

};

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3UpdatePrognosticIce : public UnitWrap::UnitTest<D>::Base
{
  void update_prognostic_ice_unit_bfb_tests() {

    constexpr Scalar nmltratio     = C::nmltratio;
    constexpr Scalar dt            = 1.8000E+03;
    constexpr bool   do_predict_nc = true;
    constexpr Scalar latvap        = C::LatVap;
    constexpr Scalar latice        = C::LatIce;

    //baseline generated data is input to the following
    P3UpdatePrognosticIceData pupidc[max_pack_size] = {

      {4.9078E-19, 1.5312E-09, 4.4387E-09, 3.7961E+06, 1.7737E-04, 0.0000E+00, 3.8085E-08, 5.1281E+04, 1.9251E-15,
       3.4778E-04, 3.5801E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.1386E-07, 0.0000E+00, 0.0000E+00, 2.7053E-02,
       0.0000E+00, 1.9209E-10, 1.0686E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       4.5312E+02, 2.8720E+02, 5.0000E-03, 6.4286E-05, 1.2344E+08, 7.3684E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       6.4286E-05, 1.0000E-02},

      {2.1097E-18, 2.7648E-09, 3.8261E-09, 3.7754E+06, 6.8685E-04, 0.0000E+00, 4.1018E-08, 5.1227E+04, 4.8876E-15,
       1.3468E-03, 2.8059E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 7.1049E-07, 0.0000E+00, 0.0000E+00, 2.4547E-02,
       0.0000E+00, 2.8615E-10, 1.0741E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       3.4890E+02, 2.8642E+02, 5.0000E-03, 7.1429E-05, 1.2345E+08, 7.8947E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.1429E-05, 1.0000E-02},

      {8.9820E-18, 4.2529E-09, 2.9520E-09, 3.7537E+06, 2.6598E-03, 0.0000E+00, 4.3700E-08, 5.1171E+04, 1.4266E-14,
       5.2153E-03, 1.9880E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 9.0244E-07, 0.0000E+00, 0.0000E+00, 2.1083E-02,
       0.0000E+00, 3.7631E-10, 1.0796E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.8656E+02, 2.8565E+02, 5.0000E-03, 7.8571E-05, 1.2345E+08, 8.4211E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.8571E-05, 1.0000E-02},

      {3.7942E-17, 6.0115E-09, 1.8004E-09, 3.7310E+06, 1.0300E-02, 0.0000E+00, 4.6119E-08, 5.1112E+04, 4.4518E-14,
       2.0196E-02, 1.1226E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 1.0879E-06, 0.0000E+00, 0.0000E+00, 1.7646E-02,
       0.0000E+00, 4.5891E-10, 1.0853E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.4570E+02, 2.8489E+02, 5.0000E-03, 8.5714E-05, 1.2345E+08, 8.9474E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       8.5714E-05, 1.0000E-02},

      {4.9078E-19, 1.5312E-09, 4.4387E-09, 3.7961E+06, 1.7737E-04, 0.0000E+00, 3.8085E-08, 5.1281E+04, 1.9251E-15,
       3.4778E-04, 3.5801E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.1386E-07, 0.0000E+00, 0.0000E+00, 2.7053E-02,
       0.0000E+00, 1.9209E-10, 1.0686E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       4.5312E+02, 2.8720E+02, 5.0000E-03, 6.4286E-05, 1.2344E+08, 7.3684E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       6.4286E-05, 1.0000E-02},

      {2.1097E-18, 2.7648E-09, 3.8261E-09, 3.7754E+06, 6.8685E-04, 0.0000E+00, 4.1018E-08, 5.1227E+04, 4.8876E-15,
       1.3468E-03, 2.8059E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 7.1049E-07, 0.0000E+00, 0.0000E+00, 2.4547E-02,
       0.0000E+00, 2.8615E-10, 1.0741E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       3.4890E+02, 2.8642E+02, 5.0000E-03, 7.1429E-05, 1.2345E+08, 7.8947E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.1429E-05, 1.0000E-02},

      {8.9820E-18, 4.2529E-09, 2.9520E-09, 3.7537E+06, 2.6598E-03, 0.0000E+00, 4.3700E-08, 5.1171E+04, 1.4266E-14,
       5.2153E-03, 1.9880E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 9.0244E-07, 0.0000E+00, 0.0000E+00, 2.1083E-02,
       0.0000E+00, 3.7631E-10, 1.0796E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.8656E+02, 2.8565E+02, 5.0000E-03, 7.8571E-05, 1.2345E+08, 8.4211E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.8571E-05, 1.0000E-02},

      {3.7942E-17, 6.0115E-09, 1.8004E-09, 3.7310E+06, 1.0300E-02, 0.0000E+00, 4.6119E-08, 5.1112E+04, 4.4518E-14,
       2.0196E-02, 1.1226E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 1.0879E-06, 0.0000E+00, 0.0000E+00, 1.7646E-02,
       0.0000E+00, 4.5891E-10, 1.0853E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.4570E+02, 2.8489E+02, 5.0000E-03, 8.5714E-05, 1.2345E+08, 8.9474E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       8.5714E-05, 1.0000E-02},

      {4.9078E-19, 1.5312E-09, 4.4387E-09, 3.7961E+06, 1.7737E-04, 0.0000E+00, 3.8085E-08, 5.1281E+04, 1.9251E-15,
       3.4778E-04, 3.5801E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.1386E-07, 0.0000E+00, 0.0000E+00, 2.7053E-02,
       0.0000E+00, 1.9209E-10, 1.0686E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       4.5312E+02, 2.8720E+02, 5.0000E-03, 6.4286E-05, 1.2344E+08, 7.3684E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       6.4286E-05, 1.0000E-02},

      {2.1097E-18, 2.7648E-09, 3.8261E-09, 3.7754E+06, 6.8685E-04, 0.0000E+00, 4.1018E-08, 5.1227E+04, 4.8876E-15,
       1.3468E-03, 2.8059E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 7.1049E-07, 0.0000E+00, 0.0000E+00, 2.4547E-02,
       0.0000E+00, 2.8615E-10, 1.0741E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       3.4890E+02, 2.8642E+02, 5.0000E-03, 7.1429E-05, 1.2345E+08, 7.8947E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.1429E-05, 1.0000E-02},

      {8.9820E-18, 4.2529E-09, 2.9520E-09, 3.7537E+06, 2.6598E-03, 0.0000E+00, 4.3700E-08, 5.1171E+04, 1.4266E-14,
       5.2153E-03, 1.9880E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 9.0244E-07, 0.0000E+00, 0.0000E+00, 2.1083E-02,
       0.0000E+00, 3.7631E-10, 1.0796E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.8656E+02, 2.8565E+02, 5.0000E-03, 7.8571E-05, 1.2345E+08, 8.4211E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.8571E-05, 1.0000E-02},

      {3.7942E-17, 6.0115E-09, 1.8004E-09, 3.7310E+06, 1.0300E-02, 0.0000E+00, 4.6119E-08, 5.1112E+04, 4.4518E-14,
       2.0196E-02, 1.1226E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 1.0879E-06, 0.0000E+00, 0.0000E+00, 1.7646E-02,
       0.0000E+00, 4.5891E-10, 1.0853E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.4570E+02, 2.8489E+02, 5.0000E-03, 8.5714E-05, 1.2345E+08, 8.9474E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       8.5714E-05, 1.0000E-02},

      {4.9078E-19, 1.5312E-09, 4.4387E-09, 3.7961E+06, 1.7737E-04, 0.0000E+00, 3.8085E-08, 5.1281E+04, 1.9251E-15,
       3.4778E-04, 3.5801E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.1386E-07, 0.0000E+00, 0.0000E+00, 2.7053E-02,
       0.0000E+00, 1.9209E-10, 1.0686E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       4.5312E+02, 2.8720E+02, 5.0000E-03, 6.4286E-05, 1.2344E+08, 7.3684E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       6.4286E-05, 1.0000E-02},

      {2.1097E-18, 2.7648E-09, 3.8261E-09, 3.7754E+06, 6.8685E-04, 0.0000E+00, 4.1018E-08, 5.1227E+04, 4.8876E-15,
       1.3468E-03, 2.8059E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 7.1049E-07, 0.0000E+00, 0.0000E+00, 2.4547E-02,
       0.0000E+00, 2.8615E-10, 1.0741E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       3.4890E+02, 2.8642E+02, 5.0000E-03, 7.1429E-05, 1.2345E+08, 7.8947E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.1429E-05, 1.0000E-02},

      {8.9820E-18, 4.2529E-09, 2.9520E-09, 3.7537E+06, 2.6598E-03, 0.0000E+00, 4.3700E-08, 5.1171E+04, 1.4266E-14,
       5.2153E-03, 1.9880E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 9.0244E-07, 0.0000E+00, 0.0000E+00, 2.1083E-02,
       0.0000E+00, 3.7631E-10, 1.0796E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.8656E+02, 2.8565E+02, 5.0000E-03, 7.8571E-05, 1.2345E+08, 8.4211E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       7.8571E-05, 1.0000E-02},

      {3.7942E-17, 6.0115E-09, 1.8004E-09, 3.7310E+06, 1.0300E-02, 0.0000E+00, 4.6119E-08, 5.1112E+04, 4.4518E-14,
       2.0196E-02, 1.1226E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 1.0879E-06, 0.0000E+00, 0.0000E+00, 1.7646E-02,
       0.0000E+00, 4.5891E-10, 1.0853E+00, latvap+latice, latice, do_predict_nc,    true,         dt, nmltratio,
       2.4570E+02, 2.8489E+02, 5.0000E-03, 8.5714E-05, 1.2345E+08, 8.9474E-06, 1.0000E+06, 1.0000E-04, 1.0000E+06,
       8.5714E-05, 1.0000E-02},
    };

    // Sync to device
    view_1d<P3UpdatePrognosticIceData> pupidc_device("pupidc", max_pack_size);
    auto pupidc_host = Kokkos::create_mirror_view(pupidc_device);

    // This copy only copies the input variables.
    std::copy(&pupidc[0], &pupidc[0] + max_pack_size, pupidc_host.data());
    Kokkos::deep_copy(pupidc_device, pupidc_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        pupidc[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend,
            qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend, ni2nr_melt_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend,
            ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend, qc2qi_berg_tend, inv_exner,
            rho_qm_cloud, th_atm, qv, qc, nc, qr, nr, qi, ni, qm, bm;
      Scalar dt;
      bool do_predict_nc;
      Smask log_wetgrowth;

      // variables with single values assigned outside of the for loop
      dt            = pupidc_device(0).dt;
      do_predict_nc = pupidc_device(0).do_predict_nc;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qc2qi_hetero_freeze_tend[s] = pupidc_device(vs).qc2qi_hetero_freeze_tend;
        qc2qi_collect_tend[s]  = pupidc_device(vs).qc2qi_collect_tend;
        qc2qr_ice_shed_tend[s]  = pupidc_device(vs).qc2qr_ice_shed_tend;
        nc_collect_tend[s]  = pupidc_device(vs).nc_collect_tend;
        nc2ni_immers_freeze_tend[s] = pupidc_device(vs).nc2ni_immers_freeze_tend;
        ncshdc[s] = pupidc_device(vs).ncshdc;
        qr2qi_collect_tend[s]  = pupidc_device(vs).qr2qi_collect_tend;
        nr_collect_tend[s]  = pupidc_device(vs).nr_collect_tend;
        qr2qi_immers_freeze_tend[s] = pupidc_device(vs).qr2qi_immers_freeze_tend;
        nr2ni_immers_freeze_tend[s] = pupidc_device(vs).nr2ni_immers_freeze_tend;
        nr_ice_shed_tend[s] = pupidc_device(vs).nr_ice_shed_tend;
        qi2qr_melt_tend[s]  = pupidc_device(vs).qi2qr_melt_tend;
        ni2nr_melt_tend[s]  = pupidc_device(vs).ni2nr_melt_tend;
        qi2qv_sublim_tend[s]  = pupidc_device(vs).qi2qv_sublim_tend;
        qv2qi_vapdep_tend[s]  = pupidc_device(vs).qv2qi_vapdep_tend;
        qv2qi_nucleat_tend[s]  = pupidc_device(vs).qv2qi_nucleat_tend;
        ni_nucleat_tend[s]  = pupidc_device(vs).ni_nucleat_tend;
        ni_selfcollect_tend[s]  = pupidc_device(vs).ni_selfcollect_tend;
        ni_sublim_tend[s]  = pupidc_device(vs).ni_sublim_tend;
        qc2qi_berg_tend[s] = pupidc_device(vs).qc2qi_berg_tend;
        inv_exner[s]  = pupidc_device(vs).inv_exner;

        rho_qm_cloud[s] = pupidc_device(vs).rho_qm_cloud;
        th_atm[s]    = pupidc_device(vs).th_atm;
        qv[s]    = pupidc_device(vs).qv;
        qc[s]    = pupidc_device(vs).qc;
        nc[s]    = pupidc_device(vs).nc;
        qr[s]    = pupidc_device(vs).qr;
        nr[s]    = pupidc_device(vs).nr;
        qi[s] = pupidc_device(vs).qi;
        ni[s] = pupidc_device(vs).ni;
        qm[s] = pupidc_device(vs).qm;
        bm[s] = pupidc_device(vs).bm;

        log_wetgrowth.set(s, pupidc_device(vs).log_wetgrowth);
      }

      Functions::update_prognostic_ice(qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend,ncshdc,
                                       qr2qi_collect_tend,   nr_collect_tend,  qr2qi_immers_freeze_tend,  nr2ni_immers_freeze_tend,  nr_ice_shed_tend,
                                       qi2qr_melt_tend,  ni2nr_melt_tend,  qi2qv_sublim_tend,  qv2qi_vapdep_tend,  qv2qi_nucleat_tend,  ni_nucleat_tend,
                                       ni_selfcollect_tend,  ni_sublim_tend,  qc2qi_berg_tend,  inv_exner,
                                       do_predict_nc, log_wetgrowth,  dt,  pupidc_device(0).nmltratio,
                                       rho_qm_cloud, th_atm, qv, qi, ni, qm,
                                       bm, qc, nc, qr, nr);

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        pupidc_device(vs).th_atm    = th_atm[s];
        pupidc_device(vs).qv        = qv[s];
        pupidc_device(vs).qc        = qc[s];
        pupidc_device(vs).nc        = nc[s];
        pupidc_device(vs).qr        = qr[s];
        pupidc_device(vs).nr        = nr[s];
        pupidc_device(vs).qi        = qi[s];
        pupidc_device(vs).ni        = ni[s];
        pupidc_device(vs).qm        = qm[s];
        pupidc_device(vs).bm        = bm[s];
      }

    });

    // Sync back to host
    Kokkos::deep_copy(pupidc_host, pupidc_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(pupidc[s].th_atm    == pupidc_host(s).th_atm);
        REQUIRE(pupidc[s].qc        == pupidc_host(s).qc);
        REQUIRE(pupidc[s].nr        == pupidc_host(s).nr);
        REQUIRE(pupidc[s].qr        == pupidc_host(s).qr);
        REQUIRE(pupidc[s].qv        == pupidc_host(s).qv);
        REQUIRE(pupidc[s].nc        == pupidc_host(s).nc);
        REQUIRE(pupidc[s].qi        == pupidc_host(s).qi);
        REQUIRE(pupidc[s].ni        == pupidc_host(s).ni);
        REQUIRE(pupidc[s].qm        == pupidc_host(s).qm);
        REQUIRE(pupidc[s].bm        == pupidc_host(s).bm );
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        pupidc_host(s).write(Base::m_fid);
      }
    }
  }

  void run_bfb() {
    update_prognostic_ice_unit_bfb_tests();
  }

}; //TestP3UpdatePrognosticIce

template <typename D>
struct UnitWrap::UnitTest<D>::TestGetTimeSpacePhysVariables : public UnitWrap::UnitTest<D>::Base
{
  void get_time_space_phys_variables_unit_bfb_tests() {
    constexpr Scalar latvap = C::LatVap;
    constexpr Scalar latice = C::LatIce;

    //baseline generated data is input to the following
    GetTimeSpacePhysVarsData gtspvd[max_pack_size] = {
      //        T_atm,       pres,        rho,       latent_heat_vapor,       latent_heat_sublim,        qv_sat_l,        qv_sat_i
      {2.9792E+02, 9.8711E+04, 1.1532E+00, latvap, latvap+latice, 2.0321E-02, 2.0321E-02},
      {2.9792E+02, 9.8711E+04, 1.1532E+00, latvap, latvap+latice, 2.0321E-02, 2.0321E-02},
      {2.9583E+02, 9.7322E+04, 1.1449E+00, latvap, latvap+latice, 1.8120E-02, 1.8120E-02},
      {2.9375E+02, 9.5933E+04, 1.1366E+00, latvap, latvap+latice, 1.6134E-02, 1.6134E-02},
      {2.8959E+02, 9.3156E+04, 1.1196E+00, latvap, latvap+latice, 1.2729E-02, 1.2729E-02},
      {2.8750E+02, 9.1767E+04, 1.1109E+00, latvap, latvap+latice, 1.1279E-02, 1.1279E-02},
      {2.8542E+02, 9.0378E+04, 1.1020E+00, latvap, latvap+latice, 9.9759E-03, 9.9759E-03},
      {2.8334E+02, 8.8989E+04, 1.0931E+00, latvap, latvap+latice, 8.8076E-03, 8.8076E-03},
      {2.8125E+02, 8.7600E+04, 1.0840E+00, latvap, latvap+latice, 7.7615E-03, 7.7615E-03},
      {2.7917E+02, 8.6211E+04, 1.0748E+00, latvap, latvap+latice, 6.8265E-03, 6.8265E-03},
      {2.7709E+02, 8.4822E+04, 1.0654E+00, latvap, latvap+latice, 5.9921E-03, 5.9921E-03},
      {2.7501E+02, 8.3433E+04, 1.0559E+00, latvap, latvap+latice, 5.2488E-03, 5.2488E-03},
      {2.7292E+02, 8.2044E+04, 1.0463E+00, latvap, latvap+latice, 4.5879E-03, 4.5766E-03},
      {2.7084E+02, 8.0656E+04, 1.0365E+00, latvap, latvap+latice, 4.0015E-03, 3.9112E-03},
      {2.6876E+02, 7.9267E+04, 1.0265E+00, latvap, latvap+latice, 3.4821E-03, 3.3349E-03},
      {2.6667E+02, 7.7878E+04, 1.0164E+00, latvap, latvap+latice, 3.0231E-03, 2.8368E-03},
    };

    // Sync to device
    view_1d<GetTimeSpacePhysVarsData> gtspvd_device("gtspvd", max_pack_size);
    auto gtspvd_host = Kokkos::create_mirror_view(gtspvd_device);

    // This copy only copies the input variables.
    std::copy(&gtspvd[0], &gtspvd[0] + max_pack_size, gtspvd_host.data());
    Kokkos::deep_copy(gtspvd_device, gtspvd_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        gtspvd[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack T_atm, pres, rho, qv_sat_l, qv_sat_i, mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        T_atm[s]                = gtspvd_device(vs).T_atm;
        pres[s]                 = gtspvd_device(vs).pres;
        rho[s]                  = gtspvd_device(vs).rho;
        qv_sat_l[s]             = gtspvd_device(vs).qv_sat_l;
        qv_sat_i[s]             = gtspvd_device(vs).qv_sat_i;

        mu[s]     = gtspvd_device(vs).mu;
        dv[s]     = gtspvd_device(vs).dv;
        sc[s]     = gtspvd_device(vs).sc;
        dqsdt[s]  = gtspvd_device(vs).dqsdt;
        dqsidt[s] = gtspvd_device(vs).dqsidt;
        ab[s]     = gtspvd_device(vs).ab;
        abi[s]    = gtspvd_device(vs).abi;
        kap[s]    = gtspvd_device(vs).kap;
        eii[s]    = gtspvd_device(vs).eii;
      }

      Functions::get_time_space_phys_variables(T_atm, pres, rho, qv_sat_l, qv_sat_i, mu, dv, sc, dqsdt, dqsidt,
                                               ab, abi, kap, eii);

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        gtspvd_device(vs).T_atm                = T_atm[s];
        gtspvd_device(vs).pres                 = pres[s];
        gtspvd_device(vs).rho                  = rho[s];
        gtspvd_device(vs).qv_sat_l             = qv_sat_l[s];
        gtspvd_device(vs).qv_sat_i             = qv_sat_i[s];

        gtspvd_device(vs).mu     = mu[s];
        gtspvd_device(vs).dv     = dv[s];
        gtspvd_device(vs).sc     = sc[s];
        gtspvd_device(vs).dqsdt  = dqsdt[s];
        gtspvd_device(vs).dqsidt = dqsidt[s];
        gtspvd_device(vs).ab     = ab[s];
        gtspvd_device(vs).abi    = abi[s];
        gtspvd_device(vs).kap    = kap[s];
        gtspvd_device(vs).eii    = eii[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(gtspvd_host, gtspvd_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(gtspvd[s].mu     == gtspvd_host(s).mu);
        REQUIRE(gtspvd[s].dv     == gtspvd_host(s).dv);
        REQUIRE(gtspvd[s].sc     == gtspvd_host(s).sc);
        REQUIRE(gtspvd[s].dqsdt  == gtspvd_host(s).dqsdt);
        REQUIRE(gtspvd[s].dqsidt == gtspvd_host(s).dqsidt);
        REQUIRE(gtspvd[s].ab     == gtspvd_host(s).ab);
        REQUIRE(gtspvd[s].abi    == gtspvd_host(s).abi);
        REQUIRE(gtspvd[s].kap    == gtspvd_host(s).kap);
        REQUIRE(gtspvd[s].eii    == gtspvd_host(s).eii);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        gtspvd_host(s).write(Base::m_fid);
      }
    }
  }

  void run_bfb() {
    get_time_space_phys_variables_unit_bfb_tests();
  }
}; //TestGetTimeSpacePhysVariables

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3UpdatePrognosticLiq : public UnitWrap::UnitTest<D>::Base
{
  void update_prognostic_liquid_unit_bfb_tests() {
    constexpr Scalar latvap = C::LatVap;

    //baseline generated data is input to the following
    P3UpdatePrognosticLiqData pupldc[max_pack_size] = {

      {1.0631E-12, 1.0631E+00, 1.5833E-12, 1.5833E+00, 2.4190E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00, 4.2517E+00,
       true      , true      , 8.6718E-01, 1.0037E+00, latvap, 1.8000E+03, 2.9902E+02, 5.0000E-02, 1.0000E-06, 1.0000E+06, 1.0010E-06,
       6.3726E+05},

      {3.2784E-08, 1.8780E+07, 2.1753E-11, 1.2461E+04, 7.8657E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.8748E+04,
       true      , true      , 9.8387E-01, 1.0741E+00, latvap, 1.8000E+03, 2.9033E+02, 3.7211E-03, 5.9050E-05,-6.6723E+09,-5.9050E-05,
       -8.6159E+07},

      {3.2796E-09, 1.8778E+07, 1.8830E-12, 1.0782E+04, 6.8061E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 6.3698E+04,
       true      , true      , 9.0740E-01, 1.0293E+00, latvap, 1.8000E+03, 2.9376E+02, 5.0000E-03, 5.9067E-06,-6.9543E+09, 1.0439E-04,
       -1.6967E+07},

      {6.5634E-09, 1.8778E+07, 3.8238E-12, 1.0940E+04, 6.9061E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 6.3181E+04,
       true      , true      , 9.1484E-01, 1.0339E+00, latvap, 1.8000E+03, 2.9291E+02, 5.0000E-03, 1.1821E-05,-6.9282E+09, 1.0615E-04,
       -2.8223E+07},

      {9.8516E-09, 1.8779E+07, 5.8258E-12, 1.1105E+04, 7.0101E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 6.2655E+04,
       true      , true      , 9.2251E-01, 1.0386E+00, latvap, 1.8000E+03, 2.9206E+02, 5.0000E-03, 1.7743E-05,-6.9009E+09, 1.0790E-04,
       -3.9628E+07},

      {1.3145E-08, 1.8779E+07, 7.8929E-12, 1.1276E+04, 7.1180E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 6.2122E+04,
       true      , true      , 9.3043E-01, 1.0433E+00, latvap, 1.8000E+03, 2.9123E+02, 5.0000E-03, 2.3674E-05,-6.8725E+09, 1.0963E-04,
       -5.1189E+07},

      {1.6443E-08, 1.8779E+07, 1.0029E-11, 1.1454E+04, 7.2303E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 6.1581E+04,
       true      , true      , 9.3860E-01, 1.0482E+00, latvap, 1.8000E+03, 2.9040E+02, 5.0000E-03, 2.9615E-05,-6.8428E+09, 1.1136E-04,
       -6.2915E+07},

      {1.9746E-08, 1.8779E+07, 1.2238E-11, 1.1639E+04, 7.3471E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 6.1031E+04,
       true      , true      , 9.4705E-01, 1.0531E+00, latvap, 1.8000E+03, 2.8958E+02, 5.0000E-03, 3.5565E-05,-6.8117E+09, 1.1308E-04,
       -7.4813E+07},

      {2.3047E-08, 1.8779E+07, 1.4521E-11, 1.1832E+04, 7.4688E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 6.0474E+04,
       true      , true      , 9.5579E-01, 1.0582E+00, latvap, 1.8000E+03, 2.8941E+02, 4.7949E-03, 4.1510E-05,-6.7792E+09, 1.4787E-05,
       -8.2885E+07},

      {2.6289E-08, 1.8779E+07, 1.6845E-11, 1.2033E+04, 7.5955E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.9907E+04,
       true      , true      , 9.6483E-01, 1.0634E+00, latvap, 1.8000E+03, 2.8972E+02, 4.4341E-03, 4.7350E-05,-6.7452E+09,-4.7350E-05,
       -8.3634E+07},

      {2.9533E-08, 1.8779E+07, 1.9253E-11, 1.2242E+04, 7.7277E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.9332E+04,
       true      , false     , 9.7418E-01, 1.0686E+00, latvap, 1.8000E+03, 2.9002E+02, 4.0751E-03, 5.3194E-05,-6.7096E+09,-5.3194E-05,
       -8.4862E+07},

      {3.2784E-08, 1.8780E+07, 2.1753E-11, 1.2461E+04, 7.8657E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.8748E+04,
       true      , false     , 9.8387E-01, 1.0741E+00, latvap, 1.8000E+03, 2.9033E+02, 3.7211E-03, 5.9050E-05,-6.6723E+09,-5.9050E-05,
       -8.6159E+07},

      {3.6045E-08, 1.8780E+07, 2.4356E-11, 1.2689E+04, 8.0098E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.8154E+04,
       true      , false     , 9.9391E-01, 1.0796E+00, latvap, 1.8000E+03, 2.9063E+02, 3.3756E-03, 6.4925E-05,-6.6333E+09,-6.4925E-05,
       -8.7530E+07},

      {3.9321E-08, 1.8780E+07, 2.7069E-11, 1.2928E+04, 8.1605E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.7552E+04,
       true      , false     , 1.0043E+00, 1.0853E+00, latvap, 1.8000E+03, 2.9092E+02, 3.0417E-03, 7.0827E-05,-6.5924E+09,-7.0827E-05,
       -8.8982E+07},

      {4.2614E-08, 1.8780E+07, 2.9903E-11, 1.3178E+04, 8.3182E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.6939E+04,
       true      , false      , 1.0151E+00, 1.0911E+00, latvap, 1.8000E+03, 2.9119E+02, 2.7224E-03, 7.6760E-05,-6.5494E+09,-7.6760E-05,
       -9.0523E+07},

      {4.5927E-08, 1.8780E+07, 3.2867E-11, 1.3440E+04, 8.4833E+03, 0.0000E+00, 0.0000E+00, 0.0000E+00, 5.6317E+04,
       true      , false      , 1.0263E+00, 1.0970E+00, latvap, 1.8000E+03, 2.9143E+02, 2.4202E-03, 8.2728E-05,-6.5044E+09,-8.2728E-05,
       -9.0778E+07},
    };

    // Sync to device
    view_1d<P3UpdatePrognosticLiqData> pupldc_device("pupldc", max_pack_size);
    auto pupldc_host = Kokkos::create_mirror_view(pupldc_device);

    // This copy only copies the input variables.
    std::copy(&pupldc[0], &pupldc[0] + max_pack_size, pupldc_host.data());
    Kokkos::deep_copy(pupldc_device, pupldc_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        pupldc[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend, inv_rho,
        inv_exner, th_atm, qv, qc, nc, qr, nr;
      bool do_predict_nc, do_prescribed_CCN;
      Scalar dt;

      // variables with single values assigned outside of the for loop
      dt                = pupldc_device(0).dt;
      do_predict_nc     = pupldc_device(0).do_predict_nc;
      do_prescribed_CCN = pupldc_device(0).do_prescribed_CCN;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qc2qr_accret_tend[s]     = pupldc_device(vs).qc2qr_accret_tend;
        nc_accret_tend[s]        = pupldc_device(vs).nc_accret_tend;
        qc2qr_autoconv_tend[s]   = pupldc_device(vs).qc2qr_autoconv_tend;
        nc2nr_autoconv_tend[s]   = pupldc_device(vs).nc2nr_autoconv_tend;
        ncautr[s]                = pupldc_device(vs).ncautr;
        nc_selfcollect_tend[s]   = pupldc_device(vs).nc_selfcollect_tend;
        qr2qv_evap_tend[s]       = pupldc_device(vs).qr2qv_evap_tend;
        nr_evap_tend[s]          = pupldc_device(vs).nr_evap_tend;
        nr_selfcollect_tend[s]   = pupldc_device(vs).nr_selfcollect_tend;
        inv_rho[s]               = pupldc_device(vs).inv_rho;
        inv_exner[s]                 = pupldc_device(vs).inv_exner;

        th_atm[s]  = pupldc_device(vs).th_atm;
        qv[s]      = pupldc_device(vs).qv;
        qc[s]      = pupldc_device(vs).qc;
        nc[s]      = pupldc_device(vs).nc;
        qr[s]      = pupldc_device(vs).qr;
        nr[s]      = pupldc_device(vs).nr;
      }

      Functions::update_prognostic_liquid(qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend,
                                          qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend, do_predict_nc, do_prescribed_CCN, inv_rho, inv_exner,
                                          dt, th_atm, qv, qc, nc, qr, nr);

      // Copy results back into views
      pupldc_device(0).dt            = dt;
      pupldc_device(0).do_predict_nc = do_predict_nc;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        pupldc_device(vs).qc2qr_accret_tend     = qc2qr_accret_tend[s];
        pupldc_device(vs).nc_accret_tend        = nc_accret_tend[s];
        pupldc_device(vs).qc2qr_autoconv_tend   = qc2qr_autoconv_tend[s];
        pupldc_device(vs).nc2nr_autoconv_tend   = nc2nr_autoconv_tend[s];
        pupldc_device(vs).ncautr                = ncautr[s];
        pupldc_device(vs).nc_selfcollect_tend   = nc_selfcollect_tend[s];
        pupldc_device(vs).qr2qv_evap_tend       = qr2qv_evap_tend[s];
        pupldc_device(vs).nr_evap_tend          = nr_evap_tend[s];
        pupldc_device(vs).nr_selfcollect_tend   = nr_selfcollect_tend[s];
        pupldc_device(vs).inv_rho               = inv_rho[s];
        pupldc_device(vs).inv_exner                 = inv_exner[s];

        pupldc_device(vs).th_atm  = th_atm[s];
        pupldc_device(vs).qv      = qv[s];
        pupldc_device(vs).qc      = qc[s];
        pupldc_device(vs).nc      = nc[s];
        pupldc_device(vs).qr      = qr[s];
        pupldc_device(vs).nr      = nr[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(pupldc_host, pupldc_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(pupldc[s].th_atm == pupldc_host(s).th_atm);
        REQUIRE(pupldc[s].qv     == pupldc_host(s).qv);
        REQUIRE(pupldc[s].qc     == pupldc_host(s).qc);
        REQUIRE(pupldc[s].nc     == pupldc_host(s).nc);
        REQUIRE(pupldc[s].qr     == pupldc_host(s).qr);
        REQUIRE(pupldc[s].nr     == pupldc_host(s).nr);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        pupldc_host(s).write(Base::m_fid);
      }
    }
  }

  void run_bfb() {
    update_prognostic_liquid_unit_bfb_tests();
  }

}; //TestP3UpdatePrognosticLiq

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3FunctionsImposeMaxTotalNi : public UnitWrap::UnitTest<D>::Base
{
  void impose_max_total_ni_bfb_test() {
    constexpr Scalar max_total_ni = 740.0e3;

    ImposeMaxTotalNiData dc[max_pack_size]= {
      // ni_local, max_total_ni, inv_rho_local
      {0.000E0, max_total_ni, 5.466E3},
      {3.358E4, max_total_ni, 9.691E-1},
      {0.000E0, max_total_ni, 9.105E-1},
      {0.000E3, max_total_ni, 3.371E0},

      {0.000E0, max_total_ni, 5.466E3},
      {3.358E4, max_total_ni, 9.691E-1},
      {0.000E0, max_total_ni, 9.105E-1},
      {0.000E3, max_total_ni, 3.371E0},

      {0.000E0, max_total_ni, 5.466E3},
      {3.358E4, max_total_ni, 9.691E-1},
      {0.000E0, max_total_ni, 9.105E-1},
      {0.000E3, max_total_ni, 3.371E0},

      {0.000E0, max_total_ni, 5.466E3},
      {3.358E4, max_total_ni, 9.691E-1},
      {0.000E0, max_total_ni, 9.105E-1},
      {0.000E3, max_total_ni, 3.371E0},
    };

    //Sync to device
    view_1d<ImposeMaxTotalNiData> dc_device("dc", max_pack_size);
    auto dc_host = Kokkos::create_mirror_view(dc_device);

    //This copy only copies the input variables.
    std::copy(&dc[0], &dc[0] + max_pack_size, dc_host.data());
    Kokkos::deep_copy(dc_device, dc_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        dc[i].read(Base::m_fid);
      }
    }

    //Run function from a kernal and copy results back to the host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack ni_local, inv_rho_local;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        ni_local[s]   = dc_device(vs).ni_local;
        inv_rho_local[s] = dc_device(vs).inv_rho_local;
      }

      Functions::impose_max_total_ni(ni_local, dc_device(0).max_total_ni, inv_rho_local);
      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        dc_device(vs).ni_local   = ni_local[s];
        dc_device(vs).inv_rho_local = inv_rho_local[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(dc_host, dc_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(dc[s].ni_local   == dc_host(s).ni_local);
        REQUIRE(dc[s].inv_rho_local == dc_host(s).inv_rho_local);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        dc_host(s).write(Base::m_fid);
      }
    }
  }

  void run_bfb() {
    impose_max_total_ni_bfb_test();
  }

}; // TestP3FunctionsImposeMaxTotalNi

}//namespace unit_test
}//namespace p3
}//namespace scream

namespace {

TEST_CASE("p3_conservation_test", "[p3_unit_tests]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Conservation;

  T t;
  t.run();
  t.run_bfb();
}

TEST_CASE("p3_get_time_space_phys_variables_test", "[p3_unit_tests]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGetTimeSpacePhysVariables;

  T t;
  t.run_bfb();
}

TEST_CASE("p3_update_prognostic_ice_test", "[p3_unit_tests]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3UpdatePrognosticIce;

  T t;
  t.run_bfb();
}

TEST_CASE("p3_update_prognostic_liquid_test", "[p3_unit_tests]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3UpdatePrognosticLiq;

  T t;
  t.run_bfb();
}

TEST_CASE("p3_impose_max_total_ni_test", "[p3_unit_tests]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3FunctionsImposeMaxTotalNi;

  T t;
  t.run_bfb();
}

} // namespace
