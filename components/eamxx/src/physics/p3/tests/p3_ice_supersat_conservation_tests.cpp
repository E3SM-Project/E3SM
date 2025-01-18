#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestIceSupersatConservation : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    constexpr Scalar latvap       = C::LatVap;
    constexpr Scalar latice       = C::LatIce;

    auto engine = Base::get_engine();

    IceSupersatConservationData baseline_data[max_pack_size];

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
      d.dt = baseline_data[0].dt; // hold this fixed, it is not packed data

      // C++ impl uses constants for latent_heat values. Manually set here
      // so F90 can match
      d.latent_heat_sublim = latvap+latice;
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before reads so that
    // inout data is in original state
    view_1d<IceSupersatConservationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&baseline_data[0], &baseline_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        baseline_data[i].read(Base::m_fid);
      }
    }

    // Get data from cxx. Run ice_supersat_conservation from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack cld_frac_i, qidep, qinuc, qv, qv_sat_i, t_atm, qi2qv_sublim_tend, qr2qv_evap_tend;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cld_frac_i[s] = cxx_device(vs).cld_frac_i;
        qidep[s] = cxx_device(vs).qidep;
        qinuc[s] = cxx_device(vs).qinuc;
        qv[s] = cxx_device(vs).qv;
        qv_sat_i[s] = cxx_device(vs).qv_sat_i;
        t_atm[s] = cxx_device(vs).t_atm;
        qi2qv_sublim_tend[s] = cxx_device(vs).qi2qv_sublim_tend;
        qr2qv_evap_tend[s] = cxx_device(vs).qr2qv_evap_tend;
      }

      Functions::ice_supersat_conservation(qidep, qinuc, cld_frac_i, qv, qv_sat_i, t_atm, cxx_device(offset).dt, qi2qv_sublim_tend, qr2qv_evap_tend);

      // Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_device(vs).qidep = qidep[s];
        cxx_device(vs).qinuc = qinuc[s];
      }
    });

    Kokkos::deep_copy(cxx_host, cxx_device);

    // Verify BFB results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        IceSupersatConservationData& d_f90 = baseline_data[i];
        IceSupersatConservationData& d_cxx = cxx_host[i];
        REQUIRE(d_f90.qidep == d_cxx.qidep);
        REQUIRE(d_f90.qinuc == d_cxx.qinuc);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        cxx_host(s).write(Base::m_fid);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("ice_supersat_conservation_bfb", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceSupersatConservation;

  T t;
  t.run_bfb();
}

} // empty namespace
