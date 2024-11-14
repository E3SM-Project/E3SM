#include "catch2/catch.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"
#include "share/scream_types.hpp"
#include "physics/share/physics_functions.hpp"

#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPreventLiqSupersaturation : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  //Conceptual tests for prevent_liq_supersaturation. Note many conceptual tests make sense to run on
  //random data, so are included in run_bfb rather than here.
  {

    using physics = scream::physics::Functions<Scalar, Device>;

    constexpr Scalar inv_cp       = C::INV_CP;
    constexpr Scalar latvap       = C::LatVap;
    constexpr Scalar latice       = C::LatIce;

    //Start with reasonable values
    //============================
    Spack pres(100000);
    Spack t_atm(270);
    Spack qv(1.7e-3);
    Scalar dt=60;
    Spack qidep(0);
    Spack qinuc(0); //dep and nuc only used together, so only fiddle with one
    //Note that inout variables are initialized for each test to avoid using unintended vals
    auto context= qinuc == 0; //want a pack-sized array which is all true.

    //Test that large initial endstep supersat is removed.
    //=============================
    //First, make copies of output variables and make them big
    Spack qi2qv_sublim_tend_tmp(1e-5);
    Spack qr2qv_evap_tend_tmp(1e-5);

    //Next, confirm that initial data is subsaturated at beginning of step and supersaturated at end
    //(using code copied from prevent_liq_supersaturation_impl.hpp)
    Spack qv_sinks=qidep + qinuc;
    Spack qv_sources=qi2qv_sublim_tend_tmp + qr2qv_evap_tend_tmp;
    Spack qv_endstep=qv - qv_sinks*dt + qv_sources*dt;
    Spack T_endstep=t_atm + ( (qv_sinks-qi2qv_sublim_tend_tmp)*(latvap+latice)*inv_cp
			- qr2qv_evap_tend_tmp*latvap*inv_cp )*dt;
    Spack qsl = physics::qv_sat_dry(T_endstep,pres,false,context); //"false" means NOT sat w/ respect to ice
    //just require index 0 since all entries are identical

    REQUIRE(        qv[0]<qsl[0]); //not a test of prevent_liq_supersat, just a sanity check that
    REQUIRE(qv_endstep[0]>qsl[0]); // inputs to this test are testing what we want.

    //Now update sublim and evap tends using prevent_liq_supersaturation
    Functions::prevent_liq_supersaturation(pres, t_atm, qv,
					   dt, qidep, qinuc, qi2qv_sublim_tend_tmp, qr2qv_evap_tend_tmp);

    //Finally, recompute liquid saturation
    //(using code copied from prevent_liq_supersaturation_impl.hpp)
    qv_sinks=qidep + qinuc;
    qv_sources=qi2qv_sublim_tend_tmp + qr2qv_evap_tend_tmp;
    qv_endstep=qv - qv_sinks*dt + qv_sources*dt;
    T_endstep=t_atm + ( (qv_sinks-qi2qv_sublim_tend_tmp)*(latvap+latice)*inv_cp
			- qr2qv_evap_tend_tmp*latvap*inv_cp )*dt;
    qsl = physics::qv_sat_dry(T_endstep,pres,false,context); //"false" means NOT sat w/ respect to ice
    //just require index 0 since all entries are identical

    REQUIRE(qv_endstep[0]<=qsl[0]);

    //If qv is supersaturated by itself, evap and sublim should be 0:
    //============================
    Spack qv_tmp(1e-2);
    Spack qi2qv_sublim_tend_tmp2(1e-4);
    Spack qr2qv_evap_tend_tmp2(1e-4);

    Functions::prevent_liq_supersaturation(pres, t_atm, qv_tmp,
					   dt, qidep, qinuc, qi2qv_sublim_tend_tmp2, qr2qv_evap_tend_tmp2);
    //just require index 0 since all entries are identical.
    REQUIRE( qi2qv_sublim_tend_tmp2[0] ==0 );
    REQUIRE( qr2qv_evap_tend_tmp2[0] == 0 );

  } //end run_property

  void run_bfb()
  {
    constexpr Scalar latvap = C::LatVap;
    constexpr Scalar latice = C::LatIce;

    auto engine = Base::get_engine();

    PreventLiqSupersaturationData baseline_data[max_pack_size];

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
      d.dt = baseline_data[0].dt; // Hold this fixed, this is not packed data

      // C++ impl uses constants for latent_heat values. Manually set here
      // so BASELINE can match
      d.latent_heat_vapor = latvap;
      d.latent_heat_sublim = latvap+latice;
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before
    // reads so that inout data is in original state
    view_1d<PreventLiqSupersaturationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&baseline_data[0], &baseline_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Save copy of inout vars to check that prevent_liq_supersaturation always makes them smaller
    Real qi2qv_sublim_tend_init[max_pack_size],qr2qv_evap_tend_init[max_pack_size];
    for (Int i = 0; i < max_pack_size; ++i) {
      qi2qv_sublim_tend_init[i] = cxx_host(i).qi2qv_sublim_tend;
      qr2qv_evap_tend_init[i] = cxx_host(i).qr2qv_evap_tend;
    }

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        baseline_data[i].read(Base::m_fid);
      }
    }

    // Get data from cxx. Run prevent_liq_supersaturation from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Scalar dt;
      Spack pres, qi2qv_sublim_tend, qidep, qinuc, qr2qv_evap_tend, qv, t_atm;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
	dt = cxx_device(vs).dt; //dt is scalar but PreventLiqSupersaturationData has diff val for each row.
        pres[s] = cxx_device(vs).pres;
        qi2qv_sublim_tend[s] = cxx_device(vs).qi2qv_sublim_tend;
        qidep[s] = cxx_device(vs).qidep;
        qinuc[s] = cxx_device(vs).qinuc;
        qr2qv_evap_tend[s] = cxx_device(vs).qr2qv_evap_tend;
        qv[s] = cxx_device(vs).qv;
        t_atm[s] = cxx_device(vs).t_atm;
      }

      Functions::prevent_liq_supersaturation(pres, t_atm, qv, dt, qidep, qinuc, qi2qv_sublim_tend, qr2qv_evap_tend);

      // Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_device(vs).qi2qv_sublim_tend = qi2qv_sublim_tend[s];
        cxx_device(vs).qr2qv_evap_tend = qr2qv_evap_tend[s];
      }

    });

    Kokkos::deep_copy(cxx_host, cxx_device);

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        // Verify BFB results
        PreventLiqSupersaturationData& d_baseline = baseline_data[i];
        PreventLiqSupersaturationData& d_cxx = cxx_host[i];
        REQUIRE(d_baseline.qi2qv_sublim_tend == d_cxx.qi2qv_sublim_tend);
        REQUIRE(d_baseline.qr2qv_evap_tend == d_cxx.qr2qv_evap_tend);

        //Verify tendencies are always >=0:
        REQUIRE(d_cxx.qi2qv_sublim_tend>=0);
        REQUIRE(d_cxx.qr2qv_evap_tend>=0);

        //Verify function call always makes tendencies smaller
        REQUIRE(d_cxx.qi2qv_sublim_tend<=qi2qv_sublim_tend_init[i]);
        REQUIRE(d_cxx.qr2qv_evap_tend<=qr2qv_evap_tend_init[i]);
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

TEST_CASE("prevent_liq_supersaturation_property", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPreventLiqSupersaturation;

  T t;
  t.run_property();
}

TEST_CASE("prevent_liq_supersaturation_bfb", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPreventLiqSupersaturation;

  T t;
  t.run_bfb();
}

} // empty namespace
