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
struct UnitWrap::UnitTest<D>::TestEvapSublPrecip : public UnitWrap::UnitTest<D>::Base {

  void run_property() {

    //TEST WEIGHTING TIMESCALE
    //========================
    //dt/tau ~ 0 => weight => 1. A value of exactly 0 would cause div by 0 though.
    Spack wt;
    Functions::rain_evap_tscale_weight(Spack(1e-8),wt);

    REQUIRE( wt[0] >= 0 ); //always true
    REQUIRE( wt[0] <= 1 ); //always true
    REQUIRE( 1-wt[0] < 1e-8 );

    //dt/tau->inf => weight => 0.
    Functions::rain_evap_tscale_weight(Spack(1e15),wt);
    REQUIRE( wt[0] >= 0 ); //always true
    REQUIRE( wt[0] <= 1 ); //always true
    REQUIRE( wt[0] < 1e-8 );

    //TEST EQUILIB EVAP RATE:
    //=======================
    //if A_c=0, equilibrium evap rate should be zero
    Spack tend;
    Functions::rain_evap_equilib_tend(Spack(0),Spack(1),Spack(1),Spack(1),tend);
    REQUIRE( std::abs(tend[0])<1e-8 );

    //if tau_eff==tau_r, equilibrium should be -A_c/ab
    Spack A_c(2);
    Spack ab(1.2);
    Functions::rain_evap_equilib_tend(A_c,ab,Spack(1),Spack(1),tend);
    REQUIRE( std::abs(tend[0] + A_c[0]/ab[0]) < 1e-8 );

    //TEST INSTANTANEOUS EVAP RATE:
    //=============================
    //when ssat_r=0, tend should be zero.
    Functions::rain_evap_instant_tend(Spack(0),ab,Spack(1),tend);
    REQUIRE( tend[0] == 0 );

    //when tau_r=>inf, tend should approach 0.
    Functions::rain_evap_instant_tend(Spack(-1e-3),ab,Spack(1e12),tend);
    REQUIRE( tend[0] > 0 ); //always true for ssat_r<0.
    REQUIRE( tend[0] < 1e-8 );

    //TEST ACTUAL EVAP FUNCTION:
    //==========================
    //first, establish reasonable values to use as default:
    Spack qr_incld(1e-4);
    Spack qc_incld(1e-4);
    Spack nr_incld(1e7);
    Spack qi_incld(1e-4);
    Spack cld_frac_l(0.1);
    Spack cld_frac_r(0.9);
    Spack qv(5e-4);
    Spack qv_prev(1e-3);
    Spack qv_sat_l(1e-3); //must be > qv for subsat
    Spack qv_sat_i(6e-4); //making slightly higher than qv
    Spack abi(1.2);
    Spack epsr(1./20.);
    Spack epsi_tot(1./60.);
    Spack t(287);
    Spack t_prev(285);
    Spack dqsdt(1e-3);
    Scalar dt=60;
    Spack qrtend;
    Spack nrtend;

    //if qr_incld is too small, evap rate should be zero.
    constexpr Scalar QSMALL   = C::QSMALL;
    Functions::evaporate_rain(Spack(QSMALL/2),qc_incld,nr_incld,qi_incld, //qr_incld->QSMALL/2
			      cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i,
			      ab,abi,epsr,epsi_tot,t,t_prev,dqsdt,dt,
			      qrtend,nrtend);
    REQUIRE( std::abs(qrtend[0])<1e-8 );
    REQUIRE( std::abs(nrtend[0])<1e-8 );

    //if qr_incld is small enough but not too small, evap rate should be qr_incld/dt.
    Spack qr_tiny=Spack(5e-13);
    Functions::evaporate_rain(qr_tiny,qc_incld,nr_incld,qi_incld, //qr_incld->_tiny
				cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i,
				ab,abi,epsr,epsi_tot,t,t_prev,dqsdt,dt,
				qrtend,nrtend);
      REQUIRE( std::abs(qrtend[0] - qr_tiny[0]/dt
			*(cld_frac_r[0]-cld_frac_l[0])/cld_frac_r[0])<1e-8 );
      REQUIRE( std::abs(nrtend[0] - qrtend[0]*nr_incld[0]/qr_tiny[0]) < 1e-8 );//always true
      REQUIRE( nrtend[0] <= nr_incld[0]/dt); //keep end-of-step nr positive. Should always be true.

    //if no rainy areas outside cloud, don't evap
    Functions::evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld,
			      cld_frac_r,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i, //cld_frac_l->_r
			      ab,abi,epsr,epsi_tot,t,t_prev,dqsdt,dt,
			      qrtend,nrtend);
    REQUIRE( std::abs(qrtend[0])<1e-8 );
    REQUIRE( std::abs(nrtend[0])<1e-8 );

    //no evap if supersaturated
    Functions::evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld,
			      //set qv->qv_sat_l*2 in next line to ensure supersaturated.
			      cld_frac_l,cld_frac_r,qv_sat_l*2,qv_prev,qv_sat_l,qv_sat_i,
			      ab,abi,epsr,epsi_tot,t,t_prev,dqsdt,dt,
			      qrtend,nrtend);
    REQUIRE( std::abs(qrtend[0])<1e-8 );
    REQUIRE( std::abs(nrtend[0])<1e-8 );

    //for case with lots of evap, make sure doesn't overdeplete qr_incld
    Functions::evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld,
			      //qv -> qv*0.1 to encourage lots of rain evap
			      cld_frac_l,cld_frac_r,qv*0.1,qv_prev,qv_sat_l,qv_sat_i,
			      ab,abi,epsr,epsi_tot,t,t_prev,dqsdt,dt,
			      qrtend,nrtend);
    REQUIRE( qrtend[0] <= qr_incld[0]/dt);
    REQUIRE( nrtend[0] <= nr_incld[0]/dt); //keep end-of-step nr positive. Should always be true.

  } //end run_property

  void run_bfb() {
    constexpr Scalar latvap = C::LatVap;
    constexpr Scalar latice = C::LatIce;

    //baseline generated data is input to the following
    //This subroutine has 20 args, only 18 are supplied here for invoking it as last 2 are intent-outs
    //note that dt is the same val for each row - this is needed since dt is a scalar and all rows are executed simultaneously on CPU in C++.
    //row1: above freezing, should trigger
    //row2: below freezing, should trigger but qr_incld small enough to evap all
    //row3: supersaturated, shouldn't trigger
    //row4: below freezing, should trigger
    //row5: below freezing, should trigger
    //row6: cld frac > rain area frac. shouldn't trigger.
    //row7: qr_incld=0. shouldn't trigger.
    //row8: above freezing, should trigger.
    //rows 9-16: random junk but ensured cld_frac_r>cld_frac_l and subsaturated.
    EvapRainData espd[max_pack_size] = {
    //qr_incld,     qc_incld,    nr_incld,    qi_incld,    cld_frac_l,  cld_frac_r,  qv,          qv_prev,     qv_sat_l,    qv_sat_i,    ab,          abi,         epsr,        epsi_tot,    t,           t_prev,      lat_ht_sublim, dqsdt,     dt
      {4.634940e-03,1.215335e-03,6.073270e+07,3.594486e-04,6.134229e-01,9.134229e-01,2.747871e-03,1.911238e-03,5.913313e-03,1.057645e-03,1.782748e+00,1.571392e+00,3.868229e+02,2.248689e+02,3.101180e+02,1.395063e+02,latvap+latice,5.494606e-03,6.000000e+02},
      {6.175320e-13,4.432407e-03,8.029967e+07,1.905151e-03,2.190099e-01,7.031070e-01,4.172977e-05,7.315360e-03,7.280063e-03,1.378543e-03,1.461443e+00,1.507382e+00,8.452377e+02,1.971876e+02,2.389249e+02,1.497752e+02,latvap+latice,5.107905e-03,6.000000e+02},
      {4.519798e-03,7.348916e-03,7.420725e+07,2.220971e-03,1.882608e-01,2.934182e-01,4.957590e-03,2.550256e-03,3.136926e-03,4.498115e-03,1.433526e+00,1.207516e+00,9.716844e+02,5.602546e+01,1.389465e+02,1.075863e+02,latvap+latice,6.771428e-03,6.000000e+02},
      {7.169182e-03,6.657331e-03,9.807967e+07,7.981196e-03,2.914473e-01,6.375719e-01,2.420032e-03,1.223012e-03,7.685516e-03,5.207024e-03,1.644865e+00,1.433872e+00,3.825069e+02,6.550300e+02,1.833466e+02,1.741918e+02,latvap+latice,3.792982e-03,6.000000e+02},
      {1.103118e-03,9.158125e-03,3.136196e+07,4.286154e-03,2.699078e-01,4.668103e-01,9.645460e-03,6.379119e-03,8.283285e-03,3.342400e-03,1.546698e+00,1.417916e+00,9.289270e+02,9.844129e+02,2.543202e+02,1.932996e+02,latvap+latice,2.693119e-03,6.000000e+02},
      {4.308000e-03,8.168535e-03,7.439969e+07,5.131497e-03,6.851225e-01,3.298025e-01,4.331812e-03,2.814373e-03,3.592807e-03,1.527499e-03,1.856943e+00,1.003269e+00,9.165690e+02,9.379921e+02,2.163204e+02,3.165814e+02,latvap+latice,6.801393e-03,6.000000e+02},
      {0.000000e-00,3.318968e-03,4.664041e+07,8.737282e-03,2.585907e-01,6.297295e-02,8.747418e-03,2.710437e-03,2.164895e-03,9.455725e-03,1.241506e+00,1.561393e+00,2.492674e+02,6.546182e+02,2.228772e+02,2.147968e+02,latvap+latice,5.903261e-03,6.000000e+02},
      {7.677170e-03,6.069057e-05,6.404241e+07,3.094233e-03,3.755403e-01,5.026876e-01,4.723817e-03,1.204228e-03,6.156526e-03,8.194797e-03,1.361509e+00,1.772751e+00,6.420537e+01,4.043364e+02,2.833110e+02,3.314521e+02,latvap+latice,2.996696e-03,6.000000e+02},

      {9.999294e-03,3.138400e-03,2.355097e+07,9.897893e-03,7.667177e-01,9.739270e-01,4.221430e-03,3.570130e-03,8.370033e-03,9.527208e-03,1.597218e+00,1.111438e+00,7.832357e+02,8.364566e+02,2.854867e+02,2.340771e+02,latvap+latice,7.235757e-03,6.000000e+02},
      {8.841793e-03,3.530456e-03,9.618284e+07,9.311658e-03,3.458590e-01,6.978258e-01,1.279864e-03,4.652008e-03,1.869728e-03,8.931663e-03,1.712564e+00,1.223882e+00,9.692403e+02,2.358558e+02,3.204043e+02,1.827677e+02,latvap+latice,7.646405e-03,6.000000e+02},
      {1.425612e-03,6.653411e-04,2.843806e+07,1.922560e-03,9.100262e-01,0.996264e-01,8.973183e-04,9.857420e-03,6.221419e-03,8.133433e-03,1.815337e+00,1.885506e+00,5.508742e+02,1.612139e+02,2.798523e+02,2.631136e+02,latvap+latice,4.148666e-03,6.000000e+02},
      {4.125177e-04,4.056163e-03,2.716439e+07,6.484214e-03,1.658752e-01,2.859102e-01,5.724081e-03,6.282997e-03,7.313187e-03,6.049825e-03,1.140910e+00,1.145941e+00,7.490652e+02,5.011633e+02,1.986541e+02,2.745566e+02,latvap+latice,6.784784e-03,6.000000e+02},
      {5.010628e-03,2.863789e-04,8.953841e+07,3.953058e-03,1.135952e-01,9.718675e-01,1.846157e-03,5.743094e-03,2.842649e-03,8.155366e-03,1.227867e+00,1.894249e+00,1.161776e+02,3.578576e+02,1.240083e+02,1.639791e+02,latvap+latice,4.497257e-03,6.000000e+02},
      {9.487866e-03,6.584660e-03,6.149682e+06,9.413342e-03,4.757261e-01,6.503885e-01,1.078922e-03,3.489665e-03,3.059596e-03,9.285703e-03,1.192620e+00,1.967205e+00,5.085628e+02,3.741816e+01,1.196252e+02,2.904002e+02,latvap+latice,2.566077e-03,6.000000e+02},
      {3.241928e-03,7.024929e-03,2.212493e+07,8.600485e-03,3.963690e-01,4.834201e-01,3.736511e-03,5.724475e-03,4.790239e-03,2.766218e-03,1.151150e+00,1.150516e+00,2.089426e+02,8.666450e+02,1.898220e+02,2.862496e+02,latvap+latice,7.039800e-03,6.000000e+02},
      {4.617594e-03,3.157739e-03,5.569465e+07,8.221076e-03,7.918279e-01,9.995014e-01,1.338309e-04,1.319707e-03,2.896082e-03,4.359171e-03,1.007827e+00,1.812954e+00,5.332209e+02,2.973599e+02,3.271466e+02,2.622351e+02,latvap+latice,1.407429e-03,6.000000e+02}
    };

    // Sync to device
    view_1d<EvapRainData> espd_device("espd", max_pack_size);
    auto espd_host = Kokkos::create_mirror_view(espd_device);

    // This copy only copies the input variables.
    std::copy(&espd[0], &espd[0] + max_pack_size, espd_host.data());
    Kokkos::deep_copy(espd_device, espd_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        espd[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qr_incld,qc_incld,nr_incld,qi_incld,
	cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i,
	ab,abi,epsr,epsi_tot,t,t_prev,dqsdt;

      Scalar dt;

      // Init pack outputs
      Spack qr2qv_evap_tend, nr_evap_tend;

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qr_incld[s]    = espd_device(vs).qr_incld;
        qc_incld[s]    = espd_device(vs).qc_incld;
        nr_incld[s]    = espd_device(vs).nr_incld;
        qi_incld[s] = espd_device(vs).qi_incld;
        cld_frac_l[s]       = espd_device(vs).cld_frac_l;
        cld_frac_r[s]       = espd_device(vs).cld_frac_r;
	qv[s] = espd_device(vs).qv;
	qv_prev[s] = espd_device(vs).qv_prev;
        qv_sat_l[s]         = espd_device(vs).qv_sat_l;
	qv_sat_i[s]    = espd_device(vs).qv_sat_i;
        ab[s]          = espd_device(vs).ab;
	abi[s]         = espd_device(vs).abi;
        epsr[s]        = espd_device(vs).epsr;
	epsi_tot[s]    = espd_device(vs).epsi_tot;
	t[s]           = espd_device(vs).t;
	t_prev[s]      = espd_device(vs).t_prev;
	dqsdt[s]=espd_device(vs).dqsdt;
	dt=espd_device(vs).dt;
        //qr2qv_evap_tend[s]       = espd_device(vs).qr2qv_evap_tend; //PMC shouldn't have to init output vars.
        //nr_evap_tend[s]       = espd_device(vs).nr_evap_tend;
      }

      Functions::evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld,
				cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i,
				ab,abi,epsr,epsi_tot,t,t_prev,dqsdt,dt,
				qr2qv_evap_tend,nr_evap_tend);

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        espd_device(vs).qr_incld    = qr_incld[s];
        espd_device(vs).qc_incld    = qc_incld[s];
        espd_device(vs).nr_incld    = nr_incld[s];
        espd_device(vs).qi_incld = qi_incld[s];
        espd_device(vs).cld_frac_l       = cld_frac_l[s];
        espd_device(vs).cld_frac_r       = cld_frac_r[s];
        espd_device(vs).qv_sat_l         = qv_sat_l[s];
        espd_device(vs).ab          = ab[s];
        espd_device(vs).epsr        = epsr[s];
        espd_device(vs).qv          = qv[s];
        espd_device(vs).qr2qv_evap_tend       = qr2qv_evap_tend[s];
        espd_device(vs).nr_evap_tend       = nr_evap_tend[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(espd_host, espd_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(espd[s].qr2qv_evap_tend == espd_host(s).qr2qv_evap_tend);
        REQUIRE(espd[s].nr_evap_tend == espd_host(s).nr_evap_tend);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        espd_host(s).write(Base::m_fid);
      }
    }
  } // end run_bfb


}; //TestEvapSublPrecip UnitWrap

}//namespace unit_test
}//namespace p3
}//namespace scream

namespace {

  TEST_CASE("p3_evaporate_rain_property", "p3_unit_tests")
  {
    using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestEvapSublPrecip;

    T t;
    t.run_property();
  }

  TEST_CASE("p3_evaporate_rain_test", "p3_unit_tests")
  {
    using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestEvapSublPrecip;

    T t;
    t.run_bfb();
  }

}// end anonymous namespace
