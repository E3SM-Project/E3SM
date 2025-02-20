#include "catch2/catch.hpp"

#include "physics/share/physics_functions.hpp"
#include "physics/share/physics_saturation_impl.hpp"
#include "physics_unit_tests_common.hpp"

#include "share/eamxx_types.hpp"
#include "share/eamxx_session.hpp"
#include "share/util/eamxx_utils.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_file_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iomanip>      // std::setprecision

namespace scream {
namespace physics {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestSaturation
{
  KOKKOS_FUNCTION  static void saturation_tests(
    const Scalar& temperature, const Scalar& pressure,
    Scalar& sat_ice_fp, Scalar& sat_liq_fp, Scalar& mix_ice_fr, Scalar& mix_liq_fr,
    Scalar& sat_ice_mkp, Scalar& sat_liq_mkp, Scalar& mix_ice_mkr, Scalar& mix_liq_mkr)
  {
    //Nomenclature:
    //subscript "_fp"  stands for "Flatau Pressure"
    //subscript "_fr"  stands for "Flatau mixing Ratios"
    //subscript "_mkp" stands for "Murphy Koop Pressure"
    //subscript "_mkr" stands for "Murphy Koop mixing Ratios"

    //Allow usage of saturation functions
    using physics = scream::physics::Functions<Scalar, Device>;

    //Convert Scalar inputs to Spacks because that's what polysvp1 and qv_sat expect as inputs.
    //--------------------------------------
    const Spack temps(temperature);
    const Spack pres(pressure);

    //Get values from polysvp1 and qv_sat (qv_sat calls polysvp1 here) to test against "expected" values
    //--------------------------------------
    sat_ice_fp  = physics::polysvp1(temps, true, Smask(true))[0];
    sat_liq_fp  = physics::polysvp1(temps, false, Smask(true))[0];

    //Functions<S,D>::qv_sat_dry(const Spack& t_atm, const Spack& p_atm_dry, const bool ice, const Smask& range_mask,
    //                           const SaturationFcn func_idx, const char* caller)
    mix_ice_fr = physics::qv_sat_dry(temps, pres, true, Smask(true), physics::Polysvp1)[0];
    mix_liq_fr = physics::qv_sat_dry(temps, pres, false,Smask(true), physics::Polysvp1)[0];

    //Get values from MurphyKoop_svp and qv_sat_dry (qv_sat_dry calls MurphyKoop_svp here) to test against "expected" values
    sat_ice_mkp   = physics::MurphyKoop_svp(temps, true, Smask(true))[0];
    sat_liq_mkp   = physics::MurphyKoop_svp(temps, false, Smask(true))[0];

    mix_ice_mkr  = physics::qv_sat_dry(temps, pres, true,  Smask(true), physics::MurphyKoop)[0];
    mix_liq_mkr  = physics::qv_sat_dry(temps, pres, false, Smask(true), physics::MurphyKoop)[0];
  }

  static constexpr auto atm_pres = 1e5;
  static constexpr auto tmelt = C::Tmelt;

  /*
    Originally written by Kyle Pressel, updated by Peter Caldwell on 4/5/20 and
    by Jim Foucar on 4/21/23.

    This code tests polysvp1 and qv_sat at 0 degrees C, at a very cold T, and at a very hot T
    to make sure our impl gets the same answer as Flatau et al 1992:
    (https://journals.ametsoc.org/jamc/article/31/12/1507/14870/Polynomial-Fits-to-Saturation-Vapor-Pressure)
    For 0 degrees C, polysvp values can be read directly from Flatau. For other cases, I independently
    coded up the Flatau scheme (polysvp1) in python and used it to derive the expected values. My python code is
    in https://github.com/E3SM-Project/scream-docs.git analysis-scripts/test_qv_sat.py
   */


  TestSaturation () :
    // Hardcode these for now
    params_({
        // Just Freezing Case: Test values @ 273.15K (tmelt) @ 1e5 Pa
        //---------------------------------------
        // This is the nicest test b/c polysvp1 is a polynomial fit around (T-273.15)
        // so T=273.15 collapses back to the intercept coefficient which can be read
        // directly from the RHS of table 4 of Flatau et al 1992.
        {tmelt, atm_pres},
        {243.15, atm_pres},
        {303.15, atm_pres},

        //Following values are picked from Murphy and Koop (2005)
        //Table C1 titled: "VALUES RECOMMENDED FOR CHECKING COMPUTER CODES"
        //Saturation vapor pressure (SVP) values in the table were upto only 5 significant digits.
        {150, atm_pres},
        {180, atm_pres},
        {210, atm_pres},
        {240, atm_pres},
        {273.16, atm_pres},
        {300, atm_pres},
        {243.15, 5e4}
      })
  {}

  Int generate_baseline (const std::string& filename) {
    auto fid = ekat::FILEPtr(fopen(filename.c_str(), "w"));
    EKAT_REQUIRE_MSG( fid, "generate_baseline can't write " << filename);

    for (auto p : params_) {
      ParamSet ps = p;
      Kokkos::View<OutputData*> d_dev("",1);
      Kokkos::parallel_for(1,
        KOKKOS_LAMBDA(const size_t&) {
          TestSaturation::saturation_tests(
            ps.temperature, ps.pressure,
            d_dev[0].sat_ice_fp,
            d_dev[0].sat_liq_fp,
            d_dev[0].mix_ice_fr,
            d_dev[0].mix_liq_fr,
            d_dev[0].sat_ice_mkp,
            d_dev[0].sat_liq_mkp,
            d_dev[0].mix_ice_mkr,
            d_dev[0].mix_liq_mkr);
      });
      Kokkos::fence();

      auto d_host = Kokkos::create_mirror_view(d_dev);
      Kokkos::deep_copy(d_host,d_dev);

      write(fid, d_host[0]); // Save the fields to the baseline file.
    }

    return 0;
  }

  Int run_and_cmp (const std::string& filename, const Scalar& tol) {
    auto fid = ekat::FILEPtr(fopen(filename.c_str(), "r"));
    EKAT_REQUIRE_MSG( fid, "generate_baseline can't read " << filename);
    Int nerr = 0, ne;
    int case_num = 0;
    for (auto p : params_) {
      ++case_num;
      OutputData ref;
      ParamSet ps = p;
      std::cout << "--- checking physics saturation case # " << case_num << std::endl;
      read(fid, ref);

      Kokkos::View<OutputData*> d_dev("",1);
      Kokkos::parallel_for(1,
        KOKKOS_LAMBDA(const size_t&) {
          TestSaturation::saturation_tests(
            ps.temperature, ps.pressure,
            d_dev[0].sat_ice_fp,
            d_dev[0].sat_liq_fp,
            d_dev[0].mix_ice_fr,
            d_dev[0].mix_liq_fr,
            d_dev[0].sat_ice_mkp,
            d_dev[0].sat_liq_mkp,
            d_dev[0].mix_ice_mkr,
            d_dev[0].mix_liq_mkr);
      });
      Kokkos::fence();

      auto d_host = Kokkos::create_mirror_view(d_dev);
      Kokkos::deep_copy(d_host,d_dev);

      ne = compare(tol, ref, d_host[0]);
      if (ne) std::cout << "Ref impl failed.\n";
      nerr += ne;
    }
    return nerr;
  }

#ifndef KOKKOS_ENABLE_CUDA
  // Everything below is private but on CUDA they must be public
 private:
#endif

  // Full specification for a run
  struct ParamSet {
    Scalar temperature;
    Scalar pressure;
  };

  struct OutputData {
    Scalar sat_ice_fp;
    Scalar sat_liq_fp;
    Scalar mix_ice_fr;
    Scalar mix_liq_fr;

    Scalar sat_ice_mkp;
    Scalar sat_liq_mkp;
    Scalar mix_ice_mkr;
    Scalar mix_liq_mkr;

    KOKKOS_INLINE_FUNCTION
    OutputData& operator+=(const OutputData& rhs)
    {
      sat_ice_fp += rhs.sat_ice_fp;
      sat_liq_fp += rhs.sat_liq_fp;
      mix_ice_fr += rhs.mix_ice_fr;
      mix_liq_fr += rhs.mix_liq_fr;

      sat_ice_mkp += rhs.sat_ice_mkp;
      sat_liq_mkp += rhs.sat_liq_mkp;
      mix_ice_mkr += rhs.mix_ice_mkr;
      mix_liq_mkr += rhs.mix_liq_mkr;

      return *this;
    }
  };

  std::vector<ParamSet> params_;

  static void write (const ekat::FILEPtr& fid, const OutputData& d) {
    ekat::write(&d.sat_ice_fp, 1, fid);
    ekat::write(&d.sat_liq_fp, 1, fid);
    ekat::write(&d.mix_ice_fr, 1, fid);
    ekat::write(&d.mix_liq_fr, 1, fid);

    ekat::write(&d.sat_ice_mkp, 1, fid);
    ekat::write(&d.sat_liq_mkp, 1, fid);
    ekat::write(&d.mix_ice_mkr, 1, fid);
    ekat::write(&d.mix_liq_mkr, 1, fid);
  }

  static void read (const ekat::FILEPtr& fid, OutputData& d) {
    ekat::read(&d.sat_ice_fp, 1, fid);
    ekat::read(&d.sat_liq_fp, 1, fid);
    ekat::read(&d.mix_ice_fr, 1, fid);
    ekat::read(&d.mix_liq_fr, 1, fid);

    ekat::read(&d.sat_ice_mkp, 1, fid);
    ekat::read(&d.sat_liq_mkp, 1, fid);
    ekat::read(&d.mix_ice_mkr, 1, fid);
    ekat::read(&d.mix_liq_mkr, 1, fid);
  }

  static Int compare (const Scalar& tol,
                      const OutputData& ref, const OutputData& d) {
    Int nerr = 0;

    nerr += scream::compare("sat_ice_fp", &ref.sat_ice_fp, &d.sat_ice_fp, 1, tol);
    nerr += scream::compare("sat_liq_fp", &ref.sat_liq_fp, &d.sat_liq_fp, 1, tol);
    nerr += scream::compare("mix_ice_fr", &ref.mix_ice_fr, &d.mix_ice_fr, 1, tol);
    nerr += scream::compare("mix_liq_fr", &ref.mix_liq_fr, &d.mix_liq_fr, 1, tol);

    nerr += scream::compare("sat_ice_mkp", &ref.sat_ice_mkp, &d.sat_ice_mkp, 1, tol);
    nerr += scream::compare("sat_liq_mkp", &ref.sat_liq_mkp, &d.sat_liq_mkp, 1, tol);
    nerr += scream::compare("mix_ice_mkr", &ref.mix_ice_mkr, &d.mix_ice_mkr, 1, tol);
    nerr += scream::compare("mix_liq_mkr", &ref.mix_liq_mkr, &d.mix_liq_mkr, 1, tol);

    return nerr;
  }

};

}
}
}

namespace {

void expect_another_arg (int i, int argc) {
  EKAT_REQUIRE_MSG(i != argc-1, "Expected another cmd-line arg.");
}

} // empty namespace

int main (int argc, char** argv) {
  using UnitTest = scream::physics::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>;
  using TestSaturation = UnitTest::TestSaturation;

  int nerr = 0;

  if (argc == 1) {
    std::cout <<
      argv[0] << " [options]\n"
      "Options:\n"
      "  -g                  Generate baseline file. Default False.\n"
      "  -b <baseline-file>  Path to baseline file. Required.\n"
      "  -t <tol>            Tolerance for relative error. Default machine eps (*10000 for Release).\n";
    return 1;
  }

  // Set up options with defaults
  bool generate = false;
  scream::Real tol = UnitTest::C::macheps
#ifdef NDEBUG
      * 10000
#endif
    ;
  std::string baseline_fn;

  // Parse options
  for (int i = 1; i < argc-1; ++i) {
    if (ekat::argv_matches(argv[i], "-g", "--generate")) generate = true;
    if (ekat::argv_matches(argv[i], "-t", "--tol")) {
      expect_another_arg(i, argc);
      ++i;
      tol = std::atof(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-b", "--baseline-file")) {
      expect_another_arg(i, argc);
      ++i;
      baseline_fn = argv[i];
    }
  }

  // Decorate baseline name with precision.
  baseline_fn += std::to_string(sizeof(scream::Real));

  std::vector<char*> args;
  for (int i=0; i<argc; ++i) {
    args.push_back(argv[i]);
  }

  scream::initialize_eamxx_session(args.size(), args.data()); {
    TestSaturation bln;
    if (generate) {
      std::cout << "Generating to " << baseline_fn << "\n";
      nerr += bln.generate_baseline(baseline_fn);
    } else {
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      nerr += bln.run_and_cmp(baseline_fn, tol);
    }
  } scream::finalize_eamxx_session();

  return nerr != 0 ? 1 : 0;
}
