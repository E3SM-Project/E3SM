#ifndef SHOC_UNIT_TESTS_COMMON_HPP
#define SHOC_UNIT_TESTS_COMMON_HPP

#include "shoc_functions.hpp"
#include "share/eamxx_types.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "ekat/util/ekat_file_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

/*
 * Unit test infrastructure for shoc unit tests.
 *
 * shoc entities can friend scream::shoc::unit_test::UnitWrap to give unit tests
 * access to private members.
 *
 * All unit test impls should be within an inner struct of UnitWrap::UnitTest for
 * easy access to useful types.
 */

struct UnitWrap {

  enum BASELINE_ACTION {
    NONE,
    COMPARE,
    GENERATE
  };

  template <typename D=DefaultDevice>
  struct UnitTest : public KokkosTypes<D> {

    using Device      = D;
    using MemberType  = typename KokkosTypes<Device>::MemberType;
    using TeamPolicy  = typename KokkosTypes<Device>::TeamPolicy;
    using RangePolicy = typename KokkosTypes<Device>::RangePolicy;
    using ExeSpace    = typename KokkosTypes<Device>::ExeSpace;

    template <typename S>
    using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
    template <typename S>
    using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
    template <typename S>
    using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

    template <typename S>
    using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

    using Functions          = scream::shoc::Functions<Real, Device>;
    using Scalar             = typename Functions::Scalar;
    using Spack              = typename Functions::Spack;
    using Pack               = typename Functions::Pack;
    using IntSmallPack       = typename Functions::IntSmallPack;
    using Smask              = typename Functions::Smask;
    using C                  = typename Functions::C;

    struct Base {
      std::string     m_baseline_path;
      std::string     m_test_name;
      BASELINE_ACTION m_baseline_action;
      ekat::FILEPtr   m_fid;

      Base() :
        m_baseline_path(""),
        m_test_name(Catch::getResultCapture().getCurrentTestName()),
        m_baseline_action(NONE),
        m_fid()
      {
        //Functions::shoc_init(); // just in case there is ever global shoc data
        auto& ts = ekat::TestSession::get();
        if (ts.flags["c"]) {
          m_baseline_action = COMPARE;
        }
        else if (ts.flags["g"]) {
          m_baseline_action = GENERATE;
        }
        else if (ts.flags["n"]) {
          m_baseline_action = NONE;
        }
        m_baseline_path = ts.params["b"];


        EKAT_REQUIRE_MSG( !(m_baseline_action != NONE && m_baseline_path == ""),
                          "SHOC unit test flags problem: baseline actions were requested but no baseline path was provided");

        std::string baseline_name = m_baseline_path + "/" + m_test_name;
        if (m_baseline_action == COMPARE) {
          m_fid = ekat::FILEPtr(fopen(baseline_name.c_str(), "r"));
        }
        else if (m_baseline_action == GENERATE) {
          m_fid = ekat::FILEPtr(fopen(baseline_name.c_str(), "w"));
        }
      }

      ~Base()
      {
      }

      std::mt19937_64 get_engine()
      {
        if (m_baseline_action != COMPARE) {
          // We can use any seed
          int seed;
          auto engine = setup_random_test(nullptr, &seed);
          if (m_baseline_action == GENERATE) {
            // Write the seed
            ekat::write(&seed, 1, m_fid);
          }
          return engine;
        }
        else {
          // Read the seed
          int seed;
          ekat::read(&seed, 1, m_fid);
          return setup_random_test(seed);
        }
      }
    };

    // Put struct decls here
    struct TestCalcShocVertflux;
    struct TestShocDiagObklen;
    struct TestImpCompTmpi;
    struct TestImpDpInverse;
    struct TestImpSfcFluxes;
    struct TestImpSfcStress;
    struct TestImpTkeSfcStress;
    struct TestShocUpdateDse;
    struct TestShocEnergyFixer;
    struct TestShocEnergyInt;
    struct TestShocTotEnergyFixer;
    struct TestShocEnergyDseFixer;
    struct TestShocEnergyThreshFixer;
    struct TestShocEddyDiff;
    struct TestShocGrid;
    struct TestShocCheckTke;
    struct TestShocTke;
    struct TestShocAdvSgsTke;
    struct TestShocIntColStab;
    struct TestShocIsotropicTs;
    struct TestShocShearProd;
    struct TestShocVarorCovar;
    struct TestShocLength;
    struct TestCompBruntShocLength;
    struct TestCheckShocLength;
    struct TestClipThirdMoms;
    struct TestAAdiagThirdMoms;
    struct TestW3diagThirdMoms;
    struct TestFtermInputThirdMoms;
    struct TestFtermdiagThirdMoms;
    struct TestOmegadiagThirdMoms;
    struct TestXYdiagThirdMoms;
    struct TestCompShocConvTime;
    struct TestCompShocConvVel;
    struct TestLInfShocLength;
    struct TestCompShocMixLength;
    struct TestSecondMomSrf;
    struct TestShocCompDiagThird;
    struct TestShocDiagThird;
    struct TestShocLinearInt;
    struct TestShocPdfTildetoReal;
    struct TestShocVVParameters;
    struct TestShocThlParameters;
    struct TestShocQwParameters;
    struct TestShocInPlumeCorr;
    struct TestShocAssumedPdf;
    struct TestShocPdfComputeTemp;
    struct TestShocPdfComputeQs;
    struct TestShocPdfComputeS;
    struct TestShocPdfComputeSgsLiq;
    struct TestShocPdfCompCldVar;
    struct TestShocPdfCompLiqFlux;
    struct TestShocPdfCompBuoyFlux;
    struct TestSecondMomUbycond;
    struct TestPblintdInitPot;
    struct TestDiagSecondMomentsLbycond;
    struct TestDiagSecondMoments;
    struct TestDiagSecondShocMoments;
    struct TestPblintdCldCheck;
    struct TestComputeShocVapor;
    struct TestUpdatePrognosticsImplicit;
    struct TestShocMain;
    struct TestPblintdHeight;
    struct TestVdShocDecompandSolve;
    struct TestPblintdSurfTemp;
    struct TestPblintdCheckPblh;
    struct TestPblintd;
    struct TestComputeShocTemp;
  };

};


} // namespace unit_test
} // namespace shoc
} // namespace scream

#endif
