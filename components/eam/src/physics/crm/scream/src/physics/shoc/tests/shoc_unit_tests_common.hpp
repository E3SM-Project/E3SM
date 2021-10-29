#ifndef SHOC_UNIT_TESTS_COMMON_HPP
#define SHOC_UNIT_TESTS_COMMON_HPP

#include "physics/shoc/shoc_functions.hpp"
#include "share/scream_types.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

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

    static constexpr Int max_pack_size = 16;
    static constexpr Int num_test_itrs = max_pack_size / Spack::n;

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
  };

};


} // namespace unit_test
} // namespace shoc
} // namespace scream

#endif
