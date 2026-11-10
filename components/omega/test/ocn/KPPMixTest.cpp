//===-- Test driver for OMEGA KPP mixing -------------------------*- C++
//-*-===//
///
/// \file
/// \brief Unit tests for KPP profiles, utilities, and mixing coefficients
///
//===----------------------------------------------------------------------===//

#include "KPPMix.h"
#include "Config.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "FillValues.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "KPPConstants.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include "mpi.h"

#include <memory>
#include <string>

using namespace OMEGA;

namespace {

#ifdef SINGLE_PRECISION
constexpr Real RTol    = 2.0e-5_Real;
constexpr Real ATol    = 2.0e-6_Real;
constexpr Real BLDRTol = 2.0e-4_Real;
#else
constexpr Real RTol    = 2.0e-10_Real;
constexpr Real ATol    = 2.0e-12_Real;
constexpr Real BLDRTol = RTol;
#endif

constexpr Real LayerThickness = 10.0_Real;
constexpr Real TestOSBLDepth  = 40.0_Real;
constexpr I4 TestOSBLIndex    = 3;

std::unique_ptr<Clock> TestClock;
Clock *TestClockPtr = nullptr;

void initKPPMixTest(const std::string &TestGroup) {
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();
   initLogging(DefEnv);
   LOG_INFO("------ KPP Mixing Unit Tests ------");

   Config("Omega");
   Config::readAll("omega.yml");
   // These groups inject a rejected MatchTechnique; KPPMix::init must abort,
   // so the ctest entries for them are registered as expected failures.
   if (TestGroup == "config-gradient" || TestGroup == "config-unsupported" ||
       TestGroup == "config-parabolic") {
      Config VertMixConfig("VertMix");
      Config KPPConfig("KPP");
      Error Err;
      Err += Config::getOmegaConfig()->get(VertMixConfig);
      Err += VertMixConfig.get(KPPConfig);
      CHECK_ERROR_ABORT(Err, "KPPMixTest: unable to access KPP configuration");
      const std::string MatchTechnique =
          TestGroup == "config-gradient"    ? "MatchGradient"
          : TestGroup == "config-parabolic" ? "ParabolicNonLocal"
                                            : "NotAKPPMode";
      KPPConfig.set("MatchTechnique", MatchTechnique);
   }
   IO::init(DefComm);
   Decomp::init();
   Halo::init();

   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   TestClock    = std::make_unique<Clock>(StartTime, TimeStep);
   TestClockPtr = TestClock.get();
   Field::init(TestClockPtr);
   IOStream::init(TestClockPtr);
   HorzMesh::init(TestClockPtr);
   VertCoord::init(false);
   KPPMix::init();
}

void finalizeKPPMixTest() {
   KPPMix::destroyInstance();
   IOStream::finalize();
   VertCoord::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   TestClockPtr = nullptr;
   TestClock.reset();
   MachEnv::removeAll();
}

void checkResult(const char *TestName, int NumErrors) {
   if (NumErrors != 0) {
      ABORT_ERROR("KPPMixTest: {} FAIL with {} errors", TestName, NumErrors);
   }
   LOG_INFO("KPPMixTest: {} PASS", TestName);
}

void testStabilityFunctions() {
   constexpr int NTests = 9;
   int NumErrors        = 0;

   parallelReduce(
       "KPPMixTest-StabilityFunctions", {NTests},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          Real Zeta = 0.0_Real;
          switch (ITest) {
          case 0:
             Zeta = 2.0_Real;
             break;
          case 1:
             Zeta = 0.1_Real;
             break;
          case 2:
             Zeta = 0.0_Real;
             break;
          case 3:
             Zeta = -0.1_Real;
             break;
          case 4:
             Zeta = KPP::ZetaM;
             break;
          case 5:
             Zeta = KPP::ZetaM - 1.0e-4_Real;
             break;
          case 6:
             Zeta = KPP::ZetaS;
             break;
          case 7:
             Zeta = KPP::ZetaS - 1.0e-4_Real;
             break;
          default:
             Zeta = -10.0_Real;
             break;
          }

          Real ExpectedM;
          if (Zeta >= 0.0_Real) {
             ExpectedM = 1.0_Real / (1.0_Real + 5.0_Real * Zeta);
          } else if (Zeta >= KPP::ZetaM) {
             ExpectedM = Kokkos::pow(1.0_Real - 16.0_Real * Zeta, 0.25_Real);
          } else {
             ExpectedM =
                 Kokkos::pow(KPP::AMoM - KPP::CMoM * Zeta, 1.0_Real / 3.0_Real);
          }

          Real ExpectedS;
          if (Zeta >= 0.0_Real) {
             ExpectedS = 1.0_Real / (1.0_Real + 5.0_Real * Zeta);
          } else if (Zeta >= KPP::ZetaS) {
             ExpectedS = Kokkos::sqrt(1.0_Real - 16.0_Real * Zeta);
          } else {
             ExpectedS =
                 Kokkos::pow(KPP::AMoS - KPP::CMoS * Zeta, 1.0_Real / 3.0_Real);
          }

          const Real ActualM = KPP::kppPhiInvMomentum(Zeta);
          const Real ActualS = KPP::kppPhiInvScalar(Zeta);
          if (!isApprox(ActualM, ExpectedM, RTol, ATol) || ActualM <= 0.0_Real)
             ++ErrorCount;
          if (!isApprox(ActualS, ExpectedS, RTol, ATol) || ActualS <= 0.0_Real)
             ++ErrorCount;
       },
       NumErrors);

   checkResult("stability functions", NumErrors);

   NumErrors = 0;
   parallelReduce(
       "KPPMixTest-StabilityContinuity", {2},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          const Real Transition = ITest == 0 ? KPP::ZetaM : KPP::ZetaS;
          const Real Epsilon    = 1.0e-6_Real;
          const Real Above      = ITest == 0
                                      ? KPP::kppPhiInvMomentum(Transition + Epsilon)
                                      : KPP::kppPhiInvScalar(Transition + Epsilon);
          const Real Below      = ITest == 0
                                      ? KPP::kppPhiInvMomentum(Transition - Epsilon)
                                      : KPP::kppPhiInvScalar(Transition - Epsilon);
          if (!isApprox(Above, Below, 2.0e-5_Real, 2.0e-5_Real))
             ++ErrorCount;
       },
       NumErrors);
   checkResult("stability transition continuity", NumErrors);
}

void testShapeFunctions() {
   constexpr int NTests = 7;
   int NumErrors        = 0;

   parallelReduce(
       "KPPMixTest-ShapeFunctions", {NTests},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          Real Sigma;
          switch (ITest) {
          case 0:
             Sigma = 0.25_Real;
             break;
          case 1:
             Sigma = 0.0_Real;
             break;
          case 2:
             Sigma = -0.25_Real;
             break;
          case 3:
             Sigma = -0.5_Real;
             break;
          case 4:
             Sigma = -0.75_Real;
             break;
          case 5:
             Sigma = -1.0_Real;
             break;
          default:
             Sigma = -1.25_Real;
             break;
          }

          const Real SigmaClamped =
              Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));
          const Real SigmaMu         = -SigmaClamped;
          const Real OneMinus        = 1.0_Real - SigmaMu;
          const Real ExpectedSimple  = SigmaMu * OneMinus * OneMinus;
          constexpr Real ShapeAtBase = 0.125_Real;
          const Real Smooth =
              SigmaMu * SigmaMu * (3.0_Real - 2.0_Real * SigmaMu);

          if (!isApprox(KPP::kppShapeMomentum(Sigma), ExpectedSimple, RTol,
                        ATol))
             ++ErrorCount;
          if (!isApprox(KPP::kppShapeScalar(Sigma), ExpectedSimple, RTol, ATol))
             ++ErrorCount;
          if (!isApprox(KPP::kppShapeMatched(Sigma, ShapeAtBase),
                        ExpectedSimple + ShapeAtBase * Smooth, RTol, ATol))
             ++ErrorCount;
          if (!isApprox(KPP::kppShapeMatched(Sigma, 0.0_Real), ExpectedSimple,
                        RTol, ATol))
             ++ErrorCount;
          if (!isApprox(KPP::kppSurfaceMomentumScale(Sigma),
                        KPP::HuOn * (1.0_Real + Sigma), RTol, ATol))
             ++ErrorCount;
       },
       NumErrors);

   checkResult("shape functions", NumErrors);
}

void testLangmuirFunctions() {
   int NumErrors = 0;

   parallelReduce(
       "KPPMixTest-LangmuirFunctions", {4},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          const Real Wind           = ITest == 0   ? -5.0_Real
                                      : ITest == 1 ? 0.0_Real
                                      : ITest == 2 ? 10.0_Real
                                                   : 100.0_Real;
          const Real UStar          = ITest < 2 ? 0.0_Real : 0.01_Real;
          const Real WindClamped    = Kokkos::fmax(0.0_Real, Wind);
          const Real ExpectedStokes = 0.016_Real * WindClamped;
          const Real UStarClamped   = Kokkos::fmax(KPP::MinUStar, UStar);
          const Real StokesClamped  = Kokkos::fmax(1.0e-8_Real, ExpectedStokes);
          const Real ExpectedLa = Kokkos::sqrt(UStarClamped / StokesClamped);
          const Real LaInv      = 1.0_Real / Kokkos::fmax(0.5_Real, ExpectedLa);
          const Real ExpectedEnhancement = Kokkos::fmin(
              2.0_Real,
              Kokkos::fmax(1.0_Real,
                           Kokkos::sqrt(1.0_Real + 0.5_Real * LaInv * LaInv)));

          const Real Stokes = KPP::estimateStokesDriftSL(Wind, 50.0_Real);
          const Real La     = KPP::computeLangmuirNumber(UStar, Stokes);
          const Real Enhancement =
              KPP::computeLangmuirEnhancement(Wind, UStar, 50.0_Real);
          if (!isApprox(Stokes, ExpectedStokes, RTol, ATol))
             ++ErrorCount;
          if (!isApprox(La, ExpectedLa, RTol, ATol))
             ++ErrorCount;
          if (!isApprox(Enhancement, ExpectedEnhancement, RTol, ATol) ||
              Enhancement < 1.0_Real || Enhancement > 2.0_Real)
             ++ErrorCount;
       },
       NumErrors);

   checkResult("Langmuir functions", NumErrors);
}

void testOSBLUtilities() {
   int NumErrors = 0;

   parallelReduce(
       "KPPMixTest-OSBLUtilities", {4},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          const Real IceFraction = ITest == 0   ? 0.0_Real
                                   : ITest == 1 ? KPP::IceSuppressThresh
                                   : ITest == 2
                                       ? KPP::IceSuppressThresh + 0.01_Real
                                       : 0.0_Real;
          const I4 LandIceMask   = ITest == 3 ? 1 : 0;
          const bool ExpectedSuppression =
              LandIceMask != 0 || IceFraction > KPP::IceSuppressThresh;
          if (KPP::shouldSuppressOSBL(IceFraction, LandIceMask) !=
              ExpectedSuppression)
             ++ErrorCount;

          const Real InputDepth = ITest == 0   ? 1.0_Real
                                  : ITest == 1 ? 20.0_Real
                                  : ITest == 2 ? 1.0_Real
                                               : 200.0_Real;
          Real ExpectedDepth    = Kokkos::fmax(InputDepth, 4.0_Real);
          if (IceFraction > KPP::IceSuppressThresh)
             ExpectedDepth = Kokkos::fmax(ExpectedDepth, KPP::MinOSBLUnderIce);
          ExpectedDepth = Kokkos::fmin(ExpectedDepth, 95.0_Real);
          if (!isApprox(KPP::constrainOSBLDepth(InputDepth, 4.0_Real,
                                                100.0_Real, IceFraction),
                        ExpectedDepth, RTol, ATol))
             ++ErrorCount;
       },
       NumErrors);

   checkResult("OSBL utilities", NumErrors);
}

void testTurbulentVelocityScale() {
   int NumErrors = 0;

   parallelReduce(
       "KPPMixTest-TurbulentVelocityScale", {6},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          const Real UStar        = ITest == 0   ? 0.02_Real
                                    : ITest == 1 ? 0.0_Real
                                    : ITest == 2 ? 0.02_Real
                                    : ITest == 5 ? -0.02_Real
                                                 : 0.0_Real;
          const Real B0           = ITest == 0   ? 0.0_Real
                                    : ITest == 1 ? -1.0e-7_Real
                                    : ITest == 2 ? -1.0e-7_Real
                                    : ITest == 3 ? 1.0e-7_Real
                                                 : 0.0_Real;
          const Real H            = ITest == 5 ? -50.0_Real : 50.0_Real;
          const Real UStarClamped = Kokkos::fmax(0.0_Real, UStar);
          const Real HClamped     = Kokkos::fmax(0.0_Real, H);
          const Real Momentum     = UStarClamped * UStarClamped * UStarClamped;
          const Real Buoyancy =
              KPP::ConvectiveVelFac * Kokkos::fmax(0.0_Real, -B0) * HClamped;
          const Real Expected =
              Kokkos::pow(Momentum + Buoyancy, 1.0_Real / 3.0_Real);
          const Real Actual = KPP::computeTurbVelocityScale(UStar, B0, H);
          if (!isApprox(Actual, Expected, RTol, ATol) || Actual < 0.0_Real)
             ++ErrorCount;
       },
       NumErrors);

   checkResult("turbulent velocity scale", NumErrors);
}

void testTurbScales() {
   int NumErrors = 0;

   parallelReduce(
       "KPPMixTest-TurbScales", {5},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          // 0: wind only, 1: wind + unstable, 2: wind + stable,
          // 3: free convection (u*=0, B0<0), 4: calm and stable (both zero)
          const Real UStar = (ITest == 3 || ITest == 4) ? 0.0_Real : 0.02_Real;
          const Real B0    = ITest == 0   ? 0.0_Real
                             : ITest == 1 ? -1.0e-7_Real
                             : ITest == 2 ? 1.0e-7_Real
                             : ITest == 3 ? -1.0e-7_Real
                                          : 1.0e-7_Real;
          const Real HOSBL = 50.0_Real;
          const Real SigmaLoc = KPP::SurfaceLayerExtent;

          Real WM = -1.0_Real;
          Real WS = -1.0_Real;
          KPP::kppTurbScales(UStar, B0, HOSBL, SigmaLoc, VonKar, WM, WS);

          Real ExpectedWM = 0.0_Real;
          Real ExpectedWS = 0.0_Real;
          if (UStar > 0.0_Real) {
             const Real U3 = UStar * UStar * UStar;
             const Real Zeta =
                 SigmaLoc * HOSBL * B0 * VonKar / Kokkos::max(U3, 1.0e-20_Real);
             ExpectedWM = VonKar * UStar *
                          Kokkos::max(KPP::kppPhiInvMomentum(Zeta), 0.0_Real);
             ExpectedWS = VonKar * UStar *
                          Kokkos::max(KPP::kppPhiInvScalar(Zeta), 0.0_Real);
          } else if (B0 < 0.0_Real) {
             const Real WM3 = -KPP::CMoM * SigmaLoc * HOSBL * VonKar * B0;
             const Real WS3 = -KPP::CMoS * SigmaLoc * HOSBL * VonKar * B0;
             ExpectedWM     = VonKar * Kokkos::pow(WM3, 1.0_Real / 3.0_Real);
             ExpectedWS     = VonKar * Kokkos::pow(WS3, 1.0_Real / 3.0_Real);
          }

          if (!isApprox(WM, ExpectedWM, RTol, ATol) ||
              !isApprox(WS, ExpectedWS, RTol, ATol))
             ++ErrorCount;

          // Scales must never go negative, and scalars mix at least as
          // efficiently as momentum under unstable forcing.
          if (WM < 0.0_Real || WS < 0.0_Real)
             ++ErrorCount;
          if (B0 < 0.0_Real && WS < WM)
             ++ErrorCount;

          // Calm and stable is fully quiescent.
          if (ITest == 4 && (WM != 0.0_Real || WS != 0.0_Real))
             ++ErrorCount;
       },
       NumErrors);

   checkResult("turbulent scales", NumErrors);
}

void testMatchShape() {
   int NumErrors = 0;

   parallelReduce(
       "KPPMixTest-MatchShape", {4},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          const Real Interior = 1.0e-4_Real;
          const Real HOSBL    = ITest == 1 ? 0.0_Real : 50.0_Real;
          const Real W        = ITest == 2 ? 0.0_Real : 0.01_Real;

          const Real Shape = KPP::kppMatchShape(Interior, HOSBL, W);

          // Degenerate HOSBL or velocity scale switches matching off.
          if (ITest == 1 || ITest == 2) {
             if (Shape != 0.0_Real)
                ++ErrorCount;
             return;
          }

          if (!isApprox(Shape, Interior / (HOSBL * W), RTol, ATol))
             ++ErrorCount;

          // The matched shape must reproduce the interior coefficient when
          // multiplied back by h*w, which is the whole point of matching.
          if (!isApprox(HOSBL * W * Shape, Interior, RTol, ATol))
             ++ErrorCount;
       },
       NumErrors);

   checkResult("match shape", NumErrors);
}

void testNonLocalCs() {
   int NumErrors = 0;

   parallelReduce(
       "KPPMixTest-NonLocalCs", {1},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          const Real Cs = KPP::kppNonLocalCs(VonKar, KPP::SurfaceLayerExtent);
          const Real Expected =
              10.0_Real * VonKar *
              Kokkos::pow(KPP::CMoS * VonKar * KPP::SurfaceLayerExtent,
                          1.0_Real / 3.0_Real);
          if (!isApprox(Cs, Expected, RTol, ATol))
             ++ErrorCount;

          // Large et al. (1994) quote C_s ~ 6.33 for the default constants.
          if (Kokkos::abs(Cs - 6.33_Real) > 0.05_Real)
             ++ErrorCount;
       },
       NumErrors);

   checkResult("non-local flux constant", NumErrors);
}

void testClampOSBLDepth() {
   int NumErrors = 0;

   parallelReduce(
       "KPPMixTest-ClampOSBLDepth", {5},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          const Real MinDepth = 2.0_Real;
          const Real MaxDepth = 95.0_Real;
          const Real Input    = ITest == 0   ? 40.0_Real
                                : ITest == 1 ? 1.0_Real
                                : ITest == 2 ? 200.0_Real
                                : ITest == 3 ? 1.0_Real
                                             : 1.0_Real;
          const bool ApplyIce = (ITest == 3);

          const Real Actual = KPP::kppClampOSBLDepth(
              Input, MinDepth, MaxDepth, ApplyIce, KPP::MinOSBLUnderIce);

          Real Expected = Kokkos::fmax(Input, MinDepth);
          if (ApplyIce)
             Expected = Kokkos::fmax(Expected, KPP::MinOSBLUnderIce);
          Expected = Kokkos::fmin(Expected, MaxDepth);

          if (!isApprox(Actual, Expected, RTol, ATol))
             ++ErrorCount;

          // Result must always land inside the supported range.
          if (Actual < MinDepth || Actual > MaxDepth)
             ++ErrorCount;
       },
       NumErrors);

   checkResult("OSBL depth clamping", NumErrors);
}

void testOSBLIndex() {
   int NumErrors = 0;

   // Single column of five 10 m layers, interfaces at 0,-10,...,-50 m.
   constexpr I4 NLayers = 5;
   Array2DReal ZInterface("KPPMixTest-OSBLIndexZ", 1, NLayers + 1);
   parallelFor(
       "KPPMixTest-OSBLIndexInit", {NLayers + 1},
       KOKKOS_LAMBDA(I4 K) { ZInterface(0, K) = -10.0_Real * K; });

   parallelReduce(
       "KPPMixTest-OSBLIndex", {5},
       KOKKOS_LAMBDA(int ITest, int &ErrorCount) {
          // Cases 0-3 bracket a layer; case 4 is deeper than the column.
          const Real Depth  = ITest == 0   ? 5.0_Real
                              : ITest == 1 ? 25.0_Real
                              : ITest == 2 ? 20.0_Real
                              : ITest == 3 ? 0.0_Real
                                           : 100.0_Real;
          const I4 Expected = ITest == 0   ? 0
                              : ITest == 1 ? 2
                              : ITest == 2 ? 1
                              : ITest == 3 ? 0
                                           : NLayers - 1;

          const I4 Actual =
              KPP::kppOSBLIndex(ZInterface, 0, 0, NLayers - 1, 0.0_Real, Depth);
          if (Actual != Expected)
             ++ErrorCount;
       },
       NumErrors);

   checkResult("OSBL index lookup", NumErrors);
}

// Builds a uniform column whose free surface sits at Ssh. All KPP results must
// be invariant to Ssh since depths are measured below the free surface.
void setCoefficientTestGeometry(Real Ssh = 0.0_Real) {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;

   OMEGA_SCOPE(GeomZInterface, VCoord->GeomZInterface);
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(SshCell, VCoord->SshCell);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(OSBLDepth, KPPInstance->OSBLDepth);
   OMEGA_SCOPE(OSBLDepthIndex, KPPInstance->OSBLDepthIndex);

   parallelFor(
       "KPPMixTest-SetGeometry", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell)   = 0;
          MaxLayerCell(ICell)   = NVertLayers - 1;
          SshCell(ICell)        = Ssh;
          OSBLDepth(ICell)      = TestOSBLDepth;
          OSBLDepthIndex(ICell) = TestOSBLIndex;
          for (I4 K = 0; K <= NVertLayers; ++K) {
             GeomZInterface(ICell, K) = Ssh - LayerThickness * K;
             if (K < NVertLayers) {
                GeomZMid(ICell, K) = Ssh - LayerThickness * (K + 0.5_Real);
             }
          }
       });
}

Real nonLocalNormalization() {
   return 10.0_Real * VonKar *
          Kokkos::pow(KPP::CMoS * VonKar * KPP::SurfaceLayerExtent,
                      1.0_Real / 3.0_Real);
}

void testWindOnlyCoefficients() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-Density", Mesh->NCellsSize,
                       VCoord->NVertLayers);
   Array1DReal UStar("KPPMixTest-UStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-B0", Mesh->NCellsSize);
   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);

   KPPInstance->UseEnhancedDiffusion = false;
   KPPInstance->MatchTechnique       = KPPMatchType::SimpleShapes;
   KPPInstance->computeMixingCoefficients(UStar, B0);

   const auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   const auto NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   const auto TurbVelH =
       createHostMirrorCopy(KPPInstance->TurbulentVelocityScale);

   const Real Sigma            = -0.5_Real;
   const Real Shape            = 0.125_Real;
   const Real TurbVel          = VonKar * 0.02_Real;
   const Real ExpectedMix      = TestOSBLDepth * TurbVel * Shape;
   const Real ExpectedNonLocal = nonLocalNormalization() * Shape;
   int NumErrors               = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(VertDiffH(ICell, 2), ExpectedMix, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 2), ExpectedMix, RTol, ATol) ||
          !isApprox(TurbVelH(ICell, 2), TurbVel, RTol, ATol) ||
          !isApprox(NonLocalH(ICell, 2), ExpectedNonLocal, RTol, ATol)) {
         ++NumErrors;
      }
      if (!isApprox(KPP::kppShapeMomentum(Sigma), Shape, RTol, ATol) ||
          VertDiffH(ICell, 0) != 0.0_Real || VertViscH(ICell, 0) != 0.0_Real ||
          VertDiffH(ICell, 4) != 0.0_Real || VertViscH(ICell, 4) != 0.0_Real ||
          VertDiffH(ICell, 5) != 0.0_Real || NonLocalH(ICell, 5) != 0.0_Real) {
         ++NumErrors;
      }
   }
   checkResult("wind-only coefficients", NumErrors);
}

void testConvectionOnlyCoefficients() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-ConvDensity", Mesh->NCellsSize,
                       VCoord->NVertLayers);
   Array1DReal UStar("KPPMixTest-ConvUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-ConvB0", Mesh->NCellsSize);
   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.0_Real);
   deepCopy(B0, -1.0e-7_Real);

   KPPInstance->UseEnhancedDiffusion = false;
   KPPInstance->MatchTechnique       = KPPMatchType::SimpleShapes;
   KPPInstance->computeMixingCoefficients(UStar, B0);

   const auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   const auto TurbVelH =
       createHostMirrorCopy(KPPInstance->TurbulentVelocityScale);
   const Real SigmaLoc = KPP::SurfaceLayerExtent;
   const Real WM = VonKar * Kokkos::pow(KPP::CMoM * SigmaLoc * TestOSBLDepth *
                                            VonKar * 1.0e-7_Real,
                                        1.0_Real / 3.0_Real);
   const Real WS = VonKar * Kokkos::pow(KPP::CMoS * SigmaLoc * TestOSBLDepth *
                                            VonKar * 1.0e-7_Real,
                                        1.0_Real / 3.0_Real);
   constexpr Real Shape    = 0.125_Real;
   const Real ExpectedVisc = TestOSBLDepth * WM * Shape;
   const Real ExpectedDiff = TestOSBLDepth * WS * Shape;

   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(VertViscH(ICell, 2), ExpectedVisc, RTol, ATol) ||
          !isApprox(VertDiffH(ICell, 2), ExpectedDiff, RTol, ATol) ||
          !isApprox(TurbVelH(ICell, 2), WS, RTol, ATol) ||
          !(VertDiffH(ICell, 2) > VertViscH(ICell, 2))) {
         ++NumErrors;
      }
   }
   checkResult("convection-only coefficients", NumErrors);
}

void testNonLocalProfileModes() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-ProfileDensity", Mesh->NCellsSize,
                       VCoord->NVertLayers);
   Array1DReal UStar("KPPMixTest-ProfileUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-ProfileB0", Mesh->NCellsSize);
   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);

   KPPInstance->UseEnhancedDiffusion = false;
   const Real Normalization          = nonLocalNormalization();
   int NumErrors                     = 0;

   // The non-local flux follows the scalar diffusivity shape, so at sigma=-0.5
   // it is 0.125 * C_s and it vanishes at the surface and at the OSBL base.
   KPPInstance->MatchTechnique = KPPMatchType::SimpleShapes;
   KPPInstance->computeMixingCoefficients(UStar, B0);
   auto NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (NonLocalH(ICell, 0) != 0.0_Real ||
          !isApprox(NonLocalH(ICell, 2), 0.125_Real * Normalization, RTol,
                    ATol) ||
          NonLocalH(ICell, 4) != 0.0_Real) {
         ++NumErrors;
      }
   }

   // Without interior coefficients there is nothing to match, so MatchBoth
   // must reduce exactly to SimpleShapes.
   const auto SimpleShapesH    = NonLocalH;
   KPPInstance->MatchTechnique = KPPMatchType::MatchBoth;
   KPPInstance->computeMixingCoefficients(UStar, B0);
   NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      for (I4 K = 0; K <= VCoord->NVertLayers; ++K) {
         if (NonLocalH(ICell, K) != SimpleShapesH(ICell, K)) {
            ++NumErrors;
         }
      }
   }

   checkResult("non-local profile modes", NumErrors);
}

void testMatchBothInteriorCoefficients() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-MatchDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array1DReal UStar("KPPMixTest-MatchUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-MatchB0", Mesh->NCellsSize);
   Array2DReal InteriorDiff("KPPMixTest-InteriorDiff", Mesh->NCellsSize,
                            NVertLayers + 1);
   Array2DReal InteriorVisc("KPPMixTest-InteriorVisc", Mesh->NCellsSize,
                            NVertLayers + 1);
   constexpr Real ExpectedInteriorDiff = 2.0e-3_Real;
   constexpr Real ExpectedInteriorVisc = 4.0e-3_Real;
   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);
   deepCopy(InteriorDiff, ExpectedInteriorDiff);
   deepCopy(InteriorVisc, ExpectedInteriorVisc);

   KPPInstance->UseEnhancedDiffusion = false;
   KPPInstance->MatchTechnique       = KPPMatchType::MatchBoth;
   KPPInstance->computeMixingCoefficients(UStar, B0, InteriorDiff,
                                          InteriorVisc);

   const auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   const auto NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   constexpr Real Sigma = -0.5_Real;
   constexpr Real SimpleShape   = 0.125_Real;
   constexpr Real SmoothAtSigma = 0.5_Real;
   const Real TurbVel           = VonKar * 0.02_Real;
   const Real ExpectedDiffMid   = TestOSBLDepth * TurbVel * SimpleShape +
                                SmoothAtSigma * ExpectedInteriorDiff;
   const Real ExpectedViscMid = TestOSBLDepth * TurbVel * SimpleShape +
                                SmoothAtSigma * ExpectedInteriorVisc;
   const Real MatchDiffShape = ExpectedInteriorDiff / (TestOSBLDepth * TurbVel);
   // The non-local shape is independent of MatchTechnique, so gamma still
   // follows the unmatched scalar shape and vanishes at the OSBL base.
   const Real ExpectedNonLocal =
       nonLocalNormalization() * KPP::kppShapeScalar(Sigma);

   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(VertDiffH(ICell, 2), ExpectedDiffMid, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 2), ExpectedViscMid, RTol, ATol) ||
          !isApprox(NonLocalH(ICell, 2), ExpectedNonLocal, RTol, ATol) ||
          !isApprox(VertDiffH(ICell, 4), ExpectedInteriorDiff, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 4), ExpectedInteriorVisc, RTol, ATol) ||
          !isApprox(VertDiffH(ICell, 5), ExpectedInteriorDiff, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 5), ExpectedInteriorVisc, RTol, ATol) ||
          NonLocalH(ICell, 4) != 0.0_Real || NonLocalH(ICell, 5) != 0.0_Real) {
         ++NumErrors;
      }
   }
   checkResult("MatchBoth interior coefficients", NumErrors);
}

void testEnhancedDiffusion() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-EnhancedDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array1DReal UStar("KPPMixTest-EnhancedUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-EnhancedB0", Mesh->NCellsSize);
   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);

   KPPInstance->MatchTechnique       = KPPMatchType::SimpleShapes;
   KPPInstance->UseEnhancedDiffusion = false;
   KPPInstance->computeMixingCoefficients(UStar, B0);
   auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   int NumErrors  = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (VertDiffH(ICell, 4) != 0.0_Real) {
         ++NumErrors;
      }
   }

   KPPInstance->UseEnhancedDiffusion = true;
   KPPInstance->computeMixingCoefficients(UStar, B0);
   VertDiffH      = createHostMirrorCopy(KPPInstance->VertDiff);
   auto VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   auto NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   constexpr Real OutsideDelta = 0.5_Real;
   const Real OutsideSigma     = -35.0_Real / TestOSBLDepth;
   const Real OutsideProfile =
       TestOSBLDepth * VonKar * 0.02_Real * KPP::kppShapeScalar(OutsideSigma);
   const Real ExpectedOutside = OutsideDelta * (1.0_Real - OutsideDelta) *
                                (1.0_Real - OutsideDelta) * OutsideProfile;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(VertDiffH(ICell, 4), ExpectedOutside, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 4), ExpectedOutside, RTol, ATol) ||
          NonLocalH(ICell, 4) != 0.0_Real) {
         ++NumErrors;
      }
   }

   constexpr Real InsideOSBLDepth = 32.0_Real;
   deepCopy(KPPInstance->OSBLDepth, InsideOSBLDepth);
   deepCopy(KPPInstance->OSBLDepthIndex, TestOSBLIndex);
   KPPInstance->computeMixingCoefficients(UStar, B0);
   VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);

   constexpr Real InsideDelta         = 0.7_Real;
   constexpr Real OneMinusInsideDelta = 1.0_Real - InsideDelta;
   const Real KtupSigma               = -25.0_Real / InsideOSBLDepth;
   const Real TargetSigma             = -30.0_Real / InsideOSBLDepth;
   const Real KtupProfile =
       InsideOSBLDepth * VonKar * 0.02_Real * KPP::kppShapeScalar(KtupSigma);
   const Real TargetProfile =
       InsideOSBLDepth * VonKar * 0.02_Real * KPP::kppShapeScalar(TargetSigma);
   const Real ExpectedInside =
       InsideDelta * (OneMinusInsideDelta * OneMinusInsideDelta * KtupProfile +
                      InsideDelta * InsideDelta * TargetProfile);
   const Real ExpectedInsideNonLocal = nonLocalNormalization() *
                                       KPP::kppShapeScalar(TargetSigma) *
                                       ExpectedInside / TargetProfile;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(VertDiffH(ICell, 3), ExpectedInside, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 3), ExpectedInside, RTol, ATol) ||
          !isApprox(NonLocalH(ICell, 3), ExpectedInsideNonLocal, RTol, ATol)) {
         ++NumErrors;
      }
   }

   Array2DReal InteriorDiff("KPPMixTest-EnhancedInteriorDiff", Mesh->NCellsSize,
                            NVertLayers + 1);
   Array2DReal InteriorVisc("KPPMixTest-EnhancedInteriorVisc", Mesh->NCellsSize,
                            NVertLayers + 1);
   constexpr Real InteriorDiffValue = 2.0e-3_Real;
   constexpr Real InteriorViscValue = 4.0e-3_Real;
   deepCopy(InteriorDiff, InteriorDiffValue);
   deepCopy(InteriorVisc, InteriorViscValue);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);
   deepCopy(KPPInstance->OSBLDepth, TestOSBLDepth);
   deepCopy(KPPInstance->OSBLDepthIndex, TestOSBLIndex);
   KPPInstance->MatchTechnique = KPPMatchType::MatchBoth;
   KPPInstance->computeMixingCoefficients(UStar, B0, InteriorDiff,
                                          InteriorVisc);
   VertDiffH                    = createHostMirrorCopy(KPPInstance->VertDiff);
   VertViscH                    = createHostMirrorCopy(KPPInstance->VertVisc);
   const Real InteriorKtupSigma = -35.0_Real / TestOSBLDepth;
   const Real DiffMatchShape =
       InteriorDiffValue / (TestOSBLDepth * VonKar * 0.02_Real);
   const Real ViscMatchShape =
       InteriorViscValue / (TestOSBLDepth * VonKar * 0.02_Real);
   const Real DiffKtup =
       TestOSBLDepth * VonKar * 0.02_Real *
       KPP::kppShapeMatched(InteriorKtupSigma, DiffMatchShape);
   const Real ViscKtup =
       TestOSBLDepth * VonKar * 0.02_Real *
       KPP::kppShapeMatched(InteriorKtupSigma, ViscMatchShape);
   const Real ExpectedInteriorEnhancedDiff =
       0.625_Real * InteriorDiffValue + 0.125_Real * DiffKtup;
   const Real ExpectedInteriorEnhancedVisc =
       0.625_Real * InteriorViscValue + 0.125_Real * ViscKtup;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(VertDiffH(ICell, 4), ExpectedInteriorEnhancedDiff, RTol,
                    ATol) ||
          !isApprox(VertViscH(ICell, 4), ExpectedInteriorEnhancedVisc, RTol,
                    ATol)) {
         ++NumErrors;
      }
   }

   deepCopy(UStar, 0.0_Real);
   deepCopy(B0, 0.0_Real);
   KPPInstance->MatchTechnique = KPPMatchType::SimpleShapes;
   deepCopy(KPPInstance->OSBLDepth, InsideOSBLDepth);
   deepCopy(KPPInstance->OSBLDepthIndex, TestOSBLIndex);
   KPPInstance->computeMixingCoefficients(UStar, B0);
   VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (VertDiffH(ICell, 3) != 0.0_Real || VertViscH(ICell, 3) != 0.0_Real ||
          NonLocalH(ICell, 3) != 0.0_Real) {
         ++NumErrors;
      }
   }
   checkResult("enhanced diffusion", NumErrors);
}

void testStableAndZeroForcing() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-StableDensity", Mesh->NCellsSize,
                       VCoord->NVertLayers);
   Array1DReal UStar("KPPMixTest-StableUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-StableB0", Mesh->NCellsSize);
   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 1.0e-7_Real);

   KPPInstance->UseEnhancedDiffusion = false;
   KPPInstance->MatchTechnique       = KPPMatchType::SimpleShapes;
   KPPInstance->computeMixingCoefficients(UStar, B0);
   auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   auto NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   int NumErrors  = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!(VertDiffH(ICell, 2) > 0.0_Real) ||
          NonLocalH(ICell, 2) != 0.0_Real) {
         ++NumErrors;
      }
   }

   deepCopy(UStar, 0.0_Real);
   deepCopy(B0, 0.0_Real);
   KPPInstance->computeMixingCoefficients(UStar, B0);
   VertDiffH            = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (VertDiffH(ICell, 2) != 0.0_Real || VertViscH(ICell, 2) != 0.0_Real) {
         ++NumErrors;
      }
   }
   checkResult("stable and zero forcing", NumErrors);
}

void testCoefficientVerticalDomainEdges() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-EdgeDensity", Mesh->NCellsSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-EdgeUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-EdgeB0", Mesh->NCellsSize);
   deepCopy(Density, -99.0_Real);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(OSBLDepth, KPPInstance->OSBLDepth);
   OMEGA_SCOPE(OSBLDepthIndex, KPPInstance->OSBLDepthIndex);
   parallelFor(
       "KPPMixTest-SetPartialColumn", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell)   = 2;
          MaxLayerCell(ICell)   = 4;
          OSBLDepth(ICell)      = 40.0_Real;
          OSBLDepthIndex(ICell) = 3;
       });

   KPPInstance->UseEnhancedDiffusion = false;
   KPPInstance->MatchTechnique       = KPPMatchType::SimpleShapes;
   KPPInstance->computeMixingCoefficients(UStar, B0);
   auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   auto VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   auto NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   const Real ExpectedPartial =
       40.0_Real * VonKar * 0.02_Real * KPP::kppShapeScalar(-0.5_Real);
   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (VertDiffH(ICell, 0) != 0.0_Real || VertDiffH(ICell, 1) != 0.0_Real ||
          !isApprox(VertDiffH(ICell, 2), ExpectedPartial, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 2), ExpectedPartial, RTol, ATol) ||
          !(NonLocalH(ICell, 2) > 0.0_Real) ||
          VertDiffH(ICell, 4) != 0.0_Real || VertDiffH(ICell, 5) != 0.0_Real) {
         ++NumErrors;
      }
   }

   parallelFor(
       "KPPMixTest-SetOneLayerColumn", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell)   = 2;
          MaxLayerCell(ICell)   = 2;
          OSBLDepth(ICell)      = 25.0_Real;
          OSBLDepthIndex(ICell) = 2;
       });
   KPPInstance->computeMixingCoefficients(UStar, B0);
   VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   const Real ExpectedOneLayer =
       25.0_Real * VonKar * 0.02_Real * KPP::kppShapeScalar(-0.8_Real);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(VertDiffH(ICell, 2), ExpectedOneLayer, RTol, ATol) ||
          !isApprox(VertViscH(ICell, 2), ExpectedOneLayer, RTol, ATol) ||
          !(NonLocalH(ICell, 2) > 0.0_Real) ||
          VertDiffH(ICell, 3) != 0.0_Real || VertViscH(ICell, 3) != 0.0_Real ||
          NonLocalH(ICell, 3) != 0.0_Real) {
         ++NumErrors;
      }
   }
   setCoefficientTestGeometry();
   checkResult("coefficient vertical-domain edges", NumErrors);
}

void testCoefficientInvalidWetBounds() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();

   Array2DReal Density("KPPMixTest-InvalidBoundsDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array1DReal UStar("KPPMixTest-InvalidBoundsUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-InvalidBoundsB0", Mesh->NCellsSize);
   deepCopy(Density, -99.0_Real);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, -1.0e-7_Real);
   deepCopy(KPPInstance->VertDiff, -7.0_Real);
   deepCopy(KPPInstance->VertVisc, -7.0_Real);
   deepCopy(KPPInstance->VertNonLocalFlux, -7.0_Real);
   deepCopy(KPPInstance->TurbulentVelocityScale, -7.0_Real);

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   parallelFor(
       "KPPMixTest-SetInvalidWetBounds", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 2;
          MaxLayerCell(ICell) = 1;
       });

   KPPInstance->UseEnhancedDiffusion = false;
   KPPInstance->MatchTechnique       = KPPMatchType::SimpleShapes;
   KPPInstance->computeMixingCoefficients(UStar, B0);

   const auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto VertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   const auto NonLocalH = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   const auto TurbVelH =
       createHostMirrorCopy(KPPInstance->TurbulentVelocityScale);
   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      for (I4 K = 0; K <= NVertLayers; ++K) {
         if (VertDiffH(ICell, K) != -7.0_Real ||
             VertViscH(ICell, K) != -7.0_Real ||
             NonLocalH(ICell, K) != -7.0_Real ||
             TurbVelH(ICell, K) != -7.0_Real) {
            ++NumErrors;
         }
      }
   }

   setCoefficientTestGeometry();
   checkResult("coefficient invalid wet bounds", NumErrors);
}

void testConfiguredValues() {
   Config VertMixConfig("VertMix");
   Config KPPConfig("KPP");
   Error Err;
   Err += Config::getOmegaConfig()->get(VertMixConfig);
   Err += VertMixConfig.get(KPPConfig);

   bool ExpectedEnabled        = false;
   bool ExpectedSmoothing      = false;
   bool ExpectedEnhanced       = false;
   bool ExpectedDebug          = true;
   Real ExpectedCriticalRi     = 0.0_Real;
   Real ExpectedLangmuirIce    = 0.0_Real;
   Real ExpectedMinimumOSBLIce = 0.0_Real;
   Real ExpectedMinimumOSBL    = 0.0_Real;
   std::string ExpectedMatch;
   std::string ExpectedInterp;
   Err += KPPConfig.get("Enable", ExpectedEnabled);
   Err += KPPConfig.get("UseOSBLSmoothing", ExpectedSmoothing);
   Err += KPPConfig.get("UseEnhancedDiffusion", ExpectedEnhanced);
   Err += KPPConfig.get("DebugDiagnostics", ExpectedDebug);
   Err += KPPConfig.get("CriticalBulkRichardsonNumber", ExpectedCriticalRi);
   Err += KPPConfig.get("IceFractionThresholdForLangmuir", ExpectedLangmuirIce);
   Err += KPPConfig.get("IceFractionThresholdForMinimumOSBL",
                        ExpectedMinimumOSBLIce);
   Err += KPPConfig.get("MinimumOSBLUnderSeaIce", ExpectedMinimumOSBL);
   Err += KPPConfig.get("MatchTechnique", ExpectedMatch);
   Err += KPPConfig.get("InterpType2", ExpectedInterp);
   CHECK_ERROR_ABORT(Err, "KPPMixTest: unable to read configured KPP values");

   const KPPMatchType ExpectedMatchType = ExpectedMatch == "MatchBoth"
                                              ? KPPMatchType::MatchBoth
                                              : KPPMatchType::SimpleShapes;

   const KPPMix *KPPInstance = KPPMix::getInstance();
   int NumErrors             = 0;
   if (KPPInstance->Enabled != ExpectedEnabled ||
       KPPInstance->UseOSBLSmoothing != ExpectedSmoothing ||
       KPPInstance->UseEnhancedDiffusion != ExpectedEnhanced ||
       KPPInstance->DebugDiagnostics != ExpectedDebug ||
       KPPInstance->MatchTechnique != ExpectedMatchType ||
       KPPInstance->InterpType2Str != ExpectedInterp ||
       !isApprox(KPPInstance->CriticalRichardson, ExpectedCriticalRi, RTol,
                 ATol) ||
       !isApprox(KPPInstance->IceFractionThresholdForLangmuir,
                 ExpectedLangmuirIce, RTol, ATol) ||
       !isApprox(KPPInstance->IceFractionThresholdForMinimumOSBL,
                 ExpectedMinimumOSBLIce, RTol, ATol) ||
       !isApprox(KPPInstance->MinimumOSBLUnderSeaIce, ExpectedMinimumOSBL, RTol,
                 ATol) ||
       !isApprox(KPPInstance->SurfaceLayerExtent, 0.1_Real, RTol, ATol) ||
       !KPPInstance->UseLangmuirTurbulence) {
      ++NumErrors;
   }
   checkResult("configured values and optional defaults", NumErrors);
}

void testOSBLDepth() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());

   Array2DReal Density("KPPMixTest-BLDDensity", Mesh->NCellsSize, NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-BLDNormalVelocity", Mesh->NEdgesSize,
                              NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-BLDTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-BLDUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-BLDB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-BLDBVF", Mesh->NCellsSize, NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-BLDIce", Mesh->NCellsSize);
   Array1DReal Wind;

   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);
   deepCopy(BVF, 1.0_Real);
   deepCopy(IceFraction, 0.0_Real);

   constexpr Real RiScaling = 1.0_Real - 0.5_Real * KPP::SurfaceLayerExtent;
   constexpr Real TestN     = 1.0_Real;
   const Real UnresolvedShearConstant =
       Kokkos::sqrt(0.2_Real / (KPP::CMoS * KPP::SurfaceLayerExtent)) /
       (VonKar * VonKar);
   const Real WindTurbulentScale = VonKar * 0.02_Real;
   parallelFor(
       "KPPMixTest-SetRichardsonDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          Real TargetRi = 0.0_Real;
          if (K == 1) {
             TargetRi = 0.1_Real;
          } else if (K >= 2) {
             TargetRi = 0.4_Real;
          }
          const Real ZCenter = LayerThickness * (K + 0.5_Real);
          const Real Vt2     = 1.7_Real * UnresolvedShearConstant * ZCenter *
                           TestN * WindTurbulentScale / 0.25_Real;
          const Real DeltaRho =
              TargetRi * Vt2 * RhoSw / (RiScaling * Gravity * ZCenter);
          Density(ICell, K) = RhoSw + DeltaRho;
       });

   KPPInstance->CriticalRichardson    = 0.25_Real;
   KPPInstance->SurfaceLayerExtent    = KPP::SurfaceLayerExtent;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->UseOSBLSmoothing      = false;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);

   auto BLDH      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   auto BLDIndexH = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   auto BulkRiH   = createHostMirrorCopy(KPPInstance->BulkRichardsonNumber);
   const auto BulkShearH =
       createHostMirrorCopy(KPPInstance->BulkRichardsonShear);
   const auto UnresolvedShearH =
       createHostMirrorCopy(KPPInstance->UnresolvedShear);
   const auto BuoyancyJumpH = createHostMirrorCopy(KPPInstance->BuoyancyJump);

   constexpr Real RiAbove   = 0.1_Real;
   constexpr Real RiBelow   = 0.4_Real;
   constexpr Real ZAbove    = 15.0_Real;
   constexpr Real Slope     = 0.01_Real;
   constexpr Real Quadratic = 0.002_Real;
   const Real Discriminant =
       Slope * Slope - 4.0_Real * Quadratic * (RiAbove - 0.25_Real);
   const Real ExpectedBLD =
       ZAbove + (-Slope + Kokkos::sqrt(Discriminant)) / (2.0_Real * Quadratic);
   const Real ExpectedVt2 = 1.7_Real * UnresolvedShearConstant * 25.0_Real *
                            TestN * WindTurbulentScale / 0.25_Real;
   const Real ExpectedDeltaB = 0.4_Real * ExpectedVt2 / (RiScaling * 25.0_Real);

   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BulkRiH(ICell, 2), RiAbove, BLDRTol, ATol) ||
          !isApprox(BulkRiH(ICell, 3), RiBelow, BLDRTol, ATol) ||
          !isApprox(BLDH(ICell), ExpectedBLD, BLDRTol, ATol) ||
          BLDIndexH(ICell) != 2 ||
          !isApprox(BulkShearH(ICell, 3), 0.0_Real, RTol, ATol) ||
          !isApprox(UnresolvedShearH(ICell, 3), ExpectedVt2, RTol, ATol) ||
          !isApprox(BuoyancyJumpH(ICell, 3), ExpectedDeltaB, BLDRTol, ATol)) {
         ++NumErrors;
      }
   }
   checkResult("analytic boundary-layer depth", NumErrors);

   parallelFor(
       "KPPMixTest-SetLinearRichardsonDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const Real TargetRi = K == 0   ? 0.0_Real
                                : K == 1 ? 0.1_Real
                                : K == 2 ? 0.2_Real
                                         : 0.3_Real;
          const Real ZCenter  = LayerThickness * (K + 0.5_Real);
          const Real Vt2      = 1.7_Real * UnresolvedShearConstant * ZCenter *
                           TestN * WindTurbulentScale / 0.25_Real;
          const Real DeltaRho =
              TargetRi * Vt2 * RhoSw / (RiScaling * Gravity * ZCenter);
          Density(ICell, K) = RhoSw + DeltaRho;
       });
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   BLDH      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   BLDIndexH = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   BulkRiH   = createHostMirrorCopy(KPPInstance->BulkRichardsonNumber);
   NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BulkRiH(ICell, 3), 0.2_Real, BLDRTol, ATol) ||
          !isApprox(BulkRiH(ICell, 4), 0.3_Real, BLDRTol, ATol) ||
          !isApprox(BLDH(ICell), 30.0_Real, BLDRTol, ATol) ||
          BLDIndexH(ICell) != 2) {
         ++NumErrors;
      }
   }
   checkResult("boundary-layer linear interpolation", NumErrors);

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   parallelFor(
       "KPPMixTest-SetPartialBLDBounds", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 2;
          MaxLayerCell(ICell) = 4;
       });
   parallelFor(
       "KPPMixTest-SetPartialBLDColumn", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          if (K < NVertLayers) {
             Density(ICell, K) = K >= 2 && K <= 4 ? RhoSw : -99.0_Real;
          }
          BVF(ICell, K) = K >= 3 && K <= 5 ? 1.0_Real : -99.0_Real;
       });
   parallelFor(
       "KPPMixTest-SetPartialBLDEdges", {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          const Real Value         = K >= 2 && K <= 4 ? 0.0_Real : 99.0_Real;
          NormalVelocity(IEdge, K) = Value;
          TangentialVelocity(IEdge, K) = Value;
       });
   VCoord->minMaxLayerEdge(Halo::getDefault());
   deepCopy(KPPInstance->BulkRichardsonNumber, FillValueReal);
   deepCopy(KPPInstance->BulkRichardsonShear, FillValueReal);
   deepCopy(KPPInstance->UnresolvedShear, FillValueReal);
   deepCopy(KPPInstance->BuoyancyJump, FillValueReal);
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   BLDH      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   BLDIndexH = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   BulkRiH   = createHostMirrorCopy(KPPInstance->BulkRichardsonNumber);
   NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BLDH(ICell), 45.0_Real, RTol, ATol) ||
          BLDIndexH(ICell) != 4 || BulkRiH(ICell, 1) != FillValueReal ||
          BulkRiH(ICell, 2) != 0.0_Real) {
         ++NumErrors;
      }
   }
   checkResult("boundary-layer partial wet column", NumErrors);
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());
   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);
   deepCopy(BVF, 1.0_Real);

   parallelFor(
       "KPPMixTest-RestoreRichardsonDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          Real TargetRi = 0.0_Real;
          if (K == 1) {
             TargetRi = 0.1_Real;
          } else if (K >= 2) {
             TargetRi = 0.4_Real;
          }
          const Real ZCenter = LayerThickness * (K + 0.5_Real);
          const Real Vt2     = 1.7_Real * UnresolvedShearConstant * ZCenter *
                           TestN * WindTurbulentScale / 0.25_Real;
          const Real DeltaRho =
              TargetRi * Vt2 * RhoSw / (RiScaling * Gravity * ZCenter);
          Density(ICell, K) = RhoSw + DeltaRho;
       });

   parallelFor(
       "KPPMixTest-SetResolvedShear", {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelocity(IEdge, K)     = 0.1_Real * K;
          TangentialVelocity(IEdge, K) = 0.2_Real * K;
       });
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   BLDH = createHostMirrorCopy(KPPInstance->OSBLDepth);
   const auto ShearedBulkRiH =
       createHostMirrorCopy(KPPInstance->BulkRichardsonNumber);
   const auto ResolvedShearH =
       createHostMirrorCopy(KPPInstance->BulkRichardsonShear);
   NumErrors                            = 0;
   constexpr Real ExpectedResolvedShear = 0.2_Real;
   const Real ExpectedShearedRi = RiScaling * ExpectedDeltaB * 25.0_Real /
                                  (ExpectedResolvedShear + ExpectedVt2);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(ResolvedShearH(ICell, 3), ExpectedResolvedShear, RTol,
                    ATol) ||
          !isApprox(ShearedBulkRiH(ICell, 3), ExpectedShearedRi, BLDRTol,
                    ATol) ||
          !(BLDH(ICell) > ExpectedBLD)) {
         ++NumErrors;
      }
   }
   checkResult("boundary-layer resolved shear", NumErrors);
   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);

   deepCopy(Density, RhoSw);
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   BLDH      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   BLDIndexH = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   NumErrors = 0;
   const Real DeepestMidpoint =
       LayerThickness * (static_cast<Real>(NVertLayers) - 0.5_Real);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BLDH(ICell), DeepestMidpoint, RTol, ATol) ||
          BLDIndexH(ICell) != NVertLayers - 1) {
         ++NumErrors;
      }
   }
   checkResult("boundary-layer no-crossing fallback", NumErrors);

   parallelFor(
       "KPPMixTest-SetShallowCrossingDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const Real ZCenter  = LayerThickness * (K + 0.5_Real);
          const Real TargetRi = K == 0 ? 0.0_Real : 1.0_Real;
          const Real Vt2      = 1.7_Real * UnresolvedShearConstant * ZCenter *
                           TestN * WindTurbulentScale / 0.25_Real;
          const Real DeltaRho =
              TargetRi * Vt2 * RhoSw / (RiScaling * Gravity * ZCenter);
          Density(ICell, K) = RhoSw + DeltaRho;
       });
   deepCopy(IceFraction, 0.5_Real);
   KPPInstance->IceFractionThresholdForMinimumOSBL = 0.15_Real;
   KPPInstance->MinimumOSBLUnderSeaIce             = 25.0_Real;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   BLDH      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   BLDIndexH = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BLDH(ICell), 25.0_Real, RTol, ATol) ||
          BLDIndexH(ICell) != 2) {
         ++NumErrors;
      }
   }
   checkResult("boundary-layer sea-ice minimum", NumErrors);
}

void testBoundaryLayerNonuniformThickness() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;

   Array2DReal Density("KPPMixTest-NonuniformDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-NonuniformNormalVelocity",
                              Mesh->NEdgesSize, NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-NonuniformTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-NonuniformUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-NonuniformB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-NonuniformBVF", Mesh->NCellsSize,
                   NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-NonuniformIce", Mesh->NCellsSize);
   Array1DReal Wind;

   OMEGA_SCOPE(GeomZInterface, VCoord->GeomZInterface);
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(SshCell, VCoord->SshCell);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   constexpr Real RiScaling = 1.0_Real - 0.5_Real * KPP::SurfaceLayerExtent;
   constexpr Real TestN     = 1.0_Real;
   constexpr Real TestUStar = 0.02_Real;
   const Real UnresolvedShearConstant =
       Kokkos::sqrt(0.2_Real / (KPP::CMoS * KPP::SurfaceLayerExtent)) /
       (VonKar * VonKar);
   const Real WindTurbulentScale = VonKar * TestUStar;

   deepCopy(Density, RhoSw);
   deepCopy(UStar, TestUStar);
   deepCopy(B0, 0.0_Real);
   deepCopy(BVF, TestN * TestN);
   deepCopy(IceFraction, 0.0_Real);
   parallelFor(
       "KPPMixTest-SetNonuniformVelocity", {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelocity(IEdge, K)     = K == 0   ? 0.1_Real
                                         : K == 1 ? 0.4_Real
                                         : K == 2 ? 0.1_Real
                                         : K == 3 ? 0.6_Real
                                                  : 0.6_Real;
          TangentialVelocity(IEdge, K) = K == 0   ? 0.2_Real
                                         : K == 1 ? 0.5_Real
                                         : K == 2 ? 0.2_Real
                                         : K == 3 ? 0.8_Real
                                                  : 0.8_Real;
       });

   parallelFor(
       "KPPMixTest-SetNonuniformColumn", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 0;
          MaxLayerCell(ICell) = NVertLayers - 1;
          SshCell(ICell)      = 0.0_Real;

          for (I4 K = 0; K <= NVertLayers; ++K) {
             Real Depth = 25.0_Real + 25.0_Real * (K - 4);
             if (K == 0) {
                Depth = 0.0_Real;
             } else if (K == 1) {
                Depth = 1.0_Real;
             } else if (K == 2) {
                Depth = 3.0_Real;
             } else if (K == 3) {
                Depth = 10.0_Real;
             } else if (K == 4) {
                Depth = 25.0_Real;
             }
             GeomZInterface(ICell, K) = -Depth;
             if (K < NVertLayers) {
                Real NextDepth = 25.0_Real + 25.0_Real * (K - 3);
                if (K == 0) {
                   NextDepth = 1.0_Real;
                } else if (K == 1) {
                   NextDepth = 3.0_Real;
                } else if (K == 2) {
                   NextDepth = 10.0_Real;
                } else if (K == 3) {
                   NextDepth = 25.0_Real;
                }
                GeomZMid(ICell, K) = -0.5_Real * (Depth + NextDepth);
             }
          }

          constexpr Real Shear1 = 0.18_Real;
          constexpr Real Shear3 = 0.25_Real;
          const Real Vt2Layer1 = 1.7_Real * UnresolvedShearConstant * 2.0_Real *
                                 TestN * WindTurbulentScale / 0.25_Real;
          const Real Vt2Layer2 = 1.7_Real * UnresolvedShearConstant * 6.5_Real *
                                 TestN * WindTurbulentScale / 0.25_Real;
          const Real Vt2Layer3 = 1.7_Real * UnresolvedShearConstant *
                                 17.5_Real * TestN * WindTurbulentScale /
                                 0.25_Real;
          const Real DeltaRho1 = 0.05_Real * (Shear1 + Vt2Layer1) * RhoSw /
                                 (RiScaling * Gravity * 2.0_Real);
          const Real DeltaRho2 =
              0.10_Real * Vt2Layer2 * RhoSw / (RiScaling * Gravity * 6.5_Real);
          const Real DeltaRho3 = 0.40_Real * (Shear3 + Vt2Layer3) * RhoSw /
                                 (RiScaling * Gravity * 17.5_Real);
          Density(ICell, 0) = RhoSw;
          Density(ICell, 1) = RhoSw + DeltaRho1;
          Density(ICell, 2) = RhoSw + DeltaRho2;

          // At k=3, the 2.5 m surface layer contains the unequal 1 m and
          // 2 m layers. Construct rho(3) relative to that weighted mean.
          const Real WeightedSurfaceDensity =
              (Density(ICell, 0) + 2.0_Real * Density(ICell, 1)) / 3.0_Real;
          Density(ICell, 3) = WeightedSurfaceDensity + DeltaRho3;
       });
   VCoord->minMaxLayerEdge(Halo::getDefault());

   KPPInstance->CriticalRichardson    = 0.25_Real;
   KPPInstance->SurfaceLayerExtent    = KPP::SurfaceLayerExtent;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->UseOSBLSmoothing      = false;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);

   const auto BLDH      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   const auto BLDIndexH = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   const auto BulkRiH = createHostMirrorCopy(KPPInstance->BulkRichardsonNumber);
   const auto BulkShearH =
       createHostMirrorCopy(KPPInstance->BulkRichardsonShear);
   const auto BuoyancyJumpH = createHostMirrorCopy(KPPInstance->BuoyancyJump);

   constexpr Real ZPrevious  = 2.0_Real;
   constexpr Real ZAbove     = 6.5_Real;
   constexpr Real ZBelow     = 17.5_Real;
   constexpr Real RiPrevious = 0.05_Real;
   constexpr Real RiAbove    = 0.10_Real;
   constexpr Real RiBelow    = 0.40_Real;
   const Real Slope          = (RiAbove - RiPrevious) / (ZAbove - ZPrevious);
   const Real H              = ZBelow - ZAbove;
   const Real Quadratic      = (RiBelow - RiAbove - Slope * H) / (H * H);
   const Real Discriminant =
       Slope * Slope - 4.0_Real * Quadratic * (RiAbove - 0.25_Real);
   const Real ExpectedBLD =
       ZAbove + (-Slope + Kokkos::sqrt(Discriminant)) / (2.0_Real * Quadratic);
   const Real ExpectedDeltaB =
       0.40_Real *
       (0.25_Real + 1.7_Real * UnresolvedShearConstant * ZBelow * TestN *
                        WindTurbulentScale / 0.25_Real) /
       (RiScaling * ZBelow);

   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BulkRiH(ICell, 2), RiPrevious, BLDRTol, ATol) ||
          !isApprox(BulkRiH(ICell, 3), RiAbove, BLDRTol, ATol) ||
          !isApprox(BulkRiH(ICell, 4), RiBelow, BLDRTol, ATol) ||
          !isApprox(BulkShearH(ICell, 2), 0.18_Real, RTol, ATol) ||
          !isApprox(BulkShearH(ICell, 4), 0.25_Real, RTol, ATol) ||
          !isApprox(BuoyancyJumpH(ICell, 4), ExpectedDeltaB, BLDRTol, ATol) ||
          !isApprox(BLDH(ICell), ExpectedBLD, BLDRTol, ATol) ||
          BLDIndexH(ICell) != 3) {
         ++NumErrors;
      }
   }

   checkResult("boundary-layer nonuniform thickness", NumErrors);
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());
}

void testBoundaryLayerHorizontalThicknessVariation() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;

   Array2DReal Density("KPPMixTest-HorizontalDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-HorizontalNormalVelocity",
                              Mesh->NEdgesSize, NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-HorizontalTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-HorizontalUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-HorizontalB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-HorizontalBVF", Mesh->NCellsSize,
                   NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-HorizontalIce", Mesh->NCellsSize);
   Array1DReal Wind;

   OMEGA_SCOPE(GeomZInterface, VCoord->GeomZInterface);
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(SshCell, VCoord->SshCell);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);
   deepCopy(BVF, 1.0_Real);
   deepCopy(IceFraction, 0.0_Real);
   parallelFor(
       "KPPMixTest-SetHorizontalVelocity", {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelocity(IEdge, K)     = K == 0   ? 0.1_Real
                                         : K == 1 ? 0.4_Real
                                                  : 0.6_Real;
          TangentialVelocity(IEdge, K) = K == 0   ? 0.2_Real
                                         : K == 1 ? 0.5_Real
                                                  : 0.8_Real;
       });

   parallelFor(
       "KPPMixTest-SetHorizontalColumn", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 0;
          MaxLayerCell(ICell) = NVertLayers - 1;

          const Real Ssh         = 0.1_Real * (ICell % 2);
          const Real DepthOffset = 0.25_Real * (ICell % 3);
          SshCell(ICell)         = Ssh;
          for (I4 K = 0; K <= NVertLayers; ++K) {
             Real Depth = 25.0_Real + 25.0_Real * (K - 4);
             if (K == 0) {
                Depth = 0.0_Real;
             } else if (K == 1) {
                Depth = 1.0_Real;
             } else if (K == 2) {
                Depth = 3.0_Real + DepthOffset;
             } else if (K == 3) {
                Depth = 10.0_Real;
             } else if (K == 4) {
                Depth = 25.0_Real;
             }
             GeomZInterface(ICell, K) = Ssh - Depth;
             if (K < NVertLayers) {
                Real NextDepth = 25.0_Real + 25.0_Real * (K - 3);
                if (K == 0) {
                   NextDepth = 1.0_Real;
                } else if (K == 1) {
                   NextDepth = 3.0_Real + DepthOffset;
                } else if (K == 2) {
                   NextDepth = 10.0_Real;
                } else if (K == 3) {
                   NextDepth = 25.0_Real;
                }
                GeomZMid(ICell, K) = Ssh - 0.5_Real * (Depth + NextDepth);
             }
          }
       });
   VCoord->minMaxLayerEdge(Halo::getDefault());

   KPPInstance->CriticalRichardson    = 0.25_Real;
   KPPInstance->SurfaceLayerExtent    = KPP::SurfaceLayerExtent;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->UseOSBLSmoothing      = false;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);

   const auto BulkShearH =
       createHostMirrorCopy(KPPInstance->BulkRichardsonShear);

   int NumErrors              = 0;
   int NumDiscriminatingCells = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      Real WeightedShear = 0.0_Real;
      Real WeightSum     = 0.0_Real;
      for (I4 J = 0; J < Mesh->NEdgesOnCellH(ICell); ++J) {
         const I4 JCell = Mesh->CellsOnCellH(ICell, J);
         // CellsOnCell uses NCellsAll as a sentinel for missing/boundary
         // neighbors; KPP excludes these edges, so the reference must too.
         if (JCell >= Mesh->NCellsAll) {
            continue;
         }
         const I4 IEdge        = Mesh->EdgesOnCellH(ICell, J);
         const Real OffsetCell = 0.25_Real * (ICell % 3);
         const Real OffsetNbr  = 0.25_Real * (JCell % 3);
         const Real Thick1     = 2.0_Real + 0.5_Real * (OffsetCell + OffsetNbr);
         const Real SumThickness = 1.0_Real + Thick1;
         const Real UnAvg = (0.1_Real + 0.4_Real * Thick1) / SumThickness;
         const Real VtAvg = (0.2_Real + 0.5_Real * Thick1) / SumThickness;
         const Real DUn   = 0.6_Real - UnAvg;
         const Real DVt   = 0.8_Real - VtAvg;
         const Real Weight =
             0.25_Real * Mesh->DcEdgeH(IEdge) * Mesh->DvEdgeH(IEdge);
         WeightedShear += Weight * (DUn * DUn + DVt * DVt);
         WeightSum += Weight;
      }
      const Real ExpectedShear = WeightedShear / WeightSum;

      const Real OldThick1       = 2.0_Real + 0.25_Real * (ICell % 3);
      const Real OldSumThickness = 1.0_Real + OldThick1;
      const Real OldUnAvg = (0.1_Real + 0.4_Real * OldThick1) / OldSumThickness;
      const Real OldVtAvg = (0.2_Real + 0.5_Real * OldThick1) / OldSumThickness;
      const Real OldDUn   = 0.6_Real - OldUnAvg;
      const Real OldDVt   = 0.8_Real - OldVtAvg;
      const Real OldShear = OldDUn * OldDUn + OldDVt * OldDVt;
      if (!isApprox(ExpectedShear, OldShear, RTol, ATol)) {
         ++NumDiscriminatingCells;
      }

      if (!isApprox(BulkShearH(ICell, 4), ExpectedShear, RTol, ATol)) {
         ++NumErrors;
      }
   }
   if (NumDiscriminatingCells == 0) {
      ++NumErrors;
   }

   checkResult("boundary-layer horizontal thickness variation", NumErrors);
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());
}

// KPP depths are measured below the free surface, so rigidly translating the
// whole column by the sea surface height must leave every KPP output unchanged.
void testSshOffsetInvariance() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;

   Array2DReal Density("KPPMixTest-SshDensity", Mesh->NCellsSize, NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-SshNormalVelocity", Mesh->NEdgesSize,
                              NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-SshTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-SshUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-SshB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-SshBVF", Mesh->NCellsSize, NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-SshIce", Mesh->NCellsSize);
   Array1DReal Wind;

   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, -1.0e-7_Real);
   deepCopy(BVF, 1.0_Real);
   deepCopy(IceFraction, 0.0_Real);

   constexpr Real RiScaling = 1.0_Real - 0.5_Real * KPP::SurfaceLayerExtent;
   constexpr Real TestN     = 1.0_Real;
   const Real UnresolvedShearConstant =
       Kokkos::sqrt(0.2_Real / (KPP::CMoS * KPP::SurfaceLayerExtent)) /
       (VonKar * VonKar);
   const Real WindTurbulentScale = VonKar * 0.02_Real;
   parallelFor(
       "KPPMixTest-SetSshDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          Real TargetRi = 0.0_Real;
          if (K == 1) {
             TargetRi = 0.1_Real;
          } else if (K >= 2) {
             TargetRi = 0.4_Real;
          }
          const Real ZCenter = LayerThickness * (K + 0.5_Real);
          const Real Vt2     = 1.7_Real * UnresolvedShearConstant * ZCenter *
                           TestN * WindTurbulentScale / 0.25_Real;
          const Real DeltaRho =
              TargetRi * Vt2 * RhoSw / (RiScaling * Gravity * ZCenter);
          Density(ICell, K) = RhoSw + DeltaRho;
       });

   KPPInstance->CriticalRichardson    = 0.25_Real;
   KPPInstance->SurfaceLayerExtent    = KPP::SurfaceLayerExtent;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->UseOSBLSmoothing      = false;
   KPPInstance->UseEnhancedDiffusion  = true;
   KPPInstance->MatchTechnique        = KPPMatchType::SimpleShapes;

   auto runWithSsh = [&](Real Ssh) {
      setCoefficientTestGeometry(Ssh);
      VCoord->minMaxLayerEdge(Halo::getDefault());
      KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                    UStar, B0, BVF, IceFraction, Wind);
      KPPInstance->computeMixingCoefficients(UStar, B0);
   };

   runWithSsh(0.0_Real);
   const auto BLDRef      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   const auto BLDIndexRef = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   const auto DiffRef     = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto ViscRef     = createHostMirrorCopy(KPPInstance->VertVisc);
   const auto NonLocalRef = createHostMirrorCopy(KPPInstance->VertNonLocalFlux);

   // Large enough that an uncorrected geoid-referenced depth shifts the OSBL
   // search across the analytic crossing.
   constexpr Real TestSsh = 0.75_Real;
   runWithSsh(TestSsh);
   const auto BLDShift      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   const auto BLDIndexShift = createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   const auto DiffShift     = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto ViscShift     = createHostMirrorCopy(KPPInstance->VertVisc);
   const auto NonLocalShift =
       createHostMirrorCopy(KPPInstance->VertNonLocalFlux);

   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      // Guard against a degenerate all-zero comparison.
      if (!(BLDRef(ICell) > 0.0_Real)) {
         ++NumErrors;
         continue;
      }
      if (!isApprox(BLDShift(ICell), BLDRef(ICell), BLDRTol, ATol) ||
          BLDIndexShift(ICell) != BLDIndexRef(ICell)) {
         ++NumErrors;
         continue;
      }
      for (I4 K = 0; K <= NVertLayers; ++K) {
         if (!isApprox(DiffShift(ICell, K), DiffRef(ICell, K), BLDRTol, ATol) ||
             !isApprox(ViscShift(ICell, K), ViscRef(ICell, K), BLDRTol, ATol) ||
             !isApprox(NonLocalShift(ICell, K), NonLocalRef(ICell, K), BLDRTol,
                       ATol)) {
            ++NumErrors;
            break;
         }
      }
   }

   checkResult("sea surface height offset invariance", NumErrors);
   KPPInstance->UseOSBLSmoothing = true;
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());
}

void testBoundaryLayerEdgeFallbacks() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());

   Array2DReal Density("KPPMixTest-EdgeFallbackDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-EdgeFallbackNormalVelocity",
                              Mesh->NEdgesSize, NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-EdgeFallbackTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-EdgeFallbackUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-EdgeFallbackB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-EdgeFallbackBVF", Mesh->NCellsSize,
                   NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-EdgeFallbackIce", Mesh->NCellsSize);
   Array1DReal Wind;
   Array1DReal OriginalDcEdge("KPPMixTest-OriginalDcEdge",
                              Mesh->DcEdge.extent(0));
   deepCopy(OriginalDcEdge, Mesh->DcEdge);

   deepCopy(Density, RhoSw);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, 0.0_Real);
   deepCopy(BVF, 1.0_Real);
   deepCopy(IceFraction, 0.0_Real);
   parallelFor(
       "KPPMixTest-SetEdgeFallbackVelocity", {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelocity(IEdge, K)     = 0.1_Real * K;
          TangentialVelocity(IEdge, K) = 0.2_Real * K;
       });

   KPPInstance->CriticalRichardson    = 0.25_Real;
   KPPInstance->SurfaceLayerExtent    = KPP::SurfaceLayerExtent;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->UseOSBLSmoothing      = false;

   // Zero geometric weights leave no resolved-shear contribution.
   deepCopy(Mesh->DcEdge, 0.0_Real);
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   auto BulkShearH = createHostMirrorCopy(KPPInstance->BulkRichardsonShear);
   int NumErrors   = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BulkShearH(ICell, 3), 0.0_Real, RTol, ATol)) {
         ++NumErrors;
      }
   }
   deepCopy(Mesh->DcEdge, OriginalDcEdge);

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   deepCopy(MinLayerEdgeBot, -1);
   deepCopy(MaxLayerEdgeTop, -1);
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   BulkShearH = createHostMirrorCopy(KPPInstance->BulkRichardsonShear);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(BulkShearH(ICell, 3), 0.0_Real, RTol, ATol)) {
         ++NumErrors;
      }
   }

   VCoord->minMaxLayerEdge(Halo::getDefault());
   checkResult("boundary-layer edge fallbacks", NumErrors);
}

void testBoundaryLayerLangmuir() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());

   Array2DReal Density("KPPMixTest-LangmuirDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-LangmuirNormalVelocity",
                              Mesh->NEdgesSize, NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-LangmuirTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-LangmuirUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-LangmuirB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-LangmuirBVF", Mesh->NCellsSize, NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-LangmuirIce", Mesh->NCellsSize);
   Array1DReal Wind("KPPMixTest-LangmuirWind", Mesh->NCellsSize);

   constexpr Real TestUStar = 0.02_Real;
   constexpr Real TestB0    = -1.0e-7_Real;
   constexpr Real TestN     = 1.0_Real;
   constexpr Real RiScaling = 1.0_Real - 0.5_Real * KPP::SurfaceLayerExtent;
   const Real UnresolvedShearConstant =
       Kokkos::sqrt(0.2_Real / (KPP::CMoS * KPP::SurfaceLayerExtent)) /
       (VonKar * VonKar);

   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);
   deepCopy(UStar, TestUStar);
   deepCopy(B0, TestB0);
   deepCopy(BVF, TestN * TestN);
   deepCopy(IceFraction, 0.0_Real);
   deepCopy(Wind, 10.0_Real);
   parallelFor(
       "KPPMixTest-SetLangmuirDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const Real ZDepth  = LayerThickness * (K + 1.0_Real);
          const Real ZCenter = LayerThickness * (K + 0.5_Real);
          const Real Zeta = KPP::SurfaceLayerExtent * ZDepth * VonKar * TestB0 /
                            (TestUStar * TestUStar * TestUStar);
          const Real PhiInv = Kokkos::sqrt(1.0_Real - 16.0_Real * Zeta);
          const Real WTurb  = VonKar * TestUStar * PhiInv;
          const Real Vt2    = 1.7_Real * UnresolvedShearConstant * ZCenter *
                           TestN * WTurb / 0.25_Real;
          const Real TargetRi =
              K == 0 ? 0.0_Real : (K == 1 ? 0.1_Real : 0.26_Real);
          const Real DeltaRho =
              TargetRi * Vt2 * RhoSw / (RiScaling * Gravity * ZCenter);
          Density(ICell, K) = RhoSw + DeltaRho;
       });

   KPPInstance->CriticalRichardson    = 0.25_Real;
   KPPInstance->SurfaceLayerExtent    = KPP::SurfaceLayerExtent;
   KPPInstance->UseOSBLSmoothing      = false;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   Array1DReal DisabledBLD("KPPMixTest-LangmuirDisabledBLD", Mesh->NCellsSize);
   Array2DReal DisabledRi("KPPMixTest-LangmuirDisabledRi", Mesh->NCellsSize,
                          NVertLayers + 1);
   Array2DReal DisabledVt2("KPPMixTest-LangmuirDisabledVt2", Mesh->NCellsSize,
                           NVertLayers + 1);
   deepCopy(DisabledBLD, KPPInstance->OSBLDepth);
   deepCopy(DisabledRi, KPPInstance->BulkRichardsonNumber);
   deepCopy(DisabledVt2, KPPInstance->UnresolvedShear);
   const auto DisabledBLDH = createHostMirrorCopy(DisabledBLD);
   const auto DisabledRiH  = createHostMirrorCopy(DisabledRi);
   const auto DisabledVt2H = createHostMirrorCopy(DisabledVt2);

   KPPInstance->UseLangmuirTurbulence           = true;
   KPPInstance->IceFractionThresholdForLangmuir = 0.05_Real;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   Array1DReal EnabledBLD("KPPMixTest-LangmuirEnabledBLD", Mesh->NCellsSize);
   Array2DReal EnabledRi("KPPMixTest-LangmuirEnabledRi", Mesh->NCellsSize,
                         NVertLayers + 1);
   Array2DReal EnabledVt2("KPPMixTest-LangmuirEnabledVt2", Mesh->NCellsSize,
                          NVertLayers + 1);
   deepCopy(EnabledBLD, KPPInstance->OSBLDepth);
   deepCopy(EnabledRi, KPPInstance->BulkRichardsonNumber);
   deepCopy(EnabledVt2, KPPInstance->UnresolvedShear);
   const auto EnabledBLDH = createHostMirrorCopy(EnabledBLD);
   const auto EnabledRiH  = createHostMirrorCopy(EnabledRi);
   const auto EnabledVt2H = createHostMirrorCopy(EnabledVt2);

   deepCopy(IceFraction, 0.1_Real);
   KPPInstance->IceFractionThresholdForMinimumOSBL = 2.0_Real;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   const auto SuppressedBLDH = createHostMirrorCopy(KPPInstance->OSBLDepth);
   const auto SuppressedRiH =
       createHostMirrorCopy(KPPInstance->BulkRichardsonNumber);

   constexpr Real ZDepth   = 30.0_Real;
   constexpr Real ZCenter  = 25.0_Real;
   const Real Enhancement  = Kokkos::sqrt(3.0_Real);
   const Real DisabledZeta = KPP::SurfaceLayerExtent * ZDepth * VonKar *
                             TestB0 / (TestUStar * TestUStar * TestUStar);
   const Real EnabledZeta = DisabledZeta * Enhancement;
   const Real DisabledWTurb =
       VonKar * TestUStar * Kokkos::sqrt(1.0_Real - 16.0_Real * DisabledZeta);
   const Real EnabledWTurb =
       VonKar * TestUStar * Kokkos::sqrt(1.0_Real - 16.0_Real * EnabledZeta);
   const Real ExpectedDisabledVt2 = 1.7_Real * UnresolvedShearConstant *
                                    ZCenter * TestN * DisabledWTurb / 0.25_Real;
   const Real ExpectedEnabledVt2 = 1.7_Real * UnresolvedShearConstant *
                                   ZCenter * TestN * EnabledWTurb / 0.25_Real;
   const Real ExpectedEnabledRi =
       0.26_Real * ExpectedDisabledVt2 / ExpectedEnabledVt2;

   int RiErrors          = 0;
   int Vt2Errors         = 0;
   int DepthErrors       = 0;
   int SuppressionErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(DisabledRiH(ICell, 3), 0.26_Real, BLDRTol, ATol) ||
          !isApprox(EnabledRiH(ICell, 3), ExpectedEnabledRi, BLDRTol, ATol)) {
         ++RiErrors;
      }
      if (!isApprox(DisabledVt2H(ICell, 3), ExpectedDisabledVt2, RTol, ATol) ||
          !isApprox(EnabledVt2H(ICell, 3), ExpectedEnabledVt2, RTol, ATol)) {
         ++Vt2Errors;
      }
      if (!(EnabledBLDH(ICell) > DisabledBLDH(ICell))) {
         ++DepthErrors;
      }
      if (!isApprox(SuppressedRiH(ICell, 3), DisabledRiH(ICell, 3), RTol,
                    ATol) ||
          !isApprox(SuppressedBLDH(ICell), DisabledBLDH(ICell), RTol, ATol)) {
         ++SuppressionErrors;
      }
   }
   if (RiErrors != 0 || Vt2Errors != 0 || DepthErrors != 0 ||
       SuppressionErrors != 0) {
      LOG_ERROR(
          "Langmuir BLD failures: Ri={} Vt2={} depth={} suppression={}; "
          "cell 0 disabled Ri={} Vt2={} BLD={}, enabled Ri={} Vt2={} "
          "BLD={}, suppressed Ri={} BLD={}; expected enabled Ri={} Vt2={}",
          RiErrors, Vt2Errors, DepthErrors, SuppressionErrors,
          DisabledRiH(0, 3), DisabledVt2H(0, 3), DisabledBLDH(0),
          EnabledRiH(0, 3), EnabledVt2H(0, 3), EnabledBLDH(0),
          SuppressedRiH(0, 3), SuppressedBLDH(0), ExpectedEnabledRi,
          ExpectedEnabledVt2);
   }
   KPPInstance->UseLangmuirTurbulence = false;
   checkResult("boundary-layer Langmuir enhancement and ice suppression",
               RiErrors + Vt2Errors + DepthErrors + SuppressionErrors);
}

void testBoundaryLayerSmoothing() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());

   Array2DReal Density("KPPMixTest-SmoothingDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-SmoothingNormalVelocity",
                              Mesh->NEdgesSize, NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-SmoothingTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-SmoothingUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-SmoothingB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-SmoothingBVF", Mesh->NCellsSize,
                   NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-SmoothingIce", Mesh->NCellsSize);
   Array1DReal Wind;

   constexpr Real TestUStar = 0.02_Real;
   constexpr Real TestN     = 1.0_Real;
   constexpr Real RiScaling = 1.0_Real - 0.5_Real * KPP::SurfaceLayerExtent;
   const Real UnresolvedShearConstant =
       Kokkos::sqrt(0.2_Real / (KPP::CMoS * KPP::SurfaceLayerExtent)) /
       (VonKar * VonKar);
   const Real WTurb = VonKar * TestUStar;

   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);
   deepCopy(UStar, TestUStar);
   deepCopy(B0, 0.0_Real);
   deepCopy(BVF, TestN * TestN);
   deepCopy(IceFraction, 0.0_Real);
   parallelFor(
       "KPPMixTest-SetSmoothingDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          Real TargetRi = 0.0_Real;
          if (ICell % 2 == 0) {
             TargetRi = K == 1 ? 0.1_Real : (K >= 2 ? 0.4_Real : 0.0_Real);
          } else {
             TargetRi =
                 K == 1 ? 0.05_Real
                        : (K == 2 ? 0.1_Real : (K >= 3 ? 0.4_Real : 0.0_Real));
          }
          const Real ZCenter = LayerThickness * (K + 0.5_Real);
          const Real Vt2     = 1.7_Real * UnresolvedShearConstant * ZCenter *
                           TestN * WTurb / 0.25_Real;
          const Real DeltaRho =
              TargetRi * Vt2 * RhoSw / (RiScaling * Gravity * ZCenter);
          Density(ICell, K) = RhoSw + DeltaRho;
       });

   KPPInstance->CriticalRichardson    = 0.25_Real;
   KPPInstance->SurfaceLayerExtent    = KPP::SurfaceLayerExtent;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->UseOSBLSmoothing      = false;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   Array1DReal UnsmoothedBLD("KPPMixTest-UnsmoothedBLD", Mesh->NCellsSize);
   deepCopy(UnsmoothedBLD, KPPInstance->OSBLDepth);
   const auto UnsmoothedBLDH = createHostMirrorCopy(UnsmoothedBLD);

   KPPInstance->UseOSBLSmoothing = true;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   const auto SmoothedBLDH = createHostMirrorCopy(KPPInstance->OSBLDepth);
   const auto SmoothedIndexH =
       createHostMirrorCopy(KPPInstance->OSBLDepthIndex);

   int NumErrors       = 0;
   int NumChanged      = 0;
   const Real MinDepth = 0.5_Real * LayerThickness;
   const Real MaxDepth =
       LayerThickness * (static_cast<Real>(NVertLayers) - 0.5_Real);
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      Real WeightedDepth = 0.0_Real;
      Real AreaSum       = 0.0_Real;
      I4 ValidNeighbors  = 0;
      for (I4 J = 0; J < Mesh->NEdgesOnCellH(ICell); ++J) {
         const I4 INeighbor = Mesh->CellsOnCellH(ICell, J);
         if (INeighbor == Mesh->NCellsAll) {
            continue;
         }
         WeightedDepth +=
             2.0_Real * Mesh->AreaCellH(INeighbor) * UnsmoothedBLDH(INeighbor);
         AreaSum += 2.0_Real * Mesh->AreaCellH(INeighbor);
         ++ValidNeighbors;
      }
      if (ValidNeighbors > 0) {
         WeightedDepth +=
             UnsmoothedBLDH(ICell) * ValidNeighbors * Mesh->AreaCellH(ICell);
         AreaSum += ValidNeighbors * Mesh->AreaCellH(ICell);
      }
      Real ExpectedDepth =
          AreaSum > 0.0_Real ? WeightedDepth / AreaSum : UnsmoothedBLDH(ICell);
      ExpectedDepth = Kokkos::fmax(MinDepth, ExpectedDepth);
      ExpectedDepth = Kokkos::fmin(MaxDepth, ExpectedDepth);

      I4 ExpectedIndex = NVertLayers - 1;
      for (I4 K = 0; K < NVertLayers - 1; ++K) {
         const Real ZAbove = LayerThickness * K;
         const Real ZBelow = LayerThickness * (K + 1);
         if (ExpectedDepth >= ZAbove && ExpectedDepth <= ZBelow) {
            ExpectedIndex = K;
            break;
         }
      }
      if (!isApprox(SmoothedBLDH(ICell), ExpectedDepth, BLDRTol, ATol) ||
          SmoothedIndexH(ICell) != ExpectedIndex) {
         ++NumErrors;
      }
      if (!isApprox(SmoothedBLDH(ICell), UnsmoothedBLDH(ICell), BLDRTol,
                    ATol)) {
         ++NumChanged;
      }
   }
   if (NumChanged == 0) {
      ++NumErrors;
   }
   KPPInstance->UseOSBLSmoothing = false;
   checkResult("boundary-layer horizontal smoothing", NumErrors);
}

void testEnabledFullCall() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;
   setCoefficientTestGeometry();
   VCoord->minMaxLayerEdge(Halo::getDefault());

   Array2DReal Density("KPPMixTest-FullCallDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-FullCallNormalVelocity",
                              Mesh->NEdgesSize, NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-FullCallTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-FullCallUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-FullCallB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-FullCallBVF", Mesh->NCellsSize, NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-FullCallIce", Mesh->NCellsSize);
   Array1DReal Wind;

   parallelFor(
       "KPPMixTest-SetFullCallDensity", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          Density(ICell, K) = RhoSw + 0.01_Real * K;
       });
   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, -1.0e-7_Real);
   deepCopy(BVF, 1.0e-4_Real);
   deepCopy(IceFraction, 0.0_Real);

   KPPInstance->Enabled               = true;
   KPPInstance->UseLangmuirTurbulence = false;
   KPPInstance->UseOSBLSmoothing      = false;
   KPPInstance->UseEnhancedDiffusion  = false;
   KPPInstance->MatchTechnique        = KPPMatchType::SimpleShapes;
   KPPInstance->computeOSBLDepth(Density, NormalVelocity, TangentialVelocity,
                                 UStar, B0, BVF, IceFraction, Wind);
   KPPInstance->computeMixingCoefficients(UStar, B0);

   Array1DReal ExpectedBLD("KPPMixTest-ExpectedBLD", Mesh->NCellsSize);
   Array1DI4 ExpectedBLDIndex("KPPMixTest-ExpectedBLDIndex", Mesh->NCellsSize);
   Array2DReal ExpectedBulkRi("KPPMixTest-ExpectedBulkRi", Mesh->NCellsSize,
                              NVertLayers + 1);
   Array2DReal ExpectedVertDiff("KPPMixTest-ExpectedVertDiff", Mesh->NCellsSize,
                                NVertLayers + 1);
   Array2DReal ExpectedVertVisc("KPPMixTest-ExpectedVertVisc", Mesh->NCellsSize,
                                NVertLayers + 1);
   Array2DReal ExpectedNonLocal("KPPMixTest-ExpectedNonLocal", Mesh->NCellsSize,
                                NVertLayers + 1);
   Array2DReal ExpectedTurbVel("KPPMixTest-ExpectedTurbVel", Mesh->NCellsSize,
                               NVertLayers + 1);
   deepCopy(ExpectedBLD, KPPInstance->OSBLDepth);
   deepCopy(ExpectedBLDIndex, KPPInstance->OSBLDepthIndex);
   deepCopy(ExpectedBulkRi, KPPInstance->BulkRichardsonNumber);
   deepCopy(ExpectedVertDiff, KPPInstance->VertDiff);
   deepCopy(ExpectedVertVisc, KPPInstance->VertVisc);
   deepCopy(ExpectedNonLocal, KPPInstance->VertNonLocalFlux);
   deepCopy(ExpectedTurbVel, KPPInstance->TurbulentVelocityScale);

   deepCopy(KPPInstance->OSBLDepth, -1.0_Real);
   deepCopy(KPPInstance->OSBLDepthIndex, -1);
   deepCopy(KPPInstance->BulkRichardsonNumber, -1.0_Real);
   deepCopy(KPPInstance->VertDiff, -1.0_Real);
   deepCopy(KPPInstance->VertVisc, -1.0_Real);
   deepCopy(KPPInstance->VertNonLocalFlux, -1.0_Real);
   deepCopy(KPPInstance->TurbulentVelocityScale, -1.0_Real);
   deepCopy(KPPInstance->PotentialDensity, -1.0_Real);

   KPPInstance->computeKPPMix(Density, NormalVelocity, TangentialVelocity,
                              UStar, B0, BVF, IceFraction, Wind);

   const auto ExpectedBLDH      = createHostMirrorCopy(ExpectedBLD);
   const auto ExpectedBLDIndexH = createHostMirrorCopy(ExpectedBLDIndex);
   const auto ExpectedBulkRiH   = createHostMirrorCopy(ExpectedBulkRi);
   const auto ExpectedVertDiffH = createHostMirrorCopy(ExpectedVertDiff);
   const auto ExpectedVertViscH = createHostMirrorCopy(ExpectedVertVisc);
   const auto ExpectedNonLocalH = createHostMirrorCopy(ExpectedNonLocal);
   const auto ExpectedTurbVelH  = createHostMirrorCopy(ExpectedTurbVel);
   const auto ActualBLDH        = createHostMirrorCopy(KPPInstance->OSBLDepth);
   const auto ActualBLDIndexH =
       createHostMirrorCopy(KPPInstance->OSBLDepthIndex);
   const auto ActualBulkRiH =
       createHostMirrorCopy(KPPInstance->BulkRichardsonNumber);
   const auto ActualVertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto ActualVertViscH = createHostMirrorCopy(KPPInstance->VertVisc);
   const auto ActualNonLocalH =
       createHostMirrorCopy(KPPInstance->VertNonLocalFlux);
   const auto ActualTurbVelH =
       createHostMirrorCopy(KPPInstance->TurbulentVelocityScale);
   const auto RetainedDensityH =
       createHostMirrorCopy(KPPInstance->PotentialDensity);
   const auto InputDensityH = createHostMirrorCopy(Density);

   int NumErrors = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (!isApprox(ActualBLDH(ICell), ExpectedBLDH(ICell), RTol, ATol) ||
          ActualBLDIndexH(ICell) != ExpectedBLDIndexH(ICell)) {
         ++NumErrors;
      }
      for (I4 K = 0; K <= NVertLayers; ++K) {
         if (!isApprox(ActualBulkRiH(ICell, K), ExpectedBulkRiH(ICell, K), RTol,
                       ATol) ||
             !isApprox(ActualVertDiffH(ICell, K), ExpectedVertDiffH(ICell, K),
                       RTol, ATol) ||
             !isApprox(ActualVertViscH(ICell, K), ExpectedVertViscH(ICell, K),
                       RTol, ATol) ||
             !isApprox(ActualNonLocalH(ICell, K), ExpectedNonLocalH(ICell, K),
                       RTol, ATol) ||
             !isApprox(ActualTurbVelH(ICell, K), ExpectedTurbVelH(ICell, K),
                       RTol, ATol)) {
            ++NumErrors;
         }
         if (K < NVertLayers &&
             !isApprox(RetainedDensityH(ICell, K), InputDensityH(ICell, K),
                       RTol, ATol)) {
            ++NumErrors;
         }
      }
   }
   checkResult("enabled full call", NumErrors);
}

void testDisabledFullCall() {
   const HorzMesh *Mesh = HorzMesh::getDefault();
   VertCoord *VCoord    = VertCoord::getDefault();
   KPPMix *KPPInstance  = KPPMix::getInstance();
   const I4 NVertLayers = VCoord->NVertLayers;

   Array2DReal Density("KPPMixTest-DisabledDensity", Mesh->NCellsSize,
                       NVertLayers);
   Array2DReal NormalVelocity("KPPMixTest-DisabledNormalVelocity",
                              Mesh->NEdgesSize, NVertLayers);
   Array2DReal TangentialVelocity("KPPMixTest-DisabledTangentialVelocity",
                                  Mesh->NEdgesSize, NVertLayers);
   Array1DReal UStar("KPPMixTest-DisabledUStar", Mesh->NCellsSize);
   Array1DReal B0("KPPMixTest-DisabledB0", Mesh->NCellsSize);
   Array2DReal BVF("KPPMixTest-DisabledBVF", Mesh->NCellsSize, NVertLayers + 1);
   Array1DReal IceFraction("KPPMixTest-DisabledIce", Mesh->NCellsSize);
   deepCopy(Density, RhoSw);
   deepCopy(NormalVelocity, 0.0_Real);
   deepCopy(TangentialVelocity, 0.0_Real);
   deepCopy(UStar, 0.02_Real);
   deepCopy(B0, -1.0e-7_Real);
   deepCopy(BVF, 0.0_Real);
   deepCopy(IceFraction, 0.0_Real);
   deepCopy(KPPInstance->VertDiff, -7.0_Real);
   deepCopy(KPPInstance->OSBLDepth, -9.0_Real);

   KPPInstance->Enabled = false;
   KPPInstance->computeKPPMix(Density, NormalVelocity, TangentialVelocity,
                              UStar, B0, BVF, IceFraction);
   const auto VertDiffH = createHostMirrorCopy(KPPInstance->VertDiff);
   const auto BLDH      = createHostMirrorCopy(KPPInstance->OSBLDepth);
   int NumErrors        = 0;
   for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      if (VertDiffH(ICell, 2) != -7.0_Real || BLDH(ICell) != -9.0_Real) {
         ++NumErrors;
      }
   }
   KPPInstance->Enabled = true;
   checkResult("disabled full call", NumErrors);
}

} // namespace

int main(int argc, char *argv[]) {
   const std::string TestGroup = argc > 1 ? argv[1] : "all";

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   initKPPMixTest(TestGroup);

   if (TestGroup == "profiles" || TestGroup == "all") {
      testStabilityFunctions();
      testShapeFunctions();
      testLangmuirFunctions();
      testOSBLUtilities();
      testTurbulentVelocityScale();
      testTurbScales();
      testMatchShape();
      testNonLocalCs();
      testClampOSBLDepth();
      testOSBLIndex();
   }
   if (TestGroup == "osbl" || TestGroup == "all") {
      testOSBLDepth();
      testBoundaryLayerNonuniformThickness();
      testBoundaryLayerHorizontalThicknessVariation();
      testSshOffsetInvariance();
      testBoundaryLayerEdgeFallbacks();
      testBoundaryLayerLangmuir();
      testBoundaryLayerSmoothing();
   }
   if (TestGroup == "vmix" || TestGroup == "all") {
      testWindOnlyCoefficients();
      testConvectionOnlyCoefficients();
      testNonLocalProfileModes();
      testMatchBothInteriorCoefficients();
      testEnhancedDiffusion();
      testStableAndZeroForcing();
      testCoefficientVerticalDomainEdges();
      testCoefficientInvalidWetBounds();
   }
   if (TestGroup == "integration" || TestGroup == "all") {
      testConfiguredValues();
      testEnabledFullCall();
      testDisabledFullCall();
   }
   if (TestGroup == "config-gradient" || TestGroup == "config-unsupported" ||
       TestGroup == "config-parabolic") {
      ABORT_ERROR("KPPMixTest: KPPMix::init should have rejected the injected "
                  "MatchTechnique for group '{}'",
                  TestGroup);
   }
   if (TestGroup != "profiles" && TestGroup != "osbl" && TestGroup != "vmix" &&
       TestGroup != "integration" && TestGroup != "all") {
      ABORT_ERROR("KPPMixTest: unknown test group '{}'", TestGroup);
   }

   LOG_INFO("------ KPP {} Tests Successful ------", TestGroup);
   finalizeKPPMixTest();
   Kokkos::finalize();
   MPI_Finalize();
   return 0;
}

//===----------------------------------------------------------------------===//
