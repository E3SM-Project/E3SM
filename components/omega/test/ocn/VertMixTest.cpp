//===-- Test driver for OMEGA Vertical Mixing Coefficients -------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Vertical Mixing Coefficients
///
/// This driver tests that VertMix can be called and returns expected values
/// of diffusivity, viscosity and Brunt-Vaisala frequency, and that the
/// implicit velocity vertical mixing solve, including Rayleigh damping,
/// produces the expected velocities
///
//===-----------------------------------------------------------------------===/

#include "VertMix.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Field.h"
#include "FillValues.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "KPPMix.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "mpi.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace OMEGA;

/// Test constants and expected values
constexpr int NVertLayers   = 60;
constexpr int NVertLayersP1 = NVertLayers + 1;

/// Values to test against
const Real VertDiffExpValueN =
    1.00501; // Expected value for diffusivity for positive BVF
const Real VertViscExpValueN =
    1.0051; // Expected value for viscosity for positive BVF
const Real VertDiffExpValueP =
    0.003882748571163051; // Expected value for diffusivity for negative BVF
const Real VertViscExpValueP =
    0.003972748571163051; // Expected value for viscosity for negative BVF
const Real VertDiffBackExp =
    1.0e-5; // Expected value for background diffusivity
const Real VertViscBackExp = 1.0e-4; // Expected value for background viscosity
const Real VertConvExp =
    1.0; // Expected value for convective diffusivity/viscosity
const Real VertShearExp =
    0.00387274859; // Expected value for shear diffusivity/viscosity
const Real VertShearBaseExp =
    0.005;                   // Expected value for shear diffusivity/viscosity
const Real RiExpValue = 0.2; // Expected value for gradient Richardson number

/// Test input values
const Real BVFP = 0.1;  // Positive Brunt-Vaisala frequency in s^-2
const Real BVFN = -0.1; // Negative Brunt-Vaisala frequency in s^-2
const Real NV   = 1.0;  // Normal velocity in m/s
const Real TV   = 1.0;  // Tangential velocity in m/s
const Real RTol = 1e-7; // Relative tolerance for isApprox checks

/// The initialization routine for VertMix testing. It calls various
/// init routines, including the creation of the default decomposition.
void initVertMixTest() {

   /// Initialize the Machine Environment class - this also creates
   /// the default MachEnv. Then retrieve the default environment and
   /// some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   /// Initialize logging
   initLogging(DefEnv);
   LOG_INFO("------ Vertical Mixing Unit Tests ------");

   /// Open and read config file
   Config("Omega");
   Config::readAll("omega.yml");

   /// Initialize the default time stepper, which provides the calendar, the
   /// model clock for stream IO and the time step used by the implicit
   /// vertical mixing solve
   TimeStepper::init1();
   Clock *ModelClock = TimeStepper::getDefault()->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   /// Initialize decomposition
   Decomp::init();

   /// Initialize Halo
   Halo::init();

   /// Initialize IO streams for mesh IO
   Field::init(ModelClock);
   IOStream::init(ModelClock);

   /// Initialize mesh
   HorzMesh::init(ModelClock);

   /// Initialize vertical coordinate
   VertCoord::init(false);

   /// Initialize the tracers, auxiliary state and equation of state needed
   /// by the implicit vertical mixing solve
   Tracers::init();
   VertAdv::init();
   AuxiliaryState::init();
   Eos::init();

   /// Initialize VertMix
   VertMix::init();

   /// Retrieve VertMix
   VertMix *DefVertMix = VertMix::getInstance();
   if (!DefVertMix)
      ABORT_ERROR("VertMixTest: VertMix retrieval FAIL");
}

void testGradRichNum() {
   /// Get mesh and coordinate info
   const auto Mesh       = HorzMesh::getDefault();
   const auto VCoord     = VertCoord::getDefault();
   auto *MeshHalo        = Halo::getDefault();
   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayersP1;
   I4 NCellsSize         = Mesh->NCellsSize;
   I4 NEdgesAll          = Mesh->NEdgesAll;
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(EdgesOnCell, Mesh->EdgesOnCell);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);
   OMEGA_SCOPE(CellsOnCell, Mesh->CellsOnCell);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesAll, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesAll, NVertLayers);
   auto BruntVaisalaFreqSqCell =
       Array2DReal("BruntVaisalaFreqSqCell", NCellsSize, NVertLayersP1);
   /// Use deep copy to initialize results
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(TestVertMix->GradRichNum, 0.0);

   // Also test guard is working for skipping layer edge contribution
   // if it would access invalid edge velocity levels
   parallelFor(
       "setMinMax", {Mesh->NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 0;
          if (ICell % 2 == 0) {
             MaxLayerCell(ICell) = 10;
          } else {
             MaxLayerCell(ICell) = NVertLayers - 2;
          }
       });

   // Need to exchange halo after changing MaxLayerCell
   MeshHalo->exchangeFullArrayHalo(MaxLayerCell, OnCell);

   // Refresh edge layer ranges after overriding MaxLayerCell.
   VCoord->minMaxLayerEdge(MeshHalo);

   parallelForOuter(
       "populateArrays", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
          AreaCell(ICell) = 3.6e10_Real;

          // Set inactive layers to fill value to cause failure if used
          parallelForInner(
              Team, NVertLayers, INNER_LAMBDA(I4 K) {
                 if (K > MaxLayerCell(ICell)) {
                    GeomZMid(ICell, K) = FillValueReal;
                 } else {
                    GeomZMid(ICell, K) = -K;
                 }
              });

          // Set inactive layers to fill value to cause failure if used
          parallelForInner(
              Team, NVertLayersP1, INNER_LAMBDA(I4 K) {
                 if (K > MaxLayerCell(ICell) + 1) {
                    BruntVaisalaFreqSqCell(ICell, K) = FillValueReal;
                 } else {
                    BruntVaisalaFreqSqCell(ICell, K) = BVFP;
                 }
              });
       });

   // Rebind GradRichardsonNum so the functor captures refreshed edge ranges.
   TestVertMix->ComputeGradRichardsonNum = GradRichardsonNum(Mesh, VCoord);

   // Recapture after minMaxLayerEdge since it reallocates edge layer views.
   OMEGA_SCOPE(MaxLayerEdgeBot, VCoord->MaxLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelFor(
       "populateArrays", {NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          if (K > MaxLayerEdgeBot(IEdge)) {
             NormalVelEdge(IEdge, K) =
                 FillValueReal; // Fill value to cause failure if used
          } else {
             NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
          }

          if (K > MaxLayerEdgeTop(IEdge)) {
             TangVelEdge(IEdge, K) =
                 FillValueReal; // Fill value to cause failure if used
          } else {
             TangVelEdge(IEdge, K) = TangVelEdge(IEdge, K) + 0.5 * K;
          }

          DcEdge(IEdge) = 2.0e5_Real;
          DvEdge(IEdge) = 1.45e5_Real;
       });

   /// Compute gradient Richardson number
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);

   /// Check all array values against expected value
   int NumMismatches = 0;
   OMEGA_SCOPE(GradRichNum, TestVertMix->GradRichNum);
   parallelReduceOuter(
       "CheckGradRichNum", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;

                 // Match the production K1/K2 logic to determine whether
                 // this layer has any valid edge contributions.
                 I4 K1 = K - 1;
                 I4 K2 = K;

                 if (K == MaxLayerCell(ICell) + 1) {
                    K1 = K - 2;
                    K2 = K - 1;
                 }

                 bool HasValidEdge = false;
                 for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
                    const I4 JEdge = EdgesOnCell(ICell, J);
                    if (K1 <= MaxLayerEdgeTop(JEdge) &&
                        K2 <= MaxLayerEdgeTop(JEdge)) {
                       HasValidEdge = true;
                       break;
                    }
                 }

                 if (HasValidEdge) {
                    if (!isApprox(GradRichNum(ICell, K), RiExpValue, RTol))
                       InnerCount++;
                 } else {
                    // With all edges skipped, GradRichNum remains at sentinel
                    // scale (~1e14) from initialization in the functor.
                    if (!(GradRichNum(ICell, K) > 1.0e10_Real))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   // If test fails, print bad values and abort
   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMix: GradRichNum FAIL with {} bad values",
                  NumMismatches);
   } else {
      LOG_INFO("TestVertMix: GradRichNum PASS");
   }

   // Reset MaxLayerCell to the original values
   parallelFor(
       "resetMinMax", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 0;
          MaxLayerCell(ICell) = NVertLayers - 1;
       });
   // Refresh edge layer ranges after overriding MaxLayerCell.
   VCoord->minMaxLayerEdge(MeshHalo);

   return;
}

void testOneTwoOneFilter() {
   /// Get mesh and coordinate info
   const auto Mesh       = HorzMesh::getDefault();
   const auto VCoord     = VertCoord::getDefault();
   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayersP1;
   I4 NCellsSize         = Mesh->NCellsSize;
   I4 NChunks            = VCoord->NVertLayers / VecLength;
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto GradRichNumSmoothed =
       Array2DReal("GradRichNumSmoothed", NCellsSize, NVertLayersP1);
   auto GradRichNum = Array2DReal("GradRichNum", NCellsSize, NVertLayersP1);
   /// Use deep copy to initialize results
   deepCopy(GradRichNumSmoothed, 1.0);
   deepCopy(GradRichNum, 1.0);

   // Populate GradRichNum with alternating +1.0 and -1.0 values in vertical
   // GradRichNumSmoothed should smooth these to 0.0
   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          if (K % 2 == 0) {
             GradRichNum(ICell, K) = 1.0;
          } else {
             GradRichNum(ICell, K) = -1.0;
          }
       });

   // Apply the 1-2-1 filter to each cell
   OMEGA_SCOPE(ComputeOneTwoOneFilter, TestVertMix->ComputeOneTwoOneFilter);

   const int RiSmoothLoopsOrig = ComputeOneTwoOneFilter.RiSmoothLoops;
   // Perform only one smoothing iteration for this test
   ComputeOneTwoOneFilter.RiSmoothLoops = 1;

   parallelForOuter(
       "ApplyOneTwoOneFilter",
       LaunchConfig({Mesh->NCellsAll},
                    TeamScratch<Real>(2 * VCoord->NVertLayersP1)),
       KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
          ComputeOneTwoOneFilter(Team, GradRichNumSmoothed, ICell, GradRichNum);
       });

   // Restore original number of smoothing iterations
   ComputeOneTwoOneFilter.RiSmoothLoops = RiSmoothLoopsOrig;

   /// Check all array values against expected value
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckGradRichNum", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K > MinLayerCell(ICell) && K <= MaxLayerCell(ICell)) {
                    // Interior layers should be smoothed to 0.0
                    if (!isApprox(GradRichNumSmoothed(ICell, K), 0.0_Real,
                                  RTol))
                       InnerCount++;
                 } else {
                    // Boundary layers (K==0 or K==NVertLayers) should be
                    // the same as input
                    if (!isApprox(GradRichNumSmoothed(ICell, K),
                                  GradRichNum(ICell, K), RTol))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   // If test fails, print bad values and abort
   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMix: GradRichNumSmoothed FAIL with {} bad values",
                  NumMismatches);
   } else {
      LOG_INFO("TestVertMix: GradRichNumSmoothed PASS");
   }

   return;
}

void testBackVertMix() {
   // Get mesh and coordinate info
   const auto Mesh       = HorzMesh::getDefault();
   const auto VCoord     = VertCoord::getDefault();
   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayersP1;
   I4 NCellsSize         = Mesh->NCellsSize;
   I4 NEdgesSize         = Mesh->NEdgesSize;
   I4 NEdgesAll          = Mesh->NEdgesAll;
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesSize, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesSize, NVertLayers);
   auto BruntVaisalaFreqSqCell =
       Array2DReal("BruntVaisalaFreqSqCell", NCellsSize, NVertLayersP1);

   /// Use deep copy initialize with reference or zero values
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(BruntVaisalaFreqSqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) { GeomZMid(ICell, K) = -K; });

   parallelFor(
       "populateArrays", {NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
          TangVelEdge(IEdge, K)   = TangVelEdge(IEdge, K) + 0.5 * K;
       });

   /// Compute only background vertical viscosity and diffusivity
   TestVertMix->BackDiff                    = 1.0e-5;
   TestVertMix->BackVisc                    = 1.0e-4;
   TestVertMix->ComputeVertMixConv.Enabled  = false;
   TestVertMix->ComputeVertMixShear.Enabled = false;
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   Array2DReal BackVertVisc = TestVertMix->VertVisc;
   Array2DReal BackVertDiff = TestVertMix->VertDiff;

   /// Check Visc against expected value
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-BackgroundVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    // Surface and bottom layers should be zero
                    if (BackVertVisc(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else {
                    if (!isApprox(BackVertVisc(ICell, K), VertViscBackExp,
                                  RTol))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMixBack: VertVisc FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertViscBackExp, BackVertVisc(1, 1), NumMismatches);
   } else {
      LOG_INFO("TestVertMixBack: VertVisc PASS");
   }

   /// Check Diff against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-BackgroundDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    // Surface and bottom layers should be zero
                    if (BackVertDiff(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else {
                    if (!isApprox(BackVertDiff(ICell, K), VertDiffBackExp,
                                  RTol))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      auto BackVertDiffH = createHostMirrorCopy(BackVertDiff);
      ABORT_ERROR("TestVertMixBack: VertDiff FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertDiffBackExp, BackVertDiffH(1, 1), NumMismatches);
   } else {
      LOG_INFO("TestVertMixBack: VertDiff PASS");
   }

   return;
}

void testConvVertMix() {
   // Get mesh and coordinate info
   const auto Mesh       = HorzMesh::getDefault();
   const auto VCoord     = VertCoord::getDefault();
   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayersP1;
   I4 NCellsSize         = Mesh->NCellsSize;
   I4 NChunks            = VCoord->NVertLayers / VecLength;

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto BruntVaisalaFreqSqIn =
       Array2DReal("BruntVaisalaFreqSqIn", NCellsSize, NVertLayersP1);
   auto VertDiffOut =
       Array2DReal("VertDiffOut", Mesh->NCellsAll, NVertLayersP1);
   auto VertViscOut =
       Array2DReal("VertViscOut", Mesh->NCellsAll, NVertLayersP1);

   /// Use deep copy to initialize with the ref value
   deepCopy(BruntVaisalaFreqSqIn, 0.0);
   deepCopy(VertDiffOut, 0.0);
   deepCopy(VertViscOut, 0.0);

   // Populate arrays: positive BVF in lower half (conv off),
   // negative in upper half (conv on)
   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          if (K < 30) {
             BruntVaisalaFreqSqIn(ICell, K) = -0.2;
          } else {
             BruntVaisalaFreqSqIn(ICell, K) = 0.2;
          }
       });

   /// Compute only convective vertical viscosity and diffusivity
   OMEGA_SCOPE(ComputeVertMixConv, TestVertMix->ComputeVertMixConv);
   Array1DI4 KPPBoundaryLayerIndex; // Empty, so KPP is disabled for this test
   parallelFor(
       "ApplyVertMixConv", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
          ComputeVertMixConv(VertDiffOut, VertViscOut, ICell, KChunk,
                             BruntVaisalaFreqSqIn, KPPBoundaryLayerIndex);
       });

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check Visc against expected value
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ConvectiveVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    // Surface and bottom layers should be zero
                    if (VertViscOut(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else if (K < 30) {
                    if (!isApprox(VertViscOut(ICell, K), VertConvExp, RTol))
                       InnerCount++;
                 } else {
                    if (!isApprox(VertViscOut(ICell, K), 0.0_Real, RTol))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMixConv: VertVisc FAIL with {} bad values",
                  NumMismatches);
   } else {
      LOG_INFO("TestVertMixConv: VertVisc PASS");
   }

   /// Check Diff against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ConvectiveDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    // Surface and bottom layers should be zero
                    if (VertDiffOut(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else if (K < 30) {
                    if (!isApprox(VertDiffOut(ICell, K), VertConvExp, RTol))
                       InnerCount++;
                 } else {
                    if (!isApprox(VertDiffOut(ICell, K), 0.0_Real, RTol))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMixConv: VertDiff FAIL with {} bad values",
                  NumMismatches);
   } else {
      LOG_INFO("TestVertMixConv: VertDiff PASS");
   }

   return;
}

void testShearVertMix() {
   /// Get mesh and coordinate info
   const auto Mesh       = HorzMesh::getDefault();
   const auto VCoord     = VertCoord::getDefault();
   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayersP1;
   I4 NCellsSize         = Mesh->NCellsSize;
   I4 NChunks            = VCoord->NVertLayers / VecLength;

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto GradRichNumSmoothedIn =
       Array2DReal("GradRichNumSmoothedIn", NCellsSize, NVertLayersP1);
   auto VertDiffOut =
       Array2DReal("VertDiffOut", Mesh->NCellsAll, NVertLayersP1);
   auto VertViscOut =
       Array2DReal("VertViscOut", Mesh->NCellsAll, NVertLayersP1);

   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(GradRichNumSmoothedIn, 0.0);
   deepCopy(VertDiffOut, 0.0);
   deepCopy(VertViscOut, 0.0);

   // Populate arrays: negative Ri in upper third (base shear value),
   // positive in middle third (altered shear value), large positive
   // in lower third (no shear)
   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          if (K < 20) {
             GradRichNumSmoothedIn(ICell, K) = -0.2;
          } else if (K >= 20 && K < 40) {
             GradRichNumSmoothedIn(ICell, K) = 0.2;
          } else {
             GradRichNumSmoothedIn(ICell, K) = 10.0;
          }
       });

   /// Compute only shear vertical viscosity and diffusivity
   OMEGA_SCOPE(ComputeVertMixShear, TestVertMix->ComputeVertMixShear);
   ComputeVertMixShear.ShearExponent = 3.0;
   parallelFor(
       "ApplyVertMixShear", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
          ComputeVertMixShear(VertDiffOut, VertViscOut, ICell, KChunk,
                              GradRichNumSmoothedIn);
       });

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check Visc against expected value
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ShearVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    // Surface and bottom layers should be zero
                    if (VertViscOut(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else if (K < 20) {
                    if (!isApprox(VertViscOut(ICell, K), VertShearBaseExp,
                                  RTol))
                       InnerCount++;
                 } else if (K >= 20 && K < 40) {
                    if (!isApprox(VertViscOut(ICell, K), VertShearExp, RTol))
                       InnerCount++;
                 } else {
                    if (!isApprox(VertViscOut(ICell, K), 0.0_Real, RTol))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMixShear: VertVisc FAIL with {} bad values",
                  NumMismatches);
   } else {
      LOG_INFO("TestVertMixShear: VertVisc PASS");
   }

   /// Check Diff against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ShearVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    // Surface and bottom layers should be zero
                    if (VertDiffOut(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else if (K < 20) {
                    if (!isApprox(VertDiffOut(ICell, K), VertShearBaseExp,
                                  RTol))
                       InnerCount++;
                 } else if (K >= 20 && K < 40) {
                    if (!isApprox(VertDiffOut(ICell, K), VertShearExp, RTol))
                       InnerCount++;
                 } else {
                    if (!isApprox(VertDiffOut(ICell, K), 0.0_Real, RTol))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMixShear: VertDiff FAIL with {} bad values",
                  NumMismatches);
   } else {
      LOG_INFO("TestVertMixShear: VertDiff PASS");
   }

   return;
}

/// Test vertical mixing coefficients calculation for all cells/layers
void testTotalVertMix() {
   /// Get mesh and coordinate info
   const auto Mesh       = HorzMesh::getDefault();
   const auto VCoord     = VertCoord::getDefault();
   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayersP1;
   I4 NCellsSize         = Mesh->NCellsSize;
   I4 NEdgesAll          = Mesh->NEdgesAll;
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);
   OMEGA_SCOPE(CellsOnCell, Mesh->CellsOnCell);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesAll, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesAll, NVertLayers);
   auto BruntVaisalaFreqSqCell =
       Array2DReal("BruntVaisalaFreqSqCell", NCellsSize, NVertLayersP1);

   /// Use deep copy to initialize with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);

   // Test with positive BVF first
   deepCopy(BruntVaisalaFreqSqCell, BVFP);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);
   deepCopy(TestVertMix->GradRichNumSmoothed, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          GeomZMid(ICell, K) = -K;
          AreaCell(ICell)    = 3.6e10_Real;
       });

   parallelFor(
       "populateArrays", {NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
          TangVelEdge(IEdge, K)   = TangVelEdge(IEdge, K) + 0.5 * K;
          DcEdge(IEdge)           = 2.0e5_Real;
          DvEdge(IEdge)           = 1.45e5_Real;
       });

   /// Compute vertical viscosity and diffusivity
   TestVertMix->BackDiff                    = 1.0e-5;
   TestVertMix->BackVisc                    = 1.0e-4;
   TestVertMix->ComputeVertMixConv.Enabled  = true;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   OMEGA_SCOPE(VertDiffP, TestVertMix->VertDiff);
   OMEGA_SCOPE(VertViscP, TestVertMix->VertVisc);

   /// Check all VertDiff array values against expected value
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-TotalPosDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    // Surface and bottom layers should be zero
                    if (VertDiffP(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == KMin + 1) {
                    if (!isApprox(VertDiffP(ICell, K), VertDiffExpValueP, RTol))
                       InnerCount++;
                    // otherwise check for invalid values
                 } else {
                    if (VertDiffP(ICell, K) == 0.0 or
                        Kokkos::isnan(VertDiffP(ICell, K)) or
                        Kokkos::isinf(VertDiffP(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      auto VertDiffPH = createHostMirrorCopy(VertDiffP);
      ABORT_ERROR("TestVertMixTotal: VertDiffPositive FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertDiffExpValueP, VertDiffPH(1, 1), NumMismatches);
   } else {
      LOG_INFO("TestVertMixTotal: VertDiffPos PASS");
   }

   /// Check all VertVisc array values against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-TotalPosVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    if (VertViscP(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == KMin + 1) {
                    if (!isApprox(VertViscP(ICell, K), VertViscExpValueP, RTol))
                       InnerCount++;
                    // otherwise check for invalid values
                 } else {
                    if (VertViscP(ICell, K) == 0.0 or
                        Kokkos::isnan(VertViscP(ICell, K)) or
                        Kokkos::isinf(VertViscP(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      auto VertViscPH = createHostMirrorCopy(VertViscP);
      ABORT_ERROR("TestVertMixTotal: VertViscPositive FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertViscExpValueP, VertViscPH(1, 1), NumMismatches);
   } else {
      LOG_INFO("TestVertMixTotal: VertViscPos PASS");
   }

   // Now test with negative BVF
   deepCopy(BruntVaisalaFreqSqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   /// Compute vertical viscosity and diffusivity
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);
   OMEGA_SCOPE(VertDiffN, TestVertMix->VertDiff);
   OMEGA_SCOPE(VertViscN, TestVertMix->VertVisc);

   /// Check all VertDiff array values against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-TotalNegDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    if (VertDiffN(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == KMin + 1) {
                    if (!isApprox(VertDiffN(ICell, K), VertDiffExpValueN, RTol))
                       InnerCount++;
                    // otherwise check for invalid values
                 } else {
                    if (VertDiffN(ICell, K) == 0.0 or
                        Kokkos::isnan(VertDiffN(ICell, K)) or
                        Kokkos::isinf(VertDiffN(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      auto VertDiffNH = createHostMirrorCopy(VertDiffN);
      ABORT_ERROR("TestVertMix: VertDiffNegative FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertDiffExpValueN, VertDiffNH(1, 1), NumMismatches);
   } else {
      LOG_INFO("TestVertMixTotal: VertDiffNeg PASS");
   }

   /// Check all VertVisc array values against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-TotalNegVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == KMin || K == KMax) {
                    if (VertViscN(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == KMin + 1) {
                    if (!isApprox(VertViscN(ICell, K), VertViscExpValueN, RTol))
                       InnerCount++;
                    // otherwise check for invalid values
                 } else {
                    if (VertViscN(ICell, K) == 0.0 or
                        Kokkos::isnan(VertViscN(ICell, K)) or
                        Kokkos::isinf(VertViscN(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      auto VertViscNH = createHostMirrorCopy(VertViscN);
      ABORT_ERROR("TestVertMix: VertViscNegative FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertViscExpValueN, VertViscNH(1, 1), NumMismatches);
   } else {
      LOG_INFO("TestVertMixTotal: VertViscNeg PASS");
   }

   return;
}

// This tests that VertMix's KPP-MatchBoth interior snapshot captures
// background+shear only (excluding both KPP and convective contributions),
// and that it is left untouched when MatchBoth is not the active technique.
void testInteriorSnapshotForMatchBoth() {
   const auto Mesh          = HorzMesh::getDefault();
   const auto VCoord        = VertCoord::getDefault();
   const auto &MinLayerCell = VCoord->MinLayerCell;
   I4 NEdgesAll             = Mesh->NEdgesAll;

   VertMix *TestVertMix = VertMix::getInstance();

   // KPP contributes zero diffusivity/viscosity here; only its
   // MatchTechnique setting is exercised.
   KPPMix::init();
   KPPMix *KPPInstance         = KPPMix::getInstance();
   KPPInstance->Enabled        = true;
   KPPInstance->MatchTechnique = KPPMatchType::MatchBoth;
   deepCopy(KPPInstance->VertDiff, 0.0_Real);
   deepCopy(KPPInstance->VertVisc, 0.0_Real);
   deepCopy(KPPInstance->OSBLDepthIndex, 0);

   auto NormalVelEdge =
       Array2DReal("NormalVelEdgeInterior", NEdgesAll, NVertLayers);
   auto TangVelEdge =
       Array2DReal("TangVelEdgeInterior", NEdgesAll, NVertLayers);
   auto BruntVaisalaFreqSqCell = Array2DReal("BruntVaisalaFreqSqCellInterior",
                                             Mesh->NCellsSize, NVertLayersP1);
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   // Strongly negative BVF everywhere so convective mixing definitely
   // triggers, distinguishing the (conv-excluded) snapshot from the final
   // merged VertDiff/VertVisc.
   deepCopy(BruntVaisalaFreqSqCell, -1.0_Real);

   constexpr Real TestConvDiff              = 2.0_Real;
   TestVertMix->ComputeVertMixConv.Enabled  = true;
   TestVertMix->ComputeVertMixConv.ConvDiff = TestConvDiff;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   deepCopy(TestVertMix->InteriorVertDiff, -1.0_Real);
   deepCopy(TestVertMix->InteriorVertVisc, -1.0_Real);

   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);

   OMEGA_SCOPE(VertDiff, TestVertMix->VertDiff);
   OMEGA_SCOPE(VertVisc, TestVertMix->VertVisc);
   OMEGA_SCOPE(InteriorDiff, TestVertMix->InteriorVertDiff);
   OMEGA_SCOPE(InteriorVisc, TestVertMix->InteriorVertVisc);

   constexpr I4 TestLevel = 10; // deep enough to avoid OSBL/boundary edges

   int NumMismatches = 0;
   parallelReduce(
       "CheckInteriorSnapshot", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, int &Count) {
          const int K       = MinLayerCell(ICell) + TestLevel;
          const bool DiffOk = isApprox(
              VertDiff(ICell, K), InteriorDiff(ICell, K) + TestConvDiff, RTol);
          const bool ViscOk = isApprox(
              VertVisc(ICell, K), InteriorVisc(ICell, K) + TestConvDiff, RTol);
          const bool SnapshotFilled = InteriorDiff(ICell, K) != -1.0_Real;
          if (!DiffOk || !ViscOk || !SnapshotFilled)
             Count++;
       },
       NumMismatches);

   if (NumMismatches != 0) {
      ABORT_ERROR("TestVertMix: InteriorSnapshotForMatchBoth FAIL, {} "
                  "mismatches",
                  NumMismatches);
   } else {
      LOG_INFO("TestVertMix: InteriorSnapshotForMatchBoth (MatchBoth) PASS");
   }

   // Confirm the snapshot is skipped (left untouched) when MatchBoth is not
   // the configured technique.
   KPPInstance->MatchTechnique = KPPMatchType::SimpleShapes;
   deepCopy(TestVertMix->InteriorVertDiff, -1.0_Real);
   deepCopy(TestVertMix->InteriorVertVisc, -1.0_Real);
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);

   int NumUnexpectedUpdates = 0;
   parallelReduce(
       "CheckInteriorSnapshotSkipped", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, int &Count) {
          const int K = MinLayerCell(ICell) + TestLevel;
          if (InteriorDiff(ICell, K) != -1.0_Real ||
              InteriorVisc(ICell, K) != -1.0_Real)
             Count++;
       },
       NumUnexpectedUpdates);

   if (NumUnexpectedUpdates != 0) {
      ABORT_ERROR("TestVertMix: InteriorSnapshotSkippedForSimpleShapes FAIL, "
                  "{} mismatches",
                  NumUnexpectedUpdates);
   } else {
      LOG_INFO("TestVertMix: InteriorSnapshotSkippedForSimpleShapes PASS");
   }
}

//------------------------------------------------------------------------------
// Tests of Rayleigh damping in the implicit velocity vertical mixing solve

/// Name of the ocean state used by the implicit vertical mixing tests
const std::string ImplicitStateName = "VertMixImplicitTest";

/// Test input values for the implicit vertical mixing tests
constexpr Real ImplicitNV     = 0.3;    // Normal velocity scale in m/s
constexpr Real ImplicitDrag   = 1.0e-3; // Bottom drag coefficient
constexpr Real ImplicitViscSc = 1.0e-3; // Vertical viscosity scale in m^2/s

/// Relative tolerance for comparing implicit solve results with expected
/// values, which differ only by round-off
const Real SolveRTol = 1000 * std::numeric_limits<Real>::epsilon();

/// Initial normal velocity, which alternates in sign between edges so that
/// both signs are tested
KOKKOS_INLINE_FUNCTION Real initNormalVel(I4 IEdge, I4 K) {
   const Real Sign = (IEdge % 2 == 0) ? 1.0_Real : -1.0_Real;
   return Sign * ImplicitNV * (1.0_Real + 0.01_Real * K);
}

/// Return the time step (s) used by VertMix::applyVelVertMixImplicit
R8 getImplicitTimeStep() {
   R8 DT;
   TimeStepper::getDefault()->getTimeStep().get(DT, TimeUnits::Seconds);
   return DT;
}

/// Set the state, auxiliary state, specific volume and viscosity used by the
/// implicit velocity vertical mixing solve. Pseudo-thickness is nonuniform in
/// both the horizontal and the vertical. If UseVisc is false, the vertical
/// viscosity is zero so that the layers decouple.
OceanState *initImplicitTestState(bool UseVisc) {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   auto *MeshHalo    = Halo::getDefault();
   auto *AuxState    = AuxiliaryState::getDefault();
   VertMix *VMix     = VertMix::getInstance();
   Eos *EqState      = Eos::getInstance();
   const I4 NLayers  = VCoord->NVertLayers;

   OceanState *State = OceanState::get(ImplicitStateName);
   if (!State) {
      State = OceanState::create(ImplicitStateName, Mesh, MeshHalo, NLayers, 2);
      if (!State)
         ABORT_ERROR("TestVertMixImplicit: error creating test state");
   }

   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal NormalVel   = State->getNormalVelocity(0);
   Array2DReal SpecVol     = EqState->SpecVol;
   Array2DReal VertVisc    = VMix->VertVisc;
   Array2DReal KECell      = AuxState->KineticAux.KineticEnergyCell;

   parallelFor(
       "initImplicitCells", {Mesh->NCellsAll, NLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          PseudoThick(ICell, K) =
              1.0_Real + 0.5_Real * K + 0.1_Real * (ICell % 7);
          SpecVol(ICell, K) = (1.0_Real + 1.0e-3_Real * (K % 3)) / RhoSw;
          KECell(ICell, K)  = 0.01_Real * (1 + (ICell + K) % 3);
       });

   parallelFor(
       "initImplicitEdges", {Mesh->NEdgesAll, NLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVel(IEdge, K) = initNormalVel(IEdge, K);
       });

   if (UseVisc) {
      parallelFor(
          "initImplicitVisc", {Mesh->NCellsAll, NLayers + 1},
          KOKKOS_LAMBDA(I4 ICell, I4 K) {
             VertVisc(ICell, K) = ImplicitViscSc * (1 + (ICell + K) % 5);
          });
   } else {
      deepCopy(VertVisc, 0.0_Real);
   }

   return State;
}

/// Return a host copy of an array. Unlike createHostMirrorCopy, the result
/// never aliases the array, even when host and device memory are the same.
HostArray2DReal copyToHost(const Array2DReal &Arr) {
   const auto ArrMirror = createHostMirrorCopy(Arr);
   HostArray2DReal ArrH("ArrH", Arr.extent(0), Arr.extent(1));
   deepCopy(ArrH, ArrMirror);
   return ArrH;
}

/// Run the implicit velocity vertical mixing solve with the given options
/// and return a host copy of the resulting normal velocity. A host copy of
/// the initial normal velocity is returned in InitVel. The initial velocity
/// is computed on the device and may differ in the last bit from the same
/// expression evaluated on the host (e.g. from FMA contraction), so it must
/// be copied, not recomputed, for exact comparisons.
HostArray2DReal runImplicitVelVertMix(bool UseVisc, bool BottomDragEnabled,
                                      bool RayleighEnabled, Real RayleighCoeff,
                                      HostArray2DReal &InitVel) {
   VertMix *VMix = VertMix::getInstance();
   auto *State   = initImplicitTestState(UseVisc);
   InitVel       = copyToHost(State->getNormalVelocity(0));

   VMix->VelVertMixSetup.Enabled                   = true;
   VMix->VelVertMixSetup.ImplicitBottomDragEnabled = BottomDragEnabled;
   VMix->VelVertMixSetup.BottomDragCoeff           = ImplicitDrag;
   VMix->VelVertMixSetup.RayleighDampingEnabled    = RayleighEnabled;
   VMix->VelVertMixSetup.RayleighDampingCoeff      = RayleighCoeff;

   VMix->applyVelVertMixImplicit(State, AuxiliaryState::getDefault(), 0, 0);

   // Restore the default (disabled) state of the optional terms
   VMix->VelVertMixSetup.ImplicitBottomDragEnabled = false;
   VMix->VelVertMixSetup.BottomDragCoeff           = 0.0_Real;
   VMix->VelVertMixSetup.RayleighDampingEnabled    = false;
   VMix->VelVertMixSetup.RayleighDampingCoeff      = 0.0_Real;

   return copyToHost(State->getNormalVelocity(0));
}

/// Independently assemble the tridiagonal system for the implicit velocity
/// solve on the host and solve it with the Thomas algorithm, starting from
/// the initial normal velocity InitVel. Assumes that initImplicitTestState
/// has been called with the same UseVisc.
HostArray2DReal computeReferenceVel(const HostArray2DReal &InitVel,
                                    bool BottomDragEnabled, R8 RayleighCoeff) {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   const I4 NLayers  = VCoord->NVertLayers;
   const R8 DT       = getImplicitTimeStep();
   auto *AuxState    = AuxiliaryState::getDefault();
   const auto *State = OceanState::get(ImplicitStateName);
   const auto HH     = createHostMirrorCopy(State->getPseudoThickness(0));
   const auto SpecVH = createHostMirrorCopy(Eos::getInstance()->SpecVol);
   const auto ViscH  = createHostMirrorCopy(VertMix::getInstance()->VertVisc);
   const auto KEH =
       createHostMirrorCopy(AuxState->KineticAux.KineticEnergyCell);
   const auto CellsOnEdgeH = createHostMirrorCopy(Mesh->CellsOnEdge);
   const auto KMinH        = createHostMirrorCopy(VCoord->MinLayerEdgeBot);
   const auto KMaxH        = createHostMirrorCopy(VCoord->MaxLayerEdgeTop);

   // The state's normal velocity has already been updated by the implicit
   // solve, so the initial velocity InitVel is used instead
   HostArray2DReal RefVel("RefVel", Mesh->NEdgesSize, NLayers);
   for (I4 IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
      for (I4 K = 0; K < NLayers; ++K) {
         RefVel(IEdge, K) = InitVel(IEdge, K);
      }
   }

   std::vector<R8> Lower(NLayers), Diag(NLayers), Upper(NLayers), Rhs(NLayers);

   for (I4 IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
      const I4 C0   = CellsOnEdgeH(IEdge, 0);
      const I4 C1   = CellsOnEdgeH(IEdge, 1);
      const I4 KMin = KMinH(IEdge);
      const I4 KMax = KMaxH(IEdge);

      // Skip columns with no active layers
      if (KMax < KMin)
         continue;

      auto edgeAvg = [&](const auto &Arr, I4 K) {
         return 0.5 * (R8(Arr(C0, K)) + R8(Arr(C1, K)));
      };

      // Coupling between layers K and K+1 through the interface between them
      std::vector<R8> Coupling(NLayers, 0.0);
      for (I4 K = KMin; K < KMax; ++K) {
         const R8 HK   = edgeAvg(HH, K);
         const R8 HKp1 = edgeAvg(HH, K + 1);
         const R8 HBot = 0.5 * (HK + HKp1);
         const R8 AlpB =
             (edgeAvg(SpecVH, K) * HKp1 + edgeAvg(SpecVH, K + 1) * HK) /
             (HK + HKp1);
         const R8 Visc = edgeAvg(ViscH, K + 1);
         Coupling[K]   = DT * Visc / (R8(RhoSw) * AlpB) / HBot;
      }

      for (I4 K = KMin; K <= KMax; ++K) {
         const R8 HK = edgeAvg(HH, K);
         R8 DiagOnly = HK * (1.0 + DT * RayleighCoeff);
         if (BottomDragEnabled && K == KMax) {
            const R8 VelMag = std::sqrt(R8(KEH(C0, K)) + R8(KEH(C1, K)));
            DiagOnly += DT * R8(ImplicitDrag) * VelMag /
                        (R8(RhoSw) * edgeAvg(SpecVH, K));
         }
         const R8 GAbove = (K > KMin) ? Coupling[K - 1] : 0.0;
         const R8 GBelow = (K < KMax) ? Coupling[K] : 0.0;
         Lower[K]        = -GAbove;
         Upper[K]        = -GBelow;
         Diag[K]         = GAbove + GBelow + DiagOnly;
         Rhs[K]          = HK * R8(InitVel(IEdge, K));
      }

      // Thomas algorithm: forward elimination and back substitution
      for (I4 K = KMin + 1; K <= KMax; ++K) {
         const R8 W = Lower[K] / Diag[K - 1];
         Diag[K] -= W * Upper[K - 1];
         Rhs[K] -= W * Rhs[K - 1];
      }
      R8 Below            = Rhs[KMax] / Diag[KMax];
      RefVel(IEdge, KMax) = Below;
      for (I4 K = KMax - 1; K >= KMin; --K) {
         Below            = (Rhs[K] - Upper[K] * Below) / Diag[K];
         RefVel(IEdge, K) = Below;
      }
   }

   return RefVel;
}

/// Count mismatches between the implicit solve result and the expected
/// velocity on owned edges. Active layers must match to within SolveRTol,
/// and inactive layers must be unchanged from the initial velocity InitVel.
int countImplicitMismatches(const HostArray2DReal &Vel,
                            const HostArray2DReal &ExpectedVel,
                            const HostArray2DReal &InitVel) {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   const auto KMinH  = createHostMirrorCopy(VCoord->MinLayerEdgeBot);
   const auto KMaxH  = createHostMirrorCopy(VCoord->MaxLayerEdgeTop);

   int NumMismatches = 0;
   for (I4 IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
      for (I4 K = 0; K < VCoord->NVertLayers; ++K) {
         if (K < KMinH(IEdge) || K > KMaxH(IEdge)) {
            if (Vel(IEdge, K) != InitVel(IEdge, K))
               NumMismatches++;
         } else if (!std::isfinite(Vel(IEdge, K)) ||
                    !isApprox(Vel(IEdge, K), ExpectedVel(IEdge, K),
                              SolveRTol)) {
            NumMismatches++;
         }
      }
   }
   return NumMismatches;
}

/// With no vertical viscosity and no bottom drag, the layers decouple and
/// the implicit damping gives u^{n+1} = u^n / (1 + DT * Coeff) in every
/// active layer. The analytic factor is independent of pseudo-thickness, so
/// this also checks the pseudo-thickness weighting of the damping term. For
/// large DT * Coeff, this also checks that the decay is stable, monotone,
/// and preserves the sign of the velocity.
void testRayleighDampingDecay(const std::string &TestName, Real Coeff) {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   const R8 DT       = getImplicitTimeStep();
   const R8 Factor   = 1.0 / (1.0 + DT * Coeff);

   HostArray2DReal InitVel;
   HostArray2DReal Vel =
       runImplicitVelVertMix(false, false, true, Coeff, InitVel);

   HostArray2DReal ExpectedVel("ExpectedVel", Vel.extent(0), Vel.extent(1));
   int NumSignErrors = 0;
   for (I4 IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
      for (I4 K = 0; K < VCoord->NVertLayers; ++K) {
         const Real U0         = InitVel(IEdge, K);
         ExpectedVel(IEdge, K) = Factor * U0;
         if (Vel(IEdge, K) * U0 <= 0.0_Real ||
             std::abs(Vel(IEdge, K)) > std::abs(U0))
            NumSignErrors++;
      }
   }

   const int NumMismatches = countImplicitMismatches(Vel, ExpectedVel, InitVel);

   if (NumMismatches != 0 || NumSignErrors != 0) {
      ABORT_ERROR("TestVertMixRayleigh: {} FAIL with DT * Coeff = {}, "
                  "{} mismatches and {} sign or growth errors",
                  TestName, DT * Coeff, NumMismatches, NumSignErrors);
   } else {
      LOG_INFO("TestVertMixRayleigh: {} PASS", TestName);
   }
}

/// Disabled Rayleigh damping, and enabled damping with a zero coefficient,
/// must give bit-for-bit identical results, and both must match the system
/// without damping
void testRayleighDampingNeutral() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   HostArray2DReal InitVel;
   HostArray2DReal DisabledVel =
       runImplicitVelVertMix(true, true, false, 1.0e-4_Real, InitVel);
   HostArray2DReal RefVel = computeReferenceVel(InitVel, true, 0.0);
   HostArray2DReal ZeroCoeffInitVel;
   HostArray2DReal ZeroCoeffVel =
       runImplicitVelVertMix(true, true, true, 0.0_Real, ZeroCoeffInitVel);

   int NumDiffs = 0;
   for (I4 IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
      for (I4 K = 0; K < VCoord->NVertLayers; ++K) {
         if (DisabledVel(IEdge, K) != ZeroCoeffVel(IEdge, K))
            NumDiffs++;
      }
   }

   const int NumMismatches =
       countImplicitMismatches(DisabledVel, RefVel, InitVel);

   if (NumDiffs != 0 || NumMismatches != 0) {
      ABORT_ERROR("TestVertMixRayleigh: Neutral FAIL with {} differences "
                  "from zero coefficient and {} mismatches from reference",
                  NumDiffs, NumMismatches);
   } else {
      LOG_INFO("TestVertMixRayleigh: Neutral PASS");
   }
}

/// Rayleigh damping combined with nonuniform vertical viscosity and implicit
/// bottom drag must match an independent host solve of the same system
void testRayleighDampingCombined() {
   const Real Coeff = 1.0e-3_Real;

   HostArray2DReal InitVel;
   HostArray2DReal Vel =
       runImplicitVelVertMix(true, true, true, Coeff, InitVel);
   HostArray2DReal RefVel = computeReferenceVel(InitVel, true, Coeff);

   // Also make sure the damping changed the result
   HostArray2DReal UndampedInitVel;
   HostArray2DReal UndampedVel =
       runImplicitVelVertMix(true, true, false, Coeff, UndampedInitVel);

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   const auto KMinH  = createHostMirrorCopy(VCoord->MinLayerEdgeBot);
   const auto KMaxH  = createHostMirrorCopy(VCoord->MaxLayerEdgeTop);
   int NumUnchanged  = 0;
   for (I4 IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
      for (I4 K = KMinH(IEdge); K <= KMaxH(IEdge); ++K) {
         if (!(std::abs(Vel(IEdge, K)) < std::abs(UndampedVel(IEdge, K))))
            NumUnchanged++;
      }
   }

   const int NumMismatches = countImplicitMismatches(Vel, RefVel, InitVel);

   if (NumMismatches != 0 || NumUnchanged != 0) {
      ABORT_ERROR("TestVertMixRayleigh: Combined FAIL with {} mismatches "
                  "from reference and {} layers not damped",
                  NumMismatches, NumUnchanged);
   } else {
      LOG_INFO("TestVertMixRayleigh: Combined PASS");
   }
}

/// Run all of the Rayleigh damping tests
void testRayleighDamping() {
   const R8 DT = getImplicitTimeStep();

   // Typical coefficient for the first segment of a dynamic adjustment
   testRayleighDampingDecay("Decay", 1.0e-4_Real);

   // DT * Coeff = 1000, far outside the range an explicit scheme could
   // tolerate
   testRayleighDampingDecay("LargeCoeff", Real(1.0e3 / DT));

   testRayleighDampingNeutral();
   testRayleighDampingCombined();
}

/// Finalize and clean up all test infrastructure
void finalizeVertMixTest() {
   KPPMix::destroyInstance();
   VertMix::destroyInstance();
   Eos::destroyInstance();
   OceanState::clear();
   AuxiliaryState::clear();
   VertAdv::clear();
   Tracers::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   Halo::clear();
   VertCoord::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

// the main tests (all in one to have the same log):
// --> one tests the gradient richardson number calculation
// --> next tests the 1-2-1 filter (smoothing)
// --> next tests the vertical diffusivity and viscosity
// with only background on
// --> next tests the vertical diffusivity and viscosity
// for only convective
// --> next tests the vertical diffusivity and viscosity
// for only shear
// --> next tests the linear superposition of the
// background, convective, and shear contributions
// --> finally tests Rayleigh damping in the implicit
// velocity vertical mixing solve
void vertMixTest() {

   // initialize vertical mix and other infrastructure
   initVertMixTest();

   // test each vertical mix option
   testGradRichNum();
   testOneTwoOneFilter();
   testBackVertMix();
   testConvVertMix();
   testShearVertMix();
   testTotalVertMix();
   testInteriorSnapshotForMatchBoth();

   // test Rayleigh damping in the implicit velocity vertical mixing solve
   testRayleighDamping();

   // clean up
   finalizeVertMixTest();
}

// The test driver for VertMix testing
int main(int argc, char *argv[]) {

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   vertMixTest();

   LOG_INFO("------ Vertical Mixing Unit Tests Successful ------");
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   // if we made it here, it is successful
   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
