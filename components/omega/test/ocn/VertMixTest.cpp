//===-- Test driver for OMEGA Vertical Mixing Coefficients -------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Vertical Mixing Coefficients
///
/// This driver tests that VertMix can be called and returns expected values
/// of diffusivity, viscosity and Brunt-Vaisala frequency
///
//===-----------------------------------------------------------------------===/

#include "VertMix.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "VertCoord.h"
#include "mpi.h"

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

   // Initialize the IO system
   IO::init(DefComm);

   /// Initialize decomposition
   Decomp::init();

   /// Initialize Halo
   Halo::init();

   /// Create dummy model clock for stream IO
   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   Clock ModelClockTmp(StartTime, TimeStep);
   Clock *ModelClock = &ModelClockTmp;

   /// Initialize IO streams for mesh IO
   Field::init(ModelClock);
   IOStream::init(ModelClock);

   /// Initialize mesh
   HorzMesh::init(ModelClock);

   /// Initialize vertical coordinate
   VertCoord::init(false);

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
   deepCopy(BruntVaisalaFreqSqCell, BVFP);
   deepCopy(TestVertMix->GradRichNum, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          GeomZMid(ICell, K)  = -K;
          NEdgesOnCell(ICell) = 5;
          AreaCell(ICell)     = 3.6e10_Real;
       });

   // Also test guard is working for skipping layer edge contribution
   // if it would access invalid edge velocity levels
   parallelFor(
       "setMinMax", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 0;
          if (ICell % 2 == 0) {
             MaxLayerCell(ICell) = 10;
          } else {
             MaxLayerCell(ICell) = NVertLayers - 2;
          }
       });

   // Refresh edge layer ranges after overriding MaxLayerCell.
   VCoord->minMaxLayerEdge(MeshHalo);

   // Rebind GradRichardsonNum so the functor captures refreshed edge ranges.
   TestVertMix->ComputeGradRichardsonNum = GradRichardsonNum(Mesh, VCoord);

   // Recapture after minMaxLayerEdge since it reallocates edge layer views.
   OMEGA_SCOPE(MaxLayerEdgeBot, VCoord->MaxLayerEdgeBot);

   // filling CellsOnCell with simple mapping for this test
   parallelFor(
       "populateArrays", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          CellsOnCell(ICell, 0) = ICell;
          CellsOnCell(ICell, 1) = ICell;
          CellsOnCell(ICell, 2) = ICell;
          CellsOnCell(ICell, 3) = ICell;
          CellsOnCell(ICell, 4) = ICell;
       });

   parallelFor(
       "populateArrays", {NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          if (K > MaxLayerEdgeBot(IEdge)) {
             NormalVelEdge(IEdge, K) =
                 -9.99e30; // Fill value to cause failure if used
             TangVelEdge(IEdge, K) =
                 -9.99e30; // Fill value to cause failure if used
          } else {
             NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
             TangVelEdge(IEdge, K)   = TangVelEdge(IEdge, K) + 0.5 * K;
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
                    if (K1 <= MaxLayerEdgeBot(JEdge) &&
                        K2 <= MaxLayerEdgeBot(JEdge)) {
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

   parallelFor(
       "setMinMax", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          MinLayerCell(ICell) = 0;
          MaxLayerCell(ICell) = NVertLayers - 1;
       });

   // Apply the 1-2-1 filter to each cell
   OMEGA_SCOPE(ComputeOneTwoOneFilter, TestVertMix->ComputeOneTwoOneFilter);
   parallelFor(
       "ApplyOneTwoOneFilter", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
          ComputeOneTwoOneFilter(GradRichNumSmoothed, ICell, KChunk,
                                 GradRichNum);
       });

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
                 if (K > MinLayerCell(ICell) && K < MaxLayerCell(ICell)) {
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
   parallelFor(
       "ApplyVertMixConv", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(I4 ICell, I4 KChunk) {
          ComputeVertMixConv(VertDiffOut, VertViscOut, ICell, KChunk,
                             BruntVaisalaFreqSqIn);
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
          GeomZMid(ICell, K)  = -K;
          NEdgesOnCell(ICell) = 5;
          AreaCell(ICell)     = 3.6e10_Real;
       });

   // current mesh has some CellsOnCell value > NCellsAll, so
   // filling CellsOnCell with simple mapping for this test
   parallelFor(
       "populateArrays", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          CellsOnCell(ICell, 0) = ICell;
          CellsOnCell(ICell, 1) = ICell;
          CellsOnCell(ICell, 2) = ICell;
          CellsOnCell(ICell, 3) = ICell;
          CellsOnCell(ICell, 4) = ICell;
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

/// Finalize and clean up all test infrastructure
void finalizeVertMixTest() {
   VertMix::destroyInstance();
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
   MPI_Finalize();

   // if we made it here, it is successful
   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
