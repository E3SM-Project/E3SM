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
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "VertCoord.h"
#include "mpi.h"

using namespace OMEGA;

/// Test constants and expected values
constexpr int NVertLayers = 60;

/// Values to test against
const Real VertDiffExpValueN =
    1.0393923290498872; // Expected value for diffusivity for positive BVF
const Real VertViscExpValueN =
    1.0198269656595984; // Expected value for viscosity for positive BVF
const Real VertDiffExpValueP =
    0.0015685660274841844; // Expected value for diffusivity for negative BVF
const Real VertViscExpValueP =
    0.002332474675614262; // Expected value for viscosity for negative BVF
const Real VertDiffBackExp =
    1.0e-5; // Expected value for background diffusivity
const Real VertViscBackExp = 1.0e-4; // Expected value for background viscosity
const Real VertDiffConvExp =
    1.0; // Expected value for convective diffusivity/viscosity
const Real VertDiffShearExp =
    0.039183698912901; // Expected value for shear diffusivity
const Real VertViscShearExp =
    0.01972696565959843; // Expected value for shear viscosity

/// Test input values
const Real BVFP = 0.1;  // Positive Brunt-Vaisala frequency in s^-2
const Real BVFN = -0.1; // Negative Brunt-Vaisala frequency in s^-2
const Real NV   = 1.0;  // Normal velocity in m/s
const Real TV   = 1.0;  // Tangential velocity in m/s
const Real RTol = 1e-8; // Relative tolerance for isApprox checks

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

   /// Initialize mesh
   HorzMesh::init();

   /// Initialize vertical coordinate
   VertCoord::init(false);

   /// Initialize VertMix
   VertMix::init();

   /// Retrieve VertMix
   VertMix *DefVertMix = VertMix::getInstance();
   if (!DefVertMix)
      ABORT_ERROR("VertMixTest: VertMix retrieval FAIL");
}

void testBackVertMix() {
   // Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   I4 NEdgesSize       = Mesh->NEdgesSize;
   I4 NEdgesAll        = Mesh->NEdgesAll;
   OMEGA_SCOPE(ZMid, VCoord->ZMid);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesSize, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesSize, NVertLayers);
   auto BruntVaisalaFreqSqCell =
       Array2DReal("BruntVaisalaFreqSqCell", NCellsSize, NVertLayers);

   /// Use deep copy initialize with reference or zero values
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(BruntVaisalaFreqSqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) { ZMid(ICell, K) = -K; });

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

   /// Check total Visc against linear addition of components
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-BackgroundVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    // Surface layer should be zero
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

   if (NumMismatches != 0)
      ABORT_ERROR("TestVertMixBack: VertVisc FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertViscBackExp, BackVertVisc(1, 1), NumMismatches);

   /// Check total Diff against linear addition of components
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-BackgroundDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    // Surface layer should be zero
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
   }

   return;
}

void testConvVertMix() {

   // Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   I4 NEdgesAll        = Mesh->NEdgesAll;
   OMEGA_SCOPE(ZMid, VCoord->ZMid);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesAll, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesAll, NVertLayers);
   auto BruntVaisalaFreqSqCell =
       Array2DReal("BruntVaisalaFreqSqCell", NCellsSize, NVertLayers);

   /// Use deep copy to initialize with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(BruntVaisalaFreqSqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) { ZMid(ICell, K) = -K; });

   parallelFor(
       "populateArrays", {NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
          TangVelEdge(IEdge, K)   = TangVelEdge(IEdge, K) + 0.5 * K;
       });

   /// Compute only convective vertical viscosity and diffusivity
   TestVertMix->BackDiff                    = 0.0;
   TestVertMix->BackVisc                    = 0.0;
   TestVertMix->ComputeVertMixConv.Enabled  = true;
   TestVertMix->ComputeVertMixShear.Enabled = false;
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   Array2DReal ConvVertVisc = TestVertMix->VertVisc;
   Array2DReal ConvVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ConvectiveVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    // Surface layer should be zero
                    if (ConvVertVisc(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else {
                    if (!isApprox(ConvVertVisc(ICell, K), VertDiffConvExp,
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
      auto ConvVertViscH = createHostMirrorCopy(ConvVertVisc);
      ABORT_ERROR("TestVertMixConv: VertVisc FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertDiffConvExp, ConvVertViscH(1, 1), NumMismatches);
   }

   /// Check total Diff against linear addition of components
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ConvectiveDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    // Surface layer should be zero
                    if (ConvVertDiff(ICell, K) != 0.0_Real)
                       InnerCount++;
                 } else {
                    if (!isApprox(ConvVertDiff(ICell, K), VertDiffConvExp,
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
      auto ConvVertDiffH = createHostMirrorCopy(ConvVertDiff);
      ABORT_ERROR("TestVertMixConv: VertDiff FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertDiffConvExp, ConvVertDiffH(1, 1), NumMismatches);
   }

   return;
}

void testShearVertMix() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   I4 NEdgesAll        = Mesh->NEdgesAll;
   OMEGA_SCOPE(ZMid, VCoord->ZMid);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesAll, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesAll, NVertLayers);
   auto BruntVaisalaFreqSqCell =
       Array2DReal("BruntVaisalaFreqSqCell", NCellsSize, NVertLayers);

   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(BruntVaisalaFreqSqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          ZMid(ICell, K)      = -K;
          NEdgesOnCell(ICell) = 5;
          AreaCell(ICell)     = 3.6e10_Real;
       });

   parallelFor(
       "populateArrays", {NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
          TangVelEdge(IEdge, K)   = TangVelEdge(IEdge, K) + 0.5 * K;
          DcEdge(IEdge)           = 2.0e5_Real;
          DvEdge(IEdge)           = 1.45e5_Real;
       });

   /// Compute only shear vertical viscosity and diffusivity
   TestVertMix->BackDiff                    = 0.0;
   TestVertMix->BackVisc                    = 0.0;
   TestVertMix->ComputeVertMixConv.Enabled  = false;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqSqCell);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   Array2DReal ShearVertVisc = TestVertMix->VertVisc;
   Array2DReal ShearVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ShearVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    if (ShearVertVisc(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == 1) {
                    if (!isApprox(ShearVertVisc(ICell, K), VertViscShearExp,
                                  RTol))
                       InnerCount++;
                    // otherwise check for invalid values
                 } else {
                    if (ShearVertVisc(ICell, K) == 0.0 or
                        Kokkos::isnan(ShearVertVisc(ICell, K)) or
                        Kokkos::isinf(ShearVertVisc(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      auto ShearVertViscH = createHostMirrorCopy(ShearVertVisc);
      ABORT_ERROR("TestVertMixShear: VertVisc FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertViscShearExp, ShearVertViscH(1, 1), NumMismatches);
   }

   /// Check total Diff against linear addition of components
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-ShearVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    if (ShearVertDiff(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == 1) {
                    if (!isApprox(ShearVertDiff(ICell, K), VertDiffShearExp,
                                  RTol))
                       InnerCount++;
                    // otherwise check for invalid values
                 } else {
                    if (ShearVertDiff(ICell, K) == 0.0 or
                        Kokkos::isnan(ShearVertDiff(ICell, K)) or
                        Kokkos::isinf(ShearVertDiff(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      auto ShearVertDiffH = createHostMirrorCopy(ShearVertDiff);
      ABORT_ERROR("TestVertMixShear: VertDiff FAIL, "
                  "expected {}, got {} with {} mismatches",
                  VertDiffShearExp, ShearVertDiffH(1, 1), NumMismatches);
   }

   return;
}

/// Test vertical mixing coefficients calculation for all cells/layers
void testTotalVertMix() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   I4 NEdgesAll        = Mesh->NEdgesAll;
   OMEGA_SCOPE(ZMid, VCoord->ZMid);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);

   /// Get VertMix instance to test
   VertMix *TestVertMix = VertMix::getInstance();

   /// Create and fill ocean state arrays
   auto NormalVelEdge = Array2DReal("NormalVelEdge", NEdgesAll, NVertLayers);
   auto TangVelEdge   = Array2DReal("TangVelEdge", NEdgesAll, NVertLayers);
   auto BruntVaisalaFreqSqCell =
       Array2DReal("BruntVaisalaFreqSqCell", NCellsSize, NVertLayers);

   /// Use deep copy to initialize with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);

   // Test with positive BVF first
   deepCopy(BruntVaisalaFreqSqCell, BVFP);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          ZMid(ICell, K)      = -K;
          NEdgesOnCell(ICell) = 5;
          AreaCell(ICell)     = 3.6e10_Real;
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
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    if (VertDiffP(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == 1) {
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
   }

   /// Check all VertVisc array values against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-TotalPosVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    if (VertViscP(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == 1) {
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
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    if (VertDiffN(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == 1) {
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
   }

   /// Check all VertVisc array values against expected value
   NumMismatches = 0;
   parallelReduceOuter(
       "CheckVertMixMatrix-TotalNegVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) {
                    if (VertViscN(ICell, K) != 0.0_Real)
                       InnerCount++;
                    // K = 1 should have ref value
                 } else if (K == 1) {
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
// --> one tests the vertical diffusivity and viscosity
// with only background on
// --> next tests the vertical diffusivity and viscosity
// with only convective on
// --> next tests the vertical diffusivity and viscosity
// with only shear on
// --> next tests the linear superposition of the
// background, convective, and shear contributions
void vertMixTest() {

   // initialize vertical mix and other infrastructure
   initVertMixTest();

   // test each vertical mix option
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
