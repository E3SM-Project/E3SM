//===-- Test driver for OMEGA Vertical Mixing Coefficients ---------------*- C++
//-*-===/
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
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "mpi.h"

// added for debug
#include "AuxiliaryState.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "VertCoord.h"

using namespace OMEGA;

/// Test constants and expected values
constexpr int NVertLayers = 60;

/// Values to test against
const Real VertDiffExpValueN =
    1.0391936989129011; // Expected value for diffusivity for positive BVF
const Real VertViscExpValueN =
    1.0198269656595984; // Expected value for viscosity for positive BVF
const Real VertDiffExpValueP =
    0.0015017457509864886; // Expected value for diffusivity for negative BVF
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
I4 initVertMixTest(const std::string &mesh) {

   I4 Err = 0;

   /// Initialize the Machine Environment class - this also creates
   /// the default MachEnv. Then retrieve the default environment and
   /// some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   /// Initialize logging
   initLogging(DefEnv);

   /// Open and read config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize the IO system
   IO::init(DefComm);

   /// Initialize decomposition
   Decomp::init();

   // Initialize the vertical coordinate (phase 1)
   VertCoord::init1();

   /// Initialize mesh
   HorzMesh::init();

   /// Initialize VertMix
   VertMix::init();

   /// Retrieve VertMix
   VertMix *DefVertMix = VertMix::getInstance();
   if (DefVertMix) {
      LOG_INFO("VertMixTest: VertMix retrieval PASS");
   } else {
      Err++;
      LOG_INFO("VertMixTest: VertMix retrieval FAIL");
      return -1;
   }

   return Err;
}

int testPPBackVertMix() {
   int Err             = 0;
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   /// Get VertMix instance to test
   VertMix *TestVertMix       = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NormalVelEdge =
       Array2DReal("NormalVelEdge", Mesh->NEdgesSize, VCoord->NVertLayers);
   Array2DReal TangVelEdge =
       Array2DReal("TangVelEdge", Mesh->NEdgesSize, VCoord->NVertLayers);
   Array2DReal BruntVaisalaFreqCell = Array2DReal(
       "BruntVaisalaFreqCell", Mesh->NCellsSize, VCoord->NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(BruntVaisalaFreqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) { VCoord->ZMid(ICell, K) = -K; });

   parallelFor(
       "populateArrays", {Mesh->NEdgesAll, VCoord->NVertLayers},
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
                               BruntVaisalaFreqCell);
   Array2DReal BackVertVisc = TestVertMix->VertVisc;
   Array2DReal BackVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int numMismatches = 0;
   parallelReduce(
       "CheckVertMixMatrix-BackgroundVisc",
       {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(int i, int j, int &localCount) {
          if (j == 0) {
             // Surface layer should be zero
             if (!isApprox(BackVertVisc(i, j), 0.0_Real, RTol)) {
                localCount++;
             }
          } else {
             if (!isApprox(BackVertVisc(i, j), VertViscBackExp, RTol)) {
                localCount++;
             }
          }
          return;
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixBack: VertVisc isApprox FAIL, "
                "expected {}, got {} mismatches",
                VertViscBackExp, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixBack-Visc: PASS");
   }

   /// Check total Diff against linear addition of components
   numMismatches = 0;
   parallelReduce(
       "CheckVertMixMatrix-BackgroundDiff",
       {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(int i, int j, int &localCount) {
          if (j == 0) {
             // Surface layer should be zero
             if (!isApprox(BackVertDiff(i, j), 0.0_Real, RTol)) {
                localCount++;
             }
          } else {
             if (!isApprox(BackVertDiff(i, j), VertDiffBackExp, RTol)) {
                localCount++;
             }
          }
          return;
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixBack: VertDiff isApprox FAIL, "
                "expected {}, got {} mismatches",
                VertDiffBackExp, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixBack-Diff: PASS");
   }

   return Err;
}

int testPPConvVertMix() {
   int Err             = 0;
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   /// Get VertMix instance to test
   VertMix *TestVertMix       = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NormalVelEdge =
       Array2DReal("NormalVelEdge", Mesh->NEdgesAll, VCoord->NVertLayers);
   Array2DReal TangVelEdge =
       Array2DReal("TangVelEdge", Mesh->NEdgesAll, VCoord->NVertLayers);
   Array2DReal BruntVaisalaFreqCell = Array2DReal(
       "BruntVaisalaFreqCell", Mesh->NCellsAll, VCoord->NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(BruntVaisalaFreqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) { VCoord->ZMid(ICell, K) = -K; });

   parallelFor(
       "populateArrays", {Mesh->NEdgesAll, VCoord->NVertLayers},
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
                               BruntVaisalaFreqCell);
   Array2DReal ConvVertVisc = TestVertMix->VertVisc;
   Array2DReal ConvVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int numMismatches = 0;
   parallelReduce(
       "CheckVertMixMatrix-ConvectiveVisc",
       {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(int i, int j, int &localCount) {
          if (j == 0) {
             // Surface layer should be zero
             if (!isApprox(ConvVertVisc(i, j), 0.0_Real, RTol)) {
                localCount++;
             }
          } else {
             if (!isApprox(ConvVertVisc(i, j), VertDiffConvExp, RTol)) {
                localCount++;
             }
          }
          return;
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixConv: VertVisc isApprox FAIL, "
                "expected {}, got {} mismatches",
                VertDiffConvExp, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixConv-Visc: PASS");
   }

   /// Check total Diff against linear addition of components
   numMismatches = 0;
   parallelReduce(
       "CheckVertMixMatrix-ConvectiveDiff",
       {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(int i, int j, int &localCount) {
          if (j == 0) {
             // Surface layer should be zero
             if (!isApprox(ConvVertDiff(i, j), 0.0_Real, RTol)) {
                localCount++;
             }
          } else {
             if (!isApprox(ConvVertDiff(i, j), VertDiffConvExp, RTol)) {
                localCount++;
             }
          }
          return;
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixConv: VertDiff isApprox FAIL, "
                "expected {}, got {} mismatches",
                VertDiffConvExp, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixConv-Diff: PASS");
   }

   return Err;
}

int testPPShearVertMix() {
   int Err             = 0;
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   /// Get VertMix instance to test
   VertMix *TestVertMix       = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NormalVelEdge =
       Array2DReal("NormalVelEdge", Mesh->NEdgesAll, VCoord->NVertLayers);
   Array2DReal TangVelEdge =
       Array2DReal("TangVelEdge", Mesh->NEdgesAll, VCoord->NVertLayers);
   Array2DReal BruntVaisalaFreqCell = Array2DReal(
       "BruntVaisalaFreqCell", Mesh->NCellsAll, VCoord->NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);
   deepCopy(BruntVaisalaFreqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          VCoord->ZMid(ICell, K)    = -K;
          Mesh->NEdgesOnCell(ICell) = 5;
          Mesh->AreaCell(ICell)     = 3.6e10_Real;
       });

   parallelFor(
       "populateArrays", {Mesh->NEdgesAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
          TangVelEdge(IEdge, K)   = TangVelEdge(IEdge, K) + 0.5 * K;
          Mesh->DcEdge(IEdge)     = 2.0e5_Real;
          Mesh->DvEdge(IEdge)     = 1.45e5_Real;
       });

   /// Compute only shear vertical viscosity and diffusivity
   TestVertMix->BackDiff                    = 0.0;
   TestVertMix->BackVisc                    = 0.0;
   TestVertMix->ComputeVertMixConv.Enabled  = false;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqCell);
   Array2DReal ShearVertVisc = TestVertMix->VertVisc;
   Array2DReal ShearVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int numMismatches = 0;
   Kokkos::fence();
   auto ShearVertViscH = createHostMirrorCopy(ShearVertVisc);
   parallelReduce(
       "CheckVertMixMatrix-ShearVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int i, int &localCount) {
          // Surface layer should be zero
          if (!isApprox(ShearVertVisc(i, 0), 0.0_Real, RTol)) {
             localCount++;
          }
          if (!isApprox(ShearVertVisc(i, 1), VertViscShearExp, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixShear: VertVisc isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertViscShearExp, ShearVertViscH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixShear-Visc: PASS");
   }

   /// Check total Diff against linear addition of components
   numMismatches       = 0;
   auto ShearVertDiffH = createHostMirrorCopy(ShearVertDiff);
   parallelReduce(
       "CheckVertMixMatrix-ShearDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int i, int &localCount) {
          // Surface layer should be zero
          if (!isApprox(ShearVertDiff(i, 0), 0.0_Real, RTol)) {
             localCount++;
          }
          if (!isApprox(ShearVertDiff(i, 1), VertDiffShearExp, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixShear: VertDiff isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertDiffShearExp, ShearVertDiffH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixShear-Diff: PASS");
   }

   return Err;
}

/// Test vertical mixing coefficients calculation for all cells/layers
int testPPTotalVertMix() {
   int Err             = 0;
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   /// Get VertMix instance to test
   VertMix *TestVertMix       = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NormalVelEdge =
       Array2DReal("NormalVelEdge", Mesh->NEdgesAll, VCoord->NVertLayers);
   Array2DReal TangVelEdge =
       Array2DReal("TangVelEdge", Mesh->NEdgesAll, VCoord->NVertLayers);
   Array2DReal BruntVaisalaFreqCell = Array2DReal(
       "BruntVaisalaFreqCell", Mesh->NCellsAll, VCoord->NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NormalVelEdge, NV);
   deepCopy(TangVelEdge, TV);

   // Test with positive BVF first
   deepCopy(BruntVaisalaFreqCell, BVFP);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
       "populateArrays", {Mesh->NCellsAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          VCoord->ZMid(ICell, K)    = -K;
          Mesh->NEdgesOnCell(ICell) = 5;
          Mesh->AreaCell(ICell)     = 3.6e10_Real;
       });

   parallelFor(
       "populateArrays", {Mesh->NEdgesAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(I4 IEdge, I4 K) {
          NormalVelEdge(IEdge, K) = NormalVelEdge(IEdge, K) + 0.5 * K;
          TangVelEdge(IEdge, K)   = TangVelEdge(IEdge, K) + 0.5 * K;
          Mesh->DcEdge(IEdge)     = 2.0e5_Real;
          Mesh->DvEdge(IEdge)     = 1.45e5_Real;
       });

   /// Compute vertical viscosity and diffusivity
   TestVertMix->BackDiff                    = 1.0e-5;
   TestVertMix->BackVisc                    = 1.0e-4;
   TestVertMix->ComputeVertMixConv.Enabled  = true;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqCell);
   Array2DReal VertDiffP = TestVertMix->VertDiff;
   Array2DReal VertViscP = TestVertMix->VertVisc;

   /// Check all VertDiff array values against expected value
   int numMismatches = 0;
   auto VertDiffPH   = createHostMirrorCopy(VertDiffP);
   parallelReduce(
       "CheckVertMixMatrix-TotalPosDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int i, int &localCount) {
          // Surface layer should be zero
          if (!isApprox(VertDiffP(i, 0), 0.0_Real, RTol)) {
             localCount++;
          }
          if (!isApprox(VertDiffP(i, 1), VertDiffExpValueP, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixTotal: VertDiffPositive isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertDiffExpValueP, VertDiffPH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMixTotal VertDiffPositive: PASS");
   }

   /// Check all VertVisc array values against expected value
   numMismatches   = 0;
   auto VertViscPH = createHostMirrorCopy(VertViscP);
   parallelReduce(
       "CheckVertMixMatrix-TotalPosVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int i, int &localCount) {
          // Surface layer should be zero
          if (!isApprox(VertViscP(i, 0), 0.0_Real, RTol)) {
             localCount++;
          }
          if (!isApprox(VertViscP(i, 1), VertViscExpValueP, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixTotal: VertViscPositive isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertViscExpValueP, VertViscPH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMixTotal VertViscPositive: PASS");
   }

   // Now test with negative BVF
   deepCopy(BruntVaisalaFreqCell, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   /// Compute vertical viscosity and diffusivity
   TestVertMix->computeVertMix(NormalVelEdge, TangVelEdge,
                               BruntVaisalaFreqCell);
   Array2DReal VertDiffN = TestVertMix->VertDiff;
   Array2DReal VertViscN = TestVertMix->VertVisc;

   /// Check all VertDiff array values against expected value
   numMismatches   = 0;
   auto VertDiffNH = createHostMirrorCopy(VertDiffN);
   parallelReduce(
       "CheckVertMixMatrix-TotalNegDiff", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int i, int &localCount) {
          // Surface layer should be zero
          if (!isApprox(VertDiffN(i, 0), 0.0_Real, RTol)) {
             localCount++;
          }
          if (!isApprox(VertDiffN(i, 1), VertDiffExpValueN, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMix: VertDiffNegative isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertDiffExpValueN, VertDiffNH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMix VertDiffNegative: PASS");
   }

   /// Check all VertVisc array values against expected value
   numMismatches   = 0;
   auto VertViscNH = createHostMirrorCopy(VertViscN);
   parallelReduce(
       "CheckVertMixMatrix-TotalNegVisc", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int i, int &localCount) {
          // Surface layer should be zero
          if (!isApprox(VertViscN(i, 0), 0.0_Real, RTol)) {
             localCount++;
          }
          if (!isApprox(VertViscN(i, 1), VertViscExpValueN, RTol)) {
             localCount++;
          }
       },
       numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMix: VertViscNegative isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertViscExpValueN, VertViscNH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMix VertViscNegative : PASS");
   }

   return Err;
}

/// Finalize and clean up all test infrastructure
void finalizeVertMixTest() {
   HorzMesh::clear();
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
int vertMixTest(const std::string &MeshFile = "OmegaMesh.nc") {
   int Err = initVertMixTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("VertMix: Error initializing");
   }
   const auto &Mesh = HorzMesh::getDefault();

   Err += testPPBackVertMix();
   Err += testPPConvVertMix();
   Err += testPPShearVertMix();
   Err += testPPTotalVertMix();

   if (Err == 0) {
      LOG_INFO("VertMix: Successful completion");
   }
   finalizeVertMixTest();
   return Err;
}

// The test driver for VertMix testing
int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += vertMixTest();

   VertMix::destroyInstance();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/

// Err += initState(LayerThickCell, NormalVelEdge, TangVelEdge,
//    BruntVaisalaFreqCell, Mesh);

// int initState(const Array2DReal &LayerThickCell,
//               const Array2DReal &NormalVelEdge,
//               const Array2DReal &TangVelEdge,
//               const Array2DReal &BruntVaisalaFreqCell,
//               const HorzMesh *Mesh) {
//    int Err = 0;

//   TestSetup Setup;

//   Err += setScalar(
//      KOKKOS_LAMBDA(Real X, Real Y, I4 K) { return Setup.layerThickness(X, Y);
//      }, LayerThickCell, Geom, Mesh, OnCell);

//   Err += setVectorEdge(
//      KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
//         VecField[0] = Setup.velocityX(X, Y);
//         VecField[1] = Setup.velocityY(X, Y);
//      },
//      NormalVelEdge, EdgeComponent::Normal, Geom, Mesh);

// Err += setVectorEdge(
//     KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y, I4 K) {
//        VecField[0] = NV + 0.5 * K;
//        VecField[1] = 0.0;
//     },
//     NormalVelEdge, EdgeComponent::Normal, Geom, Mesh);
