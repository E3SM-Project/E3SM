//===-- Test driver for OMEGA GSW-C library -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA GSW-C external library
///
/// This driver tests that the GSW-C library can be called
/// and returns expected value (as published in Roquet et al 2015)
//
//===-----------------------------------------------------------------------===/

#include "Eos.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "IO.h"
#include "IOStream.h"
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

#include <gswteos-10.h>

using namespace OMEGA;

/// Test constants and expected values
constexpr int NVertLayers = 60;

/// Published values (TEOS-10 and linear) to test against
const Real TeosSVExpValue =
    0.0009732819628; // Expected value for TEOS-10 specific volume
const Real LinearExpValue =
    0.0009784735812133072; // Expected value for Linear specific volume
const Real ConstantExpValue =
    1.0_Real / RhoSw; // Expected value for constant specific volume
const Real TeosBVFExpValue =
    0.020913834194283325; // Expected value for TEOS-10 squared Brunt-Vaisala
                          // frequency
const Real LinearBVFExpValue =
    0.017834796542017275; // Expected value for Linear squared Brunt-Vaisala
                          // frequency
const Real GswBVFExpValue =
    0.02081197958166906; // Expected value from GSW-C library

/// Test input values
const Real Sa = 30.0;           // Absolute Salinity in g/kg
const Real Ct = 10.0;           // Conservative Temperature in degC
const Real P  = 1000.0 * Db2Pa; // Pressure in Pa

const I4 KDisp  = 1;     // Displate parcel to K=1 for TEOS-10 eos
const Real RTol = 1e-10; // Relative tolerance for isApprox checks

/// The initialization routine for Eos testing. It calls various
/// init routines, including the creation of the default decomposition.
void initEosTest(const std::string &mesh) {

   /// Initialize the Machine Environment class - this also creates
   /// the default MachEnv. Then retrieve the default environment and
   /// some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   /// Initialize logging
   initLogging(DefEnv);
   LOG_INFO("------ EOS Unit Tests ------");

   /// Open and read config file
   Config("Omega");
   Config::readAll("omega.yml");

   /// Initialize parallel IO
   IO::init(DefComm);

   /// Initialize decomposition
   Decomp::init(mesh);

   /// Initialize Halo
   Halo::init();

   /// Create dummy model clock
   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   Clock ModelClockTmp(StartTime, TimeStep);
   Clock *ModelClock = &ModelClockTmp;

   /// Read horizontal mesh
   Field::init(ModelClock);
   IOStream::init(ModelClock);
   HorzMesh::init(ModelClock);

   /// Initialize vertical coordinate
   VertCoord::init(false);

   /// Initialize Eos
   Eos::init();

   /// Retrieve Eos
   Eos *DefEos = Eos::getInstance();
   if (!DefEos)
      ABORT_ERROR("EosTest: Eos retrieval FAIL");
}

/// Test Linear EOS calculation for all cells/layers
void testEosLinear() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   /// Compute specific volume
   TestEos->computeSpecVol(TArray, SArray, PArray);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches   = 0;
   Array2DReal SpecVol = TestEos->SpecVol;
   parallelReduceOuter(
       "CheckSpecVolMatrix-linear", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), LinearExpValue, RTol)) {
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
      auto SpecVolH = createHostMirrorCopy(TestEos->SpecVol);
      for (int I = 0; I < Mesh->NCellsAll; ++I) {
         for (int K = 0; K < NVertLayers; ++K) {
            if (!isApprox(SpecVolH(I, K), LinearExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol Linear Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         I, K, SpecVolH(I, K), LinearExpValue);
         }
      }
      ABORT_ERROR("EosTest: SpecVol Linear FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Test Linear EOS calculation with vertical displacement
void testEosLinearDisplaced() {
   /// Get mesh and coord info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVolDisplaced, 0.0);

   /// Compute displaced specific volume
   TestEos->computeSpecVolDisp(TArray, SArray, PArray, KDisp);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches            = 0;
   Array2DReal SpecVolDisplaced = TestEos->SpecVolDisplaced;
   parallelReduceOuter(
       "CheckSpecVolDispMatrix-linear", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVolDisplaced(ICell, K), LinearExpValue,
                               RTol)) {
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
      auto SpecVolDisplacedH = createHostMirrorCopy(SpecVolDisplaced);
      for (int I = 0; I < Mesh->NCellsAll; ++I) {
         for (int K = 0; K < NVertLayers; ++K) {
            if (!isApprox(SpecVolDisplacedH(I, K), LinearExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol Linear Displaced Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         I, K, SpecVolDisplacedH(I, K), LinearExpValue);
         }
      }
      ABORT_ERROR("EosTest: Linear SpecVolDisp FAIL with {} bad values ",
                  NumMismatches);
   }

   return;
}

/// Test Constant EOS calculation for all cells/layers
void testEosConstant() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsAll        = Mesh->NCellsAll;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::ConstantEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsAll, NVertLayers);
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   /// Compute specific volume
   TestEos->computeSpecVol(TArray, SArray, PArray);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all active layers against expected constant value
   int NumMismatches   = 0;
   Array2DReal SpecVol = TestEos->SpecVol;
   parallelReduceOuter(
       "CheckSpecVolMatrix-Constant", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), ConstantExpValue, RTol)) {
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
      auto SpecVolH = createHostMirrorCopy(SpecVol);
      for (int I = 0; I < NCellsAll; ++I) {
         for (int K = 0; K < NVertLayers; ++K) {
            if (!isApprox(SpecVolH(I, K), ConstantExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol Constant Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         I, K, SpecVolH(I, K), ConstantExpValue);
         }
      }
      ABORT_ERROR("EosTest: SpecVol Constant FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Test linear squared Brunt-Vaisala frequency calculation for all cells/layers
void testBruntVaisalaFreqSqLinear() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   /// Use deep copy to initialize results to zero
   deepCopy(TestEos->SpecVol, 0.0);
   deepCopy(TestEos->BruntVaisalaFreqSq, 0.0);

   // fill remaining entries with sample values that should lead to ref result
   // for K = 1.
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);
   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          if (K == 0) {
             GeomZMid(ICell, 0) = -992.1173890198451_Real;
             SArray(ICell, 0)   = Sa - 1.0_Real;
             TArray(ICell, 0)   = Ct + 15.0_Real;
             PArray(ICell, 0)   = P;
          } else if (K == 1) {
             GeomZMid(ICell, 1) = -993.1071379053125_Real;
             SArray(ICell, 1)   = Sa;
             TArray(ICell, 1)   = Ct + 10.0_Real;
             PArray(ICell, 1)   = P + 1.0_Real;
          } else if (K == 2) {
             GeomZMid(ICell, 2) = -994.0968821072275_Real;
             SArray(ICell, 2)   = Sa + 1.0_Real;
             TArray(ICell, 2)   = Ct + 5.0_Real;
             PArray(ICell, K)   = P + 2.0_Real;
          } else { // fill rest to valid junk to avoid NaNs or Inf
             GeomZMid(ICell, K) = -994.0968821072275_Real - 0.1_Real * K;
             SArray(ICell, K)   = Sa + 1.0_Real + 0.1_Real * K;
             TArray(ICell, K)   = Ct + 5.0_Real - 0.01_Real * K;
             PArray(ICell, K)   = P + 2.0_Real + 0.1_Real * K;
          }
       });

   /// Compute specific volume first
   TestEos->computeSpecVol(TArray, SArray, PArray);
   Array2DReal SpecVol = TestEos->SpecVol;

   /// Compute squared Brunt-Vaisala frequency
   TestEos->computeBruntVaisalaFreqSq(TArray, SArray, PArray, SpecVol);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches = 0;
   OMEGA_SCOPE(BruntVaisalaFreqSq, TestEos->BruntVaisalaFreqSq);
   parallelReduceOuter(
       "CheckBruntVaisalaSq-Linear", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) { // should be zero
                    if (BruntVaisalaFreqSq(ICell, K) != 0.0)
                       InnerCount++;
                 } else if (K == 1) { // should be ref value
                    if (!isApprox(BruntVaisalaFreqSq(ICell, K),
                                  LinearBVFExpValue, RTol))
                       InnerCount++;
                 } else { // just check for unreasonable values
                    if (BruntVaisalaFreqSq(ICell, K) == 0.0 or
                        Kokkos::isnan(BruntVaisalaFreqSq(ICell, K)) or
                        Kokkos::isinf(BruntVaisalaFreqSq(ICell, K)))
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
      auto BruntVaisalaFreqSqH = createHostMirrorCopy(BruntVaisalaFreqSq);
      for (int I = 0; I < Mesh->NCellsAll; ++I) {
         // top layer should be zero
         if (BruntVaisalaFreqSqH(I, 0) != 0.0)
            LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      I, 0, BruntVaisalaFreqSqH(I, 0), 0.0);
         // K = 1 should be ref value
         if (!isApprox(BruntVaisalaFreqSqH(I, 1), LinearBVFExpValue, RTol))
            LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      I, 1, BruntVaisalaFreqSqH(I, 1), LinearBVFExpValue);
         // remaining values just check for other conditions
         for (int K = 2; K < NVertLayers; ++K) {
            if (BruntVaisalaFreqSqH(I, K) == 0.0 or
                Kokkos::isnan(BruntVaisalaFreqSqH(I, K)) or
                Kokkos::isinf(BruntVaisalaFreqSqH(I, K)))
               LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                         "BruntVaisala({},{}) = {}",
                         I, K, BruntVaisalaFreqSqH(I, K));
         }
      }
      ABORT_ERROR("EosTest: BruntVaisala Linear FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Test TEOS-10 EOS calculation for all cells/layers
void testEosTeos10() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   /// Compute specific volume
   TestEos->computeSpecVol(TArray, SArray, PArray);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches   = 0;
   Array2DReal SpecVol = TestEos->SpecVol;
   parallelReduceOuter(
       "CheckSpecVolMatrix-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), TeosSVExpValue, RTol)) {
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
      auto SpecVolH = createHostMirrorCopy(SpecVol);
      for (int I = 0; I < Mesh->NCellsAll; ++I) {
         for (int K = 0; K < NVertLayers; ++K) {
            if (!isApprox(SpecVolH(I, K), LinearExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol TEOS Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         I, K, SpecVolH(I, K), LinearExpValue);
         }
      }
      ABORT_ERROR("EosTest: SpecVol TEOS FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Test TEOS-10 EOS calculation with vertical displacement
void testEosTeos10Displaced() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVolDisplaced, 0.0);

   /// Compute displaced specific volume
   TestEos->computeSpecVolDisp(TArray, SArray, PArray, KDisp);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches            = 0;
   Array2DReal SpecVolDisplaced = TestEos->SpecVolDisplaced;
   parallelReduceOuter(
       "CheckSpecVolDispMatrix-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVolDisplaced(ICell, K), TeosSVExpValue,
                               RTol)) {
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
      auto SpecVolDisplacedH = createHostMirrorCopy(SpecVolDisplaced);
      for (int I = 0; I < Mesh->NCellsAll; ++I) {
         for (int K = 0; K < NVertLayers; ++K) {
            if (!isApprox(SpecVolDisplacedH(I, K), LinearExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol Displaced TEOS Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         I, K, SpecVolDisplacedH(I, K), LinearExpValue);
         }
      }
      ABORT_ERROR("EosTest: SpecVol Displaced TEOS FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Test TEOS-10 squared Brunt-Vaisala frequency calculation for all cells/layer
void testBruntVaisalaFreqSqTeos10() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   /// Use deep copy to initialize results to zero
   deepCopy(TestEos->BruntVaisalaFreqSq, 0.0);
   deepCopy(TestEos->SpecVol, 0.0);

   /// Fill inputs with values that should lead to ref result for K=1
   OMEGA_SCOPE(GeomZMid, VCoord->GeomZMid);
   parallelFor(
       "populateArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          if (K == 0) {
             GeomZMid(ICell, 0) = -992.1173890198451_Real;
             SArray(ICell, 0)   = Sa - 1.0_Real;
             TArray(ICell, 0)   = Ct + 15.0_Real;
             PArray(ICell, 0)   = P;
          } else if (K == 1) {
             GeomZMid(ICell, 1) = -993.1071379053125_Real;
             SArray(ICell, 1)   = Sa;
             TArray(ICell, 1)   = Ct + 10.0_Real;
             PArray(ICell, 1)   = (P * Pa2Db + 1.0_Real) * Db2Pa;
          } else if (K == 2) {
             GeomZMid(ICell, 2) = -994.0968821072275_Real;
             SArray(ICell, 2)   = Sa + 1.0_Real;
             TArray(ICell, 2)   = Ct + 5.0_Real;
             PArray(ICell, K)   = (P * Pa2Db + 2.0_Real) * Db2Pa;
          } else { // fill rest with valid junk to avoid Nans and Inf
             GeomZMid(ICell, K) = -994.0968821072275_Real - 0.1_Real * K;
             SArray(ICell, K)   = Sa + 1.0_Real + 0.1_Real * K;
             TArray(ICell, K)   = Ct + 5.0_Real - 0.01_Real * K;
             PArray(ICell, K)   = (P * Pa2Db + 2.0_Real + 0.1_Real * K) * Db2Pa;
          }
       });

   /// Compute specific volume first
   TestEos->computeSpecVol(TArray, SArray, PArray);
   Array2DReal SpecVol = TestEos->SpecVol;

   /// Compute Brunt-Vaisala frequency
   TestEos->computeBruntVaisalaFreqSq(TArray, SArray, PArray, SpecVol);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches = 0;
   OMEGA_SCOPE(BruntVaisalaFreqSq, TestEos->BruntVaisalaFreqSq);
   parallelReduceOuter(
       "CheckSpecVolMatrix-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 0) { // should be zero at top
                    if (BruntVaisalaFreqSq(ICell, K) != 0.0)
                       InnerCount++;
                 } else if (K == 1) { // should be ref value
                    if (!isApprox(BruntVaisalaFreqSq(ICell, K), TeosBVFExpValue,
                                  RTol))
                       InnerCount++;
                 } else { // just check for unreasonable values
                    if (BruntVaisalaFreqSq(ICell, K) == 0.0 or
                        Kokkos::isnan(BruntVaisalaFreqSq(ICell, K)) or
                        Kokkos::isinf(BruntVaisalaFreqSq(ICell, K)))
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
      auto BruntVaisalaFreqSqH = createHostMirrorCopy(BruntVaisalaFreqSq);
      for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         // top layer should be zero
         if (BruntVaisalaFreqSqH(ICell, 0) != 0.0)
            LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      ICell, 0, BruntVaisalaFreqSqH(ICell, 0), 0.0);
         // K = 1 should be ref value
         if (!isApprox(BruntVaisalaFreqSqH(ICell, 1), TeosBVFExpValue, RTol))
            LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      ICell, 1, BruntVaisalaFreqSqH(ICell, 1), TeosBVFExpValue);
         // remaining values just check for other conditions
         for (int K = 2; K < NVertLayers; ++K) {
            if (BruntVaisalaFreqSqH(ICell, K) == 0.0 or
                Kokkos::isnan(BruntVaisalaFreqSqH(ICell, K)) or
                Kokkos::isinf(BruntVaisalaFreqSqH(ICell, K)))
               LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                         "BruntVaisala({},{}) = {}",
                         ICell, K, BruntVaisalaFreqSqH(ICell, K));
         }
      }
      ABORT_ERROR("EosTest: BruntVaisala TEOS FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Finalize and clean up all test infrastructure
void finalizeEosTest() {
   Eos::destroyInstance();
   VertCoord::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

/// Test that the external GSW-C library returns the expected specific volume
void checkValueGswcSpecVol() {
   const Real RTol = 1e-10;

   /// Get specific volume from GSW-C library
   double SpecVol = gsw_specvol(Sa, Ct, P * Pa2Db);
   /// Check the value against the expected TEOS-10 value
   bool Check = isApprox(SpecVol, TeosSVExpValue, RTol);
   if (!Check) {
      ABORT_ERROR("checkValueGswcSpecVol: SpecVol FAIL, expected {}, got {}",
                  TeosSVExpValue, SpecVol);
   }
   return;
}

/// Test that the external GSW-C library returns the expected N2
void checkValueGswcN2() {
   const Real RTol = 1e-10;

   // Number of intervals (nz)
   int Nz = 2;

   // Input arrays: length nz+1
   double Salt[4]  = {Sa - 1.0, Sa, Sa + 1.0}; // Absolute Salinity (g/kg)
   double Temp[4]  = {Ct + 15.0, Ct + 10.0,
                      Ct + 5.0}; // Conservative Temperature (deg C)
   double Press[4] = {P * Pa2Db, P * Pa2Db + 1.0,
                      P * Pa2Db + 2.0}; // Pressure (dbar)

   // Latitude (degrees north)
   double Latitude[4] = {0.0, 0.0, 0.0};

   // Output arrays: length nz
   double N2[Nz];   // Brunt–Väisälä frequency squared
   double PMid[Nz]; // Midpoint pressure

   /// Get specific volume from GSW-C library
   gsw_nsquared(Salt, Temp, Press, Latitude, Nz, N2, PMid);

   /// Check the value against the expected TEOS-10 value
   bool Check = isApprox(N2[0], GswBVFExpValue, RTol);
   if (!Check) {
      ABORT_ERROR("checkValueGswcN2: N2 FAIL, expected {}, got {}",
                  GswBVFExpValue, N2[0]);
   }
   return;
}

// the main tests (all in one to have the same log):
// Single value test:
// --> test calls the external GSW-C library
// and compares the specific volume to the published value
// Full array tests:
// --> one tests the value on a Eos with linear option
// --> next checks the value on a Eos with linear displaced option
// --> next checks the value of the linear squared Brunt Vaisala Freq.
// calculation
// --> next checks the value on a Eos with TEOS-10 option
// --> next checks the value on a Eos with TEOS-10 displaced option
// --> last checks the value of the TOES-10 squared Brunt Vaisala Freq.
// calculation
void eosTest(const std::string &MeshFile = "OmegaMesh.nc") {
   initEosTest(MeshFile);
   const auto &Mesh = HorzMesh::getDefault();

   checkValueGswcSpecVol();
   checkValueGswcN2();

   testEosLinear();
   testEosLinearDisplaced();
   testEosConstant();
   testBruntVaisalaFreqSqLinear();
   testEosTeos10();
   testEosTeos10Displaced();
   testBruntVaisalaFreqSqTeos10();

   finalizeEosTest();

   return;
}

// The test driver for Eos testing
int main(int argc, char *argv[]) {

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   eosTest();

   LOG_INFO("------ EOS Unit Tests Successful ------");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   // If we made it here, test is successful
   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
