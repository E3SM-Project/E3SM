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
    0.017382842633545097; // Expected value for Linear squared Brunt-Vaisala
                          // frequency
const Real GswBVFExpValue =
    0.02081197958166906; // Expected value from GSW-C library

/// Linear EOS coefficients, matching the Linear subsection of the Eos group in
/// the test configuration. The expected derivatives follow from
/// SpecVol = 1 / (RhoT0S0 + DRhodT * Ct + DRhodS * Sa).
const Real LinearDRhodT      = -0.2;
const Real LinearDRhodS      = 0.8;
const Real LinearDCtExpValue = -LinearDRhodT * LinearExpValue * LinearExpValue;
const Real LinearDSaExpValue = -LinearDRhodS * LinearExpValue * LinearExpValue;

/// Test input values
const Real Sa = 30.0;           // Absolute Salinity in g/kg
const Real Ct = 10.0;           // Conservative Temperature in degC
const Real P  = 1000.0 * Db2Pa; // Pressure in Pa

const I4 KDisp  = 1;     // Displate parcel to K=1 for TEOS-10 eos
const Real RTol = 1e-10; // Relative tolerance for isApprox checks

/// States spanning the oceanographic range and its corners, used by the
/// specific volume derivative checks. The fresh end is included because that
/// is where the normalized salinity is smallest and the salinity derivative
/// worst conditioned, and the full pressure range because the reference
/// profile supplies almost all of the pressure derivative.
constexpr int NSaTest      = 6;
constexpr int NCtTest      = 6;
constexpr int NPTest       = 6;
const Real SaTest[NSaTest] = {0.0, 5.0, 20.0, 30.0, 35.0, 38.5}; // g/kg
const Real CtTest[NCtTest] = {-2.0, 0.0, 4.0, 10.0, 25.0, 35.0}; // degC
const Real PTest[NPTest]   = {0.0,    100.0,  1000.0,
                              4000.0, 8000.0, 10000.0}; // dbar

/// Relative tolerance for the specific volume derivative checks against the
/// GSW-C library. Omega and GSW-C evaluate the same polynomial in different
/// arrangements, so a few ulp of disagreement is expected; anything larger
/// than this means a term has been dropped or mis-scaled.
const Real DerivRTol = 1e-12;

/// The pressure derivative needs a looser tolerance, and the reason is on the
/// GSW-C side rather than ours. GSW-C evaluates its v_P from a table of
/// coefficients that have been pre-multiplied by their pressure exponents and
/// rounded, so its result departs from the exact derivative of the 75-term
/// polynomial by about 2e-12 relative at 10000 dbar, growing with pressure.
/// The Omega implementation differentiates the full-precision coefficients and
/// agrees with the exact derivative, evaluated in 60-digit arithmetic, to
/// around 1e-16. This tolerance therefore bounds GSW-C's rounding, not ours;
/// it is still far tighter than any real mistake would produce, and the
/// finite-difference check below pins the value independently.
const Real DerivDPRTol = 1e-10;

/// The temperature derivative passes through zero near the density maximum of
/// nearly fresh water, where a relative tolerance means nothing. This absolute
/// floor sits well above the roundoff of the polynomial sum that forms it
/// (terms of order 1e-5, so roundoff of order 1e-21) and far below any value
/// of physical interest (order 1e-7).
const Real DerivDCtATol = 1e-19;

/// Likewise for the thermal expansion coefficient, which is the temperature
/// derivative divided by the specific volume (order 1e-3).
const Real AlphaATol = 1e-16;

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

/// Test depth-mean specific volume calculation
void testDepthMeanSpecificVolume() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   Eos *TestEos      = Eos::getInstance();

   Array2DReal LayerThickness("LayerThickness", Mesh->NCellsSize,
                              VCoord->NVertLayers);

   deepCopy(TestEos->SpecVol, 2._Real);
   deepCopy(LayerThickness, 3._Real);
   deepCopy(TestEos->DepthMeanSpecificVolume, 0._Real);
   deepCopy(VCoord->SurfacePressure, 100._Real);
   deepCopy(VCoord->BottomGeomDepth, 500._Real);

   VCoord->computePressure(LayerThickness, VCoord->SurfacePressure);
   VCoord->computeGeomZHeight(LayerThickness, TestEos->SpecVol);

   TestEos->computeDepthMeanSpecificVolume(LayerThickness);

   auto DepthMeanSpecificVolumeH =
       createHostMirrorCopy(TestEos->DepthMeanSpecificVolume);

   int NumMismatches = 0;
   for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
      const Real Expected = 2._Real;
      if (!isApprox(DepthMeanSpecificVolumeH(ICell), Expected, RTol)) {
         LOG_ERROR("EosTest: DepthMeanSpecificVolume Bad Value: "
                   "DepthMeanSpecificVolume({}) = {}; Expected {}",
                   ICell, DepthMeanSpecificVolumeH(ICell), Expected);
         ++NumMismatches;
      }
   }

   if (NumMismatches != 0) {
      ABORT_ERROR("EosTest: DepthMeanSpecificVolume FAIL with {} bad values",
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
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 1 || K == 0) { // should be ref value
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
         // top layer should be ref value
         if (!isApprox(BruntVaisalaFreqSqH(I, 0), LinearBVFExpValue, RTol))
            LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      I, 0, BruntVaisalaFreqSqH(I, 0), LinearBVFExpValue);
         // K = 1 should be ref value
         if (!isApprox(BruntVaisalaFreqSqH(I, 1), LinearBVFExpValue, RTol))
            LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      I, 1, BruntVaisalaFreqSqH(I, 1), LinearBVFExpValue);
         // remaining values just check for other conditions
         for (int K = 2; K < NVertLayers + 1; ++K) {
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
          const int KMax   = MaxLayerCell(ICell) + 1;
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (K == 1 || K == 0) { // should be ref value
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
         // top layer should be ref value
         if (!isApprox(BruntVaisalaFreqSqH(ICell, 0), TeosBVFExpValue, RTol))
            LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      ICell, 0, BruntVaisalaFreqSqH(ICell, 0), TeosBVFExpValue);
         // K = 1 should be ref value
         if (!isApprox(BruntVaisalaFreqSqH(ICell, 1), TeosBVFExpValue, RTol))
            LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                      "BruntVaisala({},{}) = {}; Expected {}",
                      ICell, 1, BruntVaisalaFreqSqH(ICell, 1), TeosBVFExpValue);
         // remaining values just check for other conditions
         for (int K = 2; K < NVertLayers + 1; ++K) {
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

/// Test the array-level TEOS-10 specific volume derivatives over the mesh.
///
/// The state varies with depth rather than being uniform, so that the check
/// covers a range of values and exercises the vertical chunking, and the
/// expected values are obtained layer by layer from GSW-C on the host.
///
/// This is a test of the plumbing, not of the polynomial. It runs
/// Eos::computeSpecVolAndDerivs on the device over the whole mesh and so
/// covers the dispatch on EosChoice, the chunked vertical loop and its chunk
/// boundaries, the MinLayerCell/MaxLayerCell masking, the writing of the four
/// results into the Eos member arrays, and the registration of the derivative
/// fields in the Eos group. The math itself is covered point by point, over a
/// much wider range of states, by checkValueGswcSpecVolDerivs below; GSW-C
/// appears here only as a convenient source of expected values.
void testEosTeos10Derivs() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create the ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   deepCopy(TestEos->SpecVol, 0.0);
   deepCopy(TestEos->SpecVolDCt, 0.0);
   deepCopy(TestEos->SpecVolDSa, 0.0);
   deepCopy(TestEos->SpecVolDP, 0.0);

   /// A state that gets saltier, colder and deeper with depth, spanning a
   /// realistic part of the oceanographic range over the column
   parallelFor(
       "populateDerivArrays", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          SArray(ICell, K) = Sa + 0.1_Real * K;
          TArray(ICell, K) = Ct - 0.15_Real * K;
          PArray(ICell, K) = (100.0_Real + 150.0_Real * K) * Db2Pa;
       });

   /// Compute specific volume and its derivatives
   TestEos->computeSpecVolAndDerivs(TArray, SArray, PArray);

   /// Take local handles on the Eos members for the reduction kernels
   Array2DReal SpecVol    = TestEos->SpecVol;
   Array2DReal SpecVolDCt = TestEos->SpecVolDCt;
   Array2DReal SpecVolDSa = TestEos->SpecVolDSa;
   Array2DReal SpecVolDP  = TestEos->SpecVolDP;

   /// Expected values per layer from GSW-C, computed on the host and copied to
   /// the device for the comparison
   HostArray1DReal ExpSpecVolH("ExpSpecVolH", NVertLayers);
   HostArray1DReal ExpDCtH("ExpDCtH", NVertLayers);
   HostArray1DReal ExpDSaH("ExpDSaH", NVertLayers);
   HostArray1DReal ExpDPH("ExpDPH", NVertLayers);

   for (int K = 0; K < NVertLayers; ++K) {
      const double SaVal = Sa + 0.1 * K;
      const double CtVal = Ct - 0.15 * K;
      const double PDb   = 100.0 + 150.0 * K;

      double GswDSa, GswDCt, GswDP;
      gsw_specvol_first_derivatives(SaVal, CtVal, PDb, &GswDSa, &GswDCt,
                                    &GswDP);

      ExpSpecVolH(K) = gsw_specvol(SaVal, CtVal, PDb);
      ExpDCtH(K)     = GswDCt;
      ExpDSaH(K)     = GswDSa;
      ExpDPH(K)      = GswDP;
   }

   auto ExpSpecVol = createDeviceMirrorCopy(ExpSpecVolH);
   auto ExpDCt     = createDeviceMirrorCopy(ExpDCtH);
   auto ExpDSa     = createDeviceMirrorCopy(ExpDSaH);
   auto ExpDP      = createDeviceMirrorCopy(ExpDPH);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all active cells and layers against the expected values
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckSpecVolDerivs-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), ExpSpecVol(K), DerivRTol) or
                     !isApprox(SpecVolDCt(ICell, K), ExpDCt(K), DerivRTol,
                               DerivDCtATol) or
                     !isApprox(SpecVolDSa(ICell, K), ExpDSa(K), DerivRTol) or
                     !isApprox(SpecVolDP(ICell, K), ExpDP(K), DerivDPRTol)) {
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
      auto SpecVolH    = createHostMirrorCopy(SpecVol);
      auto SpecVolDCtH = createHostMirrorCopy(SpecVolDCt);
      auto SpecVolDSaH = createHostMirrorCopy(SpecVolDSa);
      auto SpecVolDPH  = createHostMirrorCopy(SpecVolDP);
      for (int I = 0; I < Mesh->NCellsAll; ++I) {
         for (int K = 0; K < NVertLayers; ++K) {
            if (!isApprox(SpecVolH(I, K), ExpSpecVolH(K), DerivRTol))
               LOG_ERROR("EosTest: SpecVol Deriv Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         I, K, SpecVolH(I, K), ExpSpecVolH(K));
            if (!isApprox(SpecVolDCtH(I, K), ExpDCtH(K), DerivRTol,
                          DerivDCtATol))
               LOG_ERROR("EosTest: SpecVolDCt Bad Value: "
                         "SpecVolDCt({},{}) = {}; Expected {}",
                         I, K, SpecVolDCtH(I, K), ExpDCtH(K));
            if (!isApprox(SpecVolDSaH(I, K), ExpDSaH(K), DerivRTol))
               LOG_ERROR("EosTest: SpecVolDSa Bad Value: "
                         "SpecVolDSa({},{}) = {}; Expected {}",
                         I, K, SpecVolDSaH(I, K), ExpDSaH(K));
            if (!isApprox(SpecVolDPH(I, K), ExpDPH(K), DerivDPRTol))
               LOG_ERROR("EosTest: SpecVolDP Bad Value: "
                         "SpecVolDP({},{}) = {}; Expected {}",
                         I, K, SpecVolDPH(I, K), ExpDPH(K));
         }
      }
      ABORT_ERROR("EosTest: SpecVol Derivs TEOS FAIL with {} bad values",
                  NumMismatches);
   }

   /// Check that each derivative is registered as a field in the Eos group
   /// with the member array attached, so that it can be written to a stream
   const std::string DerivFldNames[3] = {TestEos->SpecVolDCtFldName,
                                         TestEos->SpecVolDSaFldName,
                                         TestEos->SpecVolDPFldName};
   const Array2DReal DerivArrays[3]   = {SpecVolDCt, SpecVolDSa, SpecVolDP};

   for (int IFld = 0; IFld < 3; ++IFld) {
      const std::string &FldName = DerivFldNames[IFld];

      if (!Field::exists(FldName)) {
         ABORT_ERROR("EosTest: SpecVol Derivs field {} does not exist",
                     FldName);
      }

      if (!FieldGroup::isFieldInGroup(FldName, TestEos->EosGroupName)) {
         ABORT_ERROR("EosTest: SpecVol Derivs field {} not in group {}",
                     FldName, TestEos->EosGroupName);
      }

      auto DerivField = Field::get(FldName);
      auto FieldData  = DerivField->getDataArray<Array2DReal>();
      if (FieldData.data() != DerivArrays[IFld].data()) {
         ABORT_ERROR("EosTest: SpecVol Derivs field {} does not alias the "
                     "Eos member array",
                     FldName);
      }
   }

   return;
}

/// Test the array-level linear EOS specific volume derivatives, which are
/// known in closed form from the configured linear coefficients
void testEosLinearDerivs() {
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
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);
   deepCopy(TestEos->SpecVolDCt, 0.0);
   deepCopy(TestEos->SpecVolDSa, 0.0);
   deepCopy(TestEos->SpecVolDP, 0.0);

   TestEos->computeSpecVolAndDerivs(TArray, SArray, PArray);

   /// Take local handles on the Eos members for the reduction kernels
   Array2DReal SpecVol    = TestEos->SpecVol;
   Array2DReal SpecVolDCt = TestEos->SpecVolDCt;
   Array2DReal SpecVolDSa = TestEos->SpecVolDSa;
   Array2DReal SpecVolDP  = TestEos->SpecVolDP;

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against the expected values
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckSpecVolDerivs-linear", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), LinearExpValue, RTol) or
                     !isApprox(SpecVolDCt(ICell, K), LinearDCtExpValue, RTol) or
                     !isApprox(SpecVolDSa(ICell, K), LinearDSaExpValue, RTol) or
                     SpecVolDP(ICell, K) != 0.0_Real) {
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
      auto SpecVolDCtH = createHostMirrorCopy(SpecVolDCt);
      auto SpecVolDSaH = createHostMirrorCopy(SpecVolDSa);
      auto SpecVolDPH  = createHostMirrorCopy(SpecVolDP);
      for (int I = 0; I < Mesh->NCellsAll; ++I) {
         for (int K = 0; K < NVertLayers; ++K) {
            if (!isApprox(SpecVolDCtH(I, K), LinearDCtExpValue, RTol))
               LOG_ERROR("EosTest: SpecVolDCt Linear Bad Value: "
                         "SpecVolDCt({},{}) = {}; Expected {}",
                         I, K, SpecVolDCtH(I, K), LinearDCtExpValue);
            if (!isApprox(SpecVolDSaH(I, K), LinearDSaExpValue, RTol))
               LOG_ERROR("EosTest: SpecVolDSa Linear Bad Value: "
                         "SpecVolDSa({},{}) = {}; Expected {}",
                         I, K, SpecVolDSaH(I, K), LinearDSaExpValue);
            if (SpecVolDPH(I, K) != 0.0_Real)
               LOG_ERROR("EosTest: SpecVolDP Linear Bad Value: "
                         "SpecVolDP({},{}) = {}; Expected 0",
                         I, K, SpecVolDPH(I, K));
         }
      }
      ABORT_ERROR("EosTest: SpecVol Derivs Linear FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Test the array-level constant EOS specific volume derivatives, all of which
/// must be identically zero
void testEosConstantDerivs() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsSize       = Mesh->NCellsSize;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::ConstantEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsSize, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsSize, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsSize, NVertLayers);
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);
   deepCopy(TestEos->SpecVolDCt, 0.0);
   deepCopy(TestEos->SpecVolDSa, 0.0);
   deepCopy(TestEos->SpecVolDP, 0.0);

   TestEos->computeSpecVolAndDerivs(TArray, SArray, PArray);

   /// Take local handles on the Eos members for the reduction kernels
   Array2DReal SpecVol    = TestEos->SpecVol;
   Array2DReal SpecVolDCt = TestEos->SpecVolDCt;
   Array2DReal SpecVolDSa = TestEos->SpecVolDSa;
   Array2DReal SpecVolDP  = TestEos->SpecVolDP;

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against the expected values
   int NumMismatches = 0;
   parallelReduceOuter(
       "CheckSpecVolDerivs-Constant", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), ConstantExpValue, RTol) or
                     SpecVolDCt(ICell, K) != 0.0_Real or
                     SpecVolDSa(ICell, K) != 0.0_Real or
                     SpecVolDP(ICell, K) != 0.0_Real) {
                    InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   if (NumMismatches != 0) {
      ABORT_ERROR("EosTest: SpecVol Derivs Constant FAIL with {} bad values",
                  NumMismatches);
   }

   return;
}

/// Test all Eos::calcCtFreezing pathways (Teos10, Linear, Constant)
void testCalcCtFreezing() {
   const Real RTol = 1e-10;

   constexpr Real SaturationFrac = 0.0;
   constexpr Real PDb            = 500.0; // pressure in dbar for GSW pathway
   constexpr Real SaLocal        = 32.0;

   const Real CtTeosExpected =
       gsw_ct_freezing_poly(SaLocal, PDb, SaturationFrac);
   const Real CtTeos =
       Eos::calcCtFreezing(EosType::Teos10Eos, SaLocal, PDb, SaturationFrac);
   if (!isApprox(CtTeos, CtTeosExpected, RTol)) {
      ABORT_ERROR("testCalcCtFreezing: Teos10 FAIL, expected {}, got {}",
                  CtTeosExpected, CtTeos);
   }

   const Real CtLinearExpected = -0.054_Real * SaLocal / Psu2Gpkg;
   const Real CtLinear =
       Eos::calcCtFreezing(EosType::LinearEos, SaLocal, PDb, SaturationFrac);
   if (!isApprox(CtLinear, CtLinearExpected, RTol)) {
      ABORT_ERROR("testCalcCtFreezing: Linear FAIL, expected {}, got {}",
                  CtLinearExpected, CtLinear);
   }

   const Real CtConstExpected = -1.9_Real;
   const Real CtConst =
       Eos::calcCtFreezing(EosType::ConstantEos, SaLocal, PDb, SaturationFrac);
   if (!isApprox(CtConst, CtConstExpected, RTol)) {
      ABORT_ERROR("testCalcCtFreezing: Constant FAIL, expected {}, got {}",
                  CtConstExpected, CtConst);
   }

   return;
}

/// Check negative salinity handling in TEOS-10 EOS functions.
void testNegativeSalinityTeos10() {
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;
   Teos10Eos Teos(VertCoord::getDefault());

   constexpr Real PDb        = 500.0;
   constexpr Real PPa        = PDb * Db2Pa;
   constexpr Real Pt         = 10.0;
   const Real PtFromCtAtZero = TestEos->calcPtFromCt(0.0_Real, Ct);
   const Real CtFromPtAtZero = TestEos->calcCtFromPt(0.0_Real, Pt);
   const Real CtFreezingAtZero =
       Eos::calcCtFreezing(EosType::Teos10Eos, 0.0_Real, PDb, 0.0_Real);
   const Real CtFreezingLinearAtZero =
       Eos::calcCtFreezing(EosType::LinearEos, 0.0_Real, PDb, 0.0_Real);
   const Real NegativeSalinity[2] = {-1.0_Real, -2.0_Real};

   for (const Real SaValue : NegativeSalinity) {
      const Real PtFromCt = TestEos->calcPtFromCt(SaValue, Ct);
      const Real CtFromPt = TestEos->calcCtFromPt(SaValue, Pt);
      const Real CtFreezing =
          Eos::calcCtFreezing(EosType::Teos10Eos, SaValue, PDb, 0.0_Real);
      const Real CtFreezingLinear =
          Eos::calcCtFreezing(EosType::LinearEos, SaValue, PDb, 0.0_Real);
      if (!isApprox(PtFromCt, PtFromCtAtZero, RTol) ||
          !isApprox(CtFromPt, CtFromPtAtZero, RTol) ||
          !isApprox(CtFreezing, CtFreezingAtZero, RTol) ||
          !isApprox(CtFreezingLinear, CtFreezingLinearAtZero, RTol)) {
         ABORT_ERROR("testNegativeSalinityTeos10: clamped conversion mismatch "
                     "at Sa={}",
                     SaValue);
      }

      Real SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP;
      Teos.calcSpecVolAndDerivsAtPoint(Ct, SaValue, PPa, SpecVol, SpecVolDCt,
                                       SpecVolDSa, SpecVolDP);
      const Real Alpha = Eos::computeAlpha(EosType::Teos10Eos, SaValue, Ct, PDb,
                                           SpecVol, 0.0_Real);
      const Real Beta  = Eos::computeBeta(EosType::Teos10Eos, SaValue, Ct, PDb,
                                          SpecVol, 0.0_Real);
      if (Kokkos::isnan(SpecVol) || Kokkos::isinf(SpecVol) ||
          Kokkos::isnan(SpecVolDCt) || Kokkos::isinf(SpecVolDCt) ||
          Kokkos::isnan(SpecVolDSa) || Kokkos::isinf(SpecVolDSa) ||
          Kokkos::isnan(SpecVolDP) || Kokkos::isinf(SpecVolDP) ||
          Kokkos::isnan(Alpha) || Kokkos::isinf(Alpha) || Kokkos::isnan(Beta) ||
          Kokkos::isinf(Beta)) {
         ABORT_ERROR("testNegativeSalinityTeos10: non-finite EOS output at "
                     "Sa={}",
                     SaValue);
      }
   }
   LOG_INFO("testNegativeSalinityTeos10: PASS");
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

/// Test that the calcCtFreezing function returns the expected value
void checkValueGswcCtFreezing() {
   const Real RTol = 1e-10;

   constexpr Real SaturationFrac = 0.0;
   constexpr Real P              = 500.0 * Db2Pa; // Convert dbar to Pa
   constexpr Real Sa             = 32.0;

   /// Get freezing temperature from GSW-C library
   double CtFreezGswc = gsw_ct_freezing_poly(Sa, P * Pa2Db, SaturationFrac);
   double CtFreez =
       Teos10Eos::calcCtFreezingTeos10(Sa, P * Pa2Db, SaturationFrac);

   /// Check the value against the GSW-C value
   bool Check = isApprox(CtFreezGswc, CtFreez, RTol);
   if (!Check) {
      ABORT_ERROR("checkValueGswcCtFreezing: CtFreez FAIL, expected {}, got {}",
                  CtFreezGswc, CtFreez);
   }
   return;
}

/// Relative difference between two values, zero when both vanish
Real relDiff(Real X, Real Y) {
   const Real Scale = std::max(std::abs(X), std::abs(Y));
   return Scale > 0.0 ? std::abs(X - Y) / Scale : 0.0;
}

/// Test the TEOS-10 specific volume derivatives against the GSW-C library over
/// a range of states.
///
/// GSW-C is used here unmodified and through its public API, as an independent
/// oracle. The Omega implementation does not derive from it: the derivatives
/// are the analytic derivatives of the Roquet et al. 2015 polynomial that the
/// Teos10Eos functor already carries.
///
/// This is the test of the polynomial itself, as opposed to
/// testEosTeos10Derivs above, which tests the array-level machinery. It calls
/// the point-wise calcSpecVolAndDerivsAtPoint on the host at every combination
/// of the salinity, temperature and pressure values in SaTest, CtTest and
/// PTest, which reach the corners of the oceanographic range -- fresh and
/// salty, freezing and warm, surface and 10000 dbar -- rather than the single
/// realistic profile the mesh test uses. A dropped or mis-scaled term shows up
/// here, and the tolerances are tight enough to say so.
void checkValueGswcSpecVolDerivs() {

   Teos10Eos TestEos(VertCoord::getDefault());

   int NumBad     = 0;
   Real WorstSv   = 0.0;
   Real WorstDCt  = 0.0;
   Real WorstDSa  = 0.0;
   Real WorstDP   = 0.0;
   int NumChecked = 0;

   for (int ISa = 0; ISa < NSaTest; ++ISa) {
      for (int ICt = 0; ICt < NCtTest; ++ICt) {
         for (int IP = 0; IP < NPTest; ++IP) {

            const Real SaVal = SaTest[ISa];
            const Real CtVal = CtTest[ICt];
            const Real PDb   = PTest[IP];

            Real SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP;
            TestEos.calcSpecVolAndDerivsAtPoint(CtVal, SaVal, PDb * Db2Pa,
                                                SpecVol, SpecVolDCt, SpecVolDSa,
                                                SpecVolDP);

            /// GSW-C takes pressure in dbar and returns the derivatives per
            /// (g/kg), per degC, and per Pa
            double GswDSa, GswDCt, GswDP;
            gsw_specvol_first_derivatives(SaVal, CtVal, PDb, &GswDSa, &GswDCt,
                                          &GswDP);
            const double GswSpecVol = gsw_specvol(SaVal, CtVal, PDb);

            WorstSv  = std::max(WorstSv, relDiff(SpecVol, GswSpecVol));
            WorstDCt = std::max(WorstDCt, relDiff(SpecVolDCt, GswDCt));
            WorstDSa = std::max(WorstDSa, relDiff(SpecVolDSa, GswDSa));
            WorstDP  = std::max(WorstDP, relDiff(SpecVolDP, GswDP));
            ++NumChecked;

            if (!isApprox(SpecVol, GswSpecVol, DerivRTol)) {
               LOG_ERROR("EosTest: SpecVol Bad Value at Sa={}, Ct={}, "
                         "P={} dbar: expected {}, got {}",
                         SaVal, CtVal, PDb, GswSpecVol, SpecVol);
               ++NumBad;
            }
            if (!isApprox(SpecVolDCt, GswDCt, DerivRTol, DerivDCtATol)) {
               LOG_ERROR("EosTest: SpecVolDCt Bad Value at Sa={}, Ct={}, "
                         "P={} dbar: expected {}, got {}",
                         SaVal, CtVal, PDb, GswDCt, SpecVolDCt);
               ++NumBad;
            }
            if (!isApprox(SpecVolDSa, GswDSa, DerivRTol)) {
               LOG_ERROR("EosTest: SpecVolDSa Bad Value at Sa={}, Ct={}, "
                         "P={} dbar: expected {}, got {}",
                         SaVal, CtVal, PDb, GswDSa, SpecVolDSa);
               ++NumBad;
            }
            if (!isApprox(SpecVolDP, GswDP, DerivDPRTol)) {
               LOG_ERROR("EosTest: SpecVolDP Bad Value at Sa={}, Ct={}, "
                         "P={} dbar: expected {}, got {}",
                         SaVal, CtVal, PDb, GswDP, SpecVolDP);
               ++NumBad;
            }
         }
      }
   }

   LOG_INFO("EosTest: TEOS-10 derivatives vs GSW-C over {} states, max "
            "relative difference: SpecVol {}, d/dCt {}, d/dSa {}, d/dP {}",
            NumChecked, WorstSv, WorstDCt, WorstDSa, WorstDP);

   if (NumBad != 0) {
      ABORT_ERROR("EosTest: SpecVol derivatives vs GSW-C FAIL with {} bad "
                  "values",
                  NumBad);
   }

   return;
}

/// Check the specific volume derivatives against centered finite differences of
/// the Omega specific volume itself.
///
/// This overlaps the GSW-C comparison above whenever that library is present
/// and correct, and that is deliberate: it pins the unit convention of the
/// Omega interface -- per degC, per (g/kg), and per Pa -- without reference to
/// GSW. A change that silently made the pressure derivative per dbar instead,
/// a factor of 1e4, would be caught here as well as there.
void checkFiniteDiffSpecVolDerivs() {

   Teos10Eos TestEos(VertCoord::getDefault());

   /// Step sizes and tolerance are set by the finite difference itself: large
   /// enough that the difference of two specific volumes is not lost to
   /// roundoff, small enough that the truncation error stays below the
   /// tolerance.
   const Real DCtStep  = 1.0e-2;  // degC
   const Real DSaStep  = 1.0e-2;  // g/kg
   const Real DPStep   = 1.0e5;   // Pa (10 dbar)
   const Real FDRTol   = 1.0e-5;  // limited by finite difference truncation
   const Real FDCtATol = 1.0e-14; // finite difference noise floor
   const Real FDSaATol = 1.0e-14;
   const Real FDPATol  = 1.0e-20;

   /// Evaluate only the specific volume at a perturbed state
   auto SpecVolAt = [&TestEos](Real CtVal, Real SaVal, Real PPa) {
      Real SpecVol, DCt, DSa, DP;
      TestEos.calcSpecVolAndDerivsAtPoint(CtVal, SaVal, PPa, SpecVol, DCt, DSa,
                                          DP);
      return SpecVol;
   };

   int NumBad = 0;

   for (int ISa = 0; ISa < NSaTest; ++ISa) {
      for (int ICt = 0; ICt < NCtTest; ++ICt) {
         for (int IP = 0; IP < NPTest; ++IP) {

            const Real SaVal = SaTest[ISa];
            const Real CtVal = CtTest[ICt];
            const Real PPa   = PTest[IP] * Db2Pa;

            Real SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP;
            TestEos.calcSpecVolAndDerivsAtPoint(
                CtVal, SaVal, PPa, SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP);

            const Real FDCt = (SpecVolAt(CtVal + DCtStep, SaVal, PPa) -
                               SpecVolAt(CtVal - DCtStep, SaVal, PPa)) /
                              (2.0_Real * DCtStep);
            const Real FDSa = (SpecVolAt(CtVal, SaVal + DSaStep, PPa) -
                               SpecVolAt(CtVal, SaVal - DSaStep, PPa)) /
                              (2.0_Real * DSaStep);
            const Real FDP = (SpecVolAt(CtVal, SaVal, PPa + DPStep) -
                              SpecVolAt(CtVal, SaVal, PPa - DPStep)) /
                             (2.0_Real * DPStep);

            if (!isApprox(SpecVolDCt, FDCt, FDRTol, FDCtATol)) {
               LOG_ERROR("EosTest: SpecVolDCt disagrees with finite difference "
                         "at Sa={}, Ct={}, P={} dbar: {} vs {}",
                         SaVal, CtVal, PTest[IP], SpecVolDCt, FDCt);
               ++NumBad;
            }
            if (!isApprox(SpecVolDSa, FDSa, FDRTol, FDSaATol)) {
               LOG_ERROR("EosTest: SpecVolDSa disagrees with finite difference "
                         "at Sa={}, Ct={}, P={} dbar: {} vs {}",
                         SaVal, CtVal, PTest[IP], SpecVolDSa, FDSa);
               ++NumBad;
            }
            if (!isApprox(SpecVolDP, FDP, FDRTol, FDPATol)) {
               LOG_ERROR("EosTest: SpecVolDP disagrees with finite difference "
                         "at Sa={}, Ct={}, P={} dbar: {} vs {}",
                         SaVal, CtVal, PTest[IP], SpecVolDP, FDP);
               ++NumBad;
            }
         }
      }
   }

   if (NumBad != 0) {
      ABORT_ERROR("EosTest: SpecVol derivatives vs finite differences FAIL "
                  "with {} bad values",
                  NumBad);
   }

   return;
}

/// Test the thermal expansion and haline contraction coefficients used by the
/// TEOS-10 Brunt-Vaisala frequency against the GSW-C library.
///
/// These are the specific volume derivatives divided by the specific volume,
/// so this covers the same polynomial from the other side. Before this check
/// existed, calcAlpha and calcBeta were exercised only through the single
/// hardcoded BruntVaisalaFreqSq value below, which is too loose to catch a
/// mistake in either of them.
void checkValueGswcAlphaBeta() {

   Teos10Eos TestEos(VertCoord::getDefault());
   Teos10BruntVaisalaFreqSq TestBvf(VertCoord::getDefault());

   int NumBad      = 0;
   Real WorstAlpha = 0.0;
   Real WorstBeta  = 0.0;

   for (int ISa = 0; ISa < NSaTest; ++ISa) {
      for (int ICt = 0; ICt < NCtTest; ++ICt) {
         for (int IP = 0; IP < NPTest; ++IP) {

            const Real SaVal = SaTest[ISa];
            const Real CtVal = CtTest[ICt];
            const Real PDb   = PTest[IP];

            Real SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP;
            TestEos.calcSpecVolAndDerivsAtPoint(CtVal, SaVal, PDb * Db2Pa,
                                                SpecVol, SpecVolDCt, SpecVolDSa,
                                                SpecVolDP);

            const Real Alpha = TestBvf.calcAlpha(SaVal, CtVal, PDb, SpecVol);
            const Real Beta  = TestBvf.calcBeta(SaVal, CtVal, PDb, SpecVol);

            double GswSpecVol, GswAlpha, GswBeta;
            gsw_specvol_alpha_beta(SaVal, CtVal, PDb, &GswSpecVol, &GswAlpha,
                                   &GswBeta);

            WorstAlpha = std::max(WorstAlpha, relDiff(Alpha, GswAlpha));
            WorstBeta  = std::max(WorstBeta, relDiff(Beta, GswBeta));

            if (!isApprox(Alpha, GswAlpha, DerivRTol, AlphaATol)) {
               LOG_ERROR("EosTest: Alpha Bad Value at Sa={}, Ct={}, "
                         "P={} dbar: expected {}, got {}",
                         SaVal, CtVal, PDb, GswAlpha, Alpha);
               ++NumBad;
            }
            if (!isApprox(Beta, GswBeta, DerivRTol)) {
               LOG_ERROR("EosTest: Beta Bad Value at Sa={}, Ct={}, "
                         "P={} dbar: expected {}, got {}",
                         SaVal, CtVal, PDb, GswBeta, Beta);
               ++NumBad;
            }
         }
      }
   }

   LOG_INFO("EosTest: alpha and beta vs GSW-C, max relative difference: "
            "alpha {}, beta {}",
            WorstAlpha, WorstBeta);

   if (NumBad != 0) {
      ABORT_ERROR("EosTest: alpha and beta vs GSW-C FAIL with {} bad values",
                  NumBad);
   }

   return;
}

/// Test that the Eos CT-from-PT helper matches GSW-C
void checkValueGswcCtFromPt() {
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   Real CtExpValue = gsw_ct_from_pt(Sa, Ct);
   Real CtTeos     = TestEos->calcCtFromPt(Sa, Ct);
   bool Check      = isApprox(CtTeos, CtExpValue, RTol);
   if (!Check) {
      ABORT_ERROR("checkValueGswcCtFromPt: Ct FAIL, expected {}, got {}",
                  CtExpValue, CtTeos);
   }
   return;
}

/// Test that the Eos PT-from-CT helper matches GSW-C
void checkValueGswcPtFromCt() {
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   Real PtExpValue = gsw_pt_from_ct(Sa, Ct);
   Real PtTeos     = TestEos->calcPtFromCt(Sa, Ct);
   bool Check      = isApprox(PtTeos, PtExpValue, RTol);
   if (!Check) {
      ABORT_ERROR("checkValueGswcPtFromCt: Pt FAIL, expected {}, got {}",
                  PtExpValue, PtTeos);
   }
   return;
}

// the main tests (all in one to have the same log):
// Single value test:
// --> test calls the external GSW-C library
// and compares the specific volume to the published value
// --> next compares the TEOS-10 specific volume and its three first
// derivatives against GSW-C over a range of states
// --> next checks those derivatives against centered finite differences of the
// Omega specific volume, which pins the unit conventions independently of GSW
// --> next compares the thermal expansion and haline contraction coefficients
// against GSW-C over the same range of states
// Full array tests:
// --> one tests the value on a Eos with linear option
// --> next checks the value on a Eos with linear displaced option
// --> next checks the value of the linear squared Brunt Vaisala Freq.
// calculation
// --> next checks the value on a Eos with TEOS-10 option
// --> next checks the value on a Eos with TEOS-10 displaced option
// --> next checks the value of the TOES-10 squared Brunt Vaisala Freq.
// calculation
// --> last checks the specific volume derivatives for each of the three EOS
// options over the whole mesh
void eosTest(const std::string &MeshFile = "OmegaMesh.nc") {
   initEosTest(MeshFile);
   const auto &Mesh = HorzMesh::getDefault();

   checkValueGswcSpecVol();
   checkValueGswcN2();
   checkValueGswcCtFreezing();
   checkValueGswcCtFromPt();
   checkValueGswcPtFromCt();
   checkValueGswcSpecVolDerivs();
   checkFiniteDiffSpecVolDerivs();
   checkValueGswcAlphaBeta();

   testEosLinear();
   testEosLinearDisplaced();
   testEosConstant();
   testDepthMeanSpecificVolume();
   testBruntVaisalaFreqSqLinear();
   testEosTeos10();
   testEosTeos10Displaced();
   testBruntVaisalaFreqSqTeos10();
   testCalcCtFreezing();
   testNegativeSalinityTeos10();
   testEosTeos10Derivs();
   testEosLinearDerivs();
   testEosConstantDerivs();

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
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   // If we made it here, test is successful
   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
