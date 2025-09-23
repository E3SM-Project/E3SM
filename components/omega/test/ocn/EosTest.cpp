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

#include <gswteos-10.h>

using namespace OMEGA;

/// Test constants and expected values
constexpr int NVertLayers = 60;

/// Published values (TEOS-10 and linear) to test against
const Real TeosSVExpValue = 0.0009732819628; // Expected value for TEOS-10 specific volume
const Real LinearExpValue =
    0.0009784735812133072; // Expected value for Linear specific volume
const Real TeosBVFExpValue = 
    0.020855053182884893; // Expected value for TEOS-10 Brunt-Vaisala frequency
const Real SimpleBVFExpValue =
    0.017833905406889013; // Expected value for Linear Brunt-Vaisala frequency
const Real GswBVFExpValue = 0.02081197958166906; // Expected value from GSW-C library

/// Test input values
double Sa       = 30.0;   // Absolute Salinity in g/kg
double Ct       = 10.0;   // Conservative Temperature in degC
double P        = 1000.0; // Pressure in dbar
const I4 KDisp  = 1;      // Displace parcel to K=1 for TEOS-10 eos
const Real RTol = 1e-10;  // Relative tolerance for isApprox checks

/// The initialization routine for Eos testing. It calls various
/// init routines, including the creation of the default decomposition.
I4 initEosTest(const std::string &mesh) {

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

   /// Initialize parallel IO
   IO::init(DefComm);

   /// Initialize decomposition
   Decomp::init(mesh);

   /// Initialize Halo
   Halo::init();

   /// Initialize mesh
   HorzMesh::init();

   /// Initialize vertical coordinate
   VertCoord::init(false);

   /// Initialize Eos
   Eos::init();

   /// Retrieve Eos
   Eos *DefEos = Eos::getInstance();
   if (DefEos) {
      LOG_INFO("EosTest: Eos retrieval PASS");
   } else {
      Err++;
      LOG_INFO("EosTest: Eos retrieval FAIL");
      return -1;
   }

   return Err;
}

/// Test Linear EOS calculation for all cells/layers
int testEosLinear() {
   int Err            = 0;
   const auto *Mesh   = HorzMesh::getDefault();
   const auto *VCoord = VertCoord::getDefault();
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", Mesh->NCellsAll, NVertLayers);
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

   auto SpecVolH = createHostMirrorCopy(SpecVol);
   if (NumMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: SpecVol Linear isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                LinearExpValue, SpecVolH(1, 1), NumMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalc Linear: PASS");
   }

   return Err;
}

/// Test Linear EOS calculation with vertical displacement
int testEosLinearDisplaced() {
   int Err            = 0;
   const auto *Mesh   = HorzMesh::getDefault();
   const auto *VCoord = VertCoord::getDefault();
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", Mesh->NCellsAll, NVertLayers);
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

   auto SpecVolDisplacedH = createHostMirrorCopy(SpecVolDisplaced);
   if (NumMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: Linear SpecVolDisp isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                LinearExpValue, SpecVolDisplacedH(1, 1), NumMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalcDisp Linear: PASS");
   }

   return Err;
}

/// Test simple Brunt-Vaisala frequency calculation for all cells/layers
int testBruntVaisalaFreqSimple() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal TArray = Array2DReal("TArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal PArray = Array2DReal("PArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal ZArray = Array2DReal("ZArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(TestEos->SpecVol, 0.0);
   deepCopy(TestEos->BruntVaisalaFreq, 0.0);

   //parallelFor(
   //       "populateArrays", {Mesh->NCellsAll, NVertLevels},
   //       KOKKOS_LAMBDA(I4 ICell, I4 K) {
   //          ZArray(ICell,K) = -K;
   //          SArray(ICell,K) = Sa + 0.1*K;
   //          TArray(ICell,K) = Ct + 0.25*K;
   //          PArray(ICell,K) = P + 0.25*K;
   //       });

   parallelFor(
          "populateArrays", {Mesh->NCellsAll, NVertLevels},
          KOKKOS_LAMBDA(I4 ICell, I4 K) {
             ZArray(ICell,0) = -992.1173890198451; ZArray(ICell,1) = -993.1071379053125; ZArray(ICell,2) = -994.0968821072275;
             SArray(ICell,0) = Sa-1.0; SArray(ICell,1) = Sa; SArray(ICell,2) = Sa+1.0;
             TArray(ICell,0) = Ct+15.0; TArray(ICell,1) = Ct+10.0; TArray(ICell,2) = Ct+5.0;
             PArray(ICell,0) = P; PArray(ICell,1) = P+1.0; PArray(ICell,2) = P+2.0;
          });

   /// Compute specific volume first
   TestEos->computeSpecVol(TArray, SArray, PArray);
   Array2DReal SpecVol = TestEos->SpecVol;

   /// Compute Brunt-Vaisala frequency
   TestEos->computeBruntVaisalaFreq(TArray, SArray, PArray, SpecVol, ZArray);
   Array2DReal BruntVaisalaFreq = TestEos->BruntVaisalaFreq;

   /// Check all array values against expected value
   int numMismatches   = 0;
   if (!isApprox(BruntVaisalaFreq(1, 1), SimpleBVFExpValue, RTol)) {
         numMismatches = 1;
      } else {
         numMismatches = 0;
      }
   //int numMismatches   = 0;
   //parallelReduce(
   //    "CheckSpecVolMatrix-Teos", {Mesh->NCellsAll, NVertLevels},
   //    KOKKOS_LAMBDA(int i, int j, int &localCount) {
   //       if (!isApprox(BruntVaisalaFreq(i, j), SimpleBVFExpValue, RTol)) {
   //          localCount++;
   //       }
   //    },
   //    numMismatches);

   auto BruntVaisalaFreqH = createHostMirrorCopy(BruntVaisalaFreq);
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: Simple BruntVaisalaFreq isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                SimpleBVFExpValue, BruntVaisalaFreqH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest BruntVaisalaFreqCalc Simple: PASS");
   }

   return Err;
}

/// Test TEOS-10 EOS calculation for all cells/layers
int testEosTeos10() {
   int Err            = 0;
   const auto *Mesh   = HorzMesh::getDefault();
   const auto *VCoord = VertCoord::getDefault();
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", Mesh->NCellsAll, NVertLayers);
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

   auto SpecVolH = createHostMirrorCopy(SpecVol);
   if (NumMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: TEOS SpecVol isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                TeosSVExpValue, SpecVolH(1, 1), NumMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalc TEOS-10: PASS");
   }

   return Err;
}

/// Test TEOS-10 EOS calculation with vertical displacement
int testEosTeos10Displaced() {
   int Err            = 0;
   const auto *Mesh   = HorzMesh::getDefault();
   const auto *VCoord = VertCoord::getDefault();
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", Mesh->NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", Mesh->NCellsAll, NVertLayers);
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

   auto SpecVolDisplacedH = createHostMirrorCopy(SpecVolDisplaced);
   if (NumMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: TEOS SpecVolDisp isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                TeosSVExpValue, SpecVolDisplacedH(1, 1), NumMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest SpecVolCalcDisp TEOS-10: PASS");
   }

   return Err;
}

/// Test TEOS-10 Brunt-Vaisala frequency calculation for all cells/levels
int testBruntVaisalaFreqTeos10() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray  = Array2DReal("SArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal TArray  = Array2DReal("TArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal PArray  = Array2DReal("PArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal ZArray  = Array2DReal("ZArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal SpArray = Array2DReal("SpArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SpArray, 0.0);
   deepCopy(TestEos->BruntVaisalaFreq, 0.0);
   deepCopy(TestEos->SpecVol, 0.0);

   //parallelFor(
   //       "populateArrays", {Mesh->NCellsAll, NVertLevels},
   //       KOKKOS_LAMBDA(I4 ICell, I4 K) {
   //          ZArray(ICell,K) = -K;
   //          SArray(ICell,K) = Sa + 0.1*K;
   //          TArray(ICell,K) = Ct + 0.25*K;
   //          PArray(ICell,K) = P + 0.25*K;
   //       });

   parallelFor(
          "populateArrays", {Mesh->NCellsAll, NVertLevels},
          KOKKOS_LAMBDA(I4 ICell, I4 K) {
             ZArray(ICell,0) = -992.1173890198451; ZArray(ICell,1) = -993.1071379053125; ZArray(ICell,2) = -994.0968821072275;
             SArray(ICell,0) = Sa-1.0; SArray(ICell,1) = Sa; SArray(ICell,2) = Sa+1.0;
             TArray(ICell,0) = Ct+15.0; TArray(ICell,1) = Ct+10.0; TArray(ICell,2) = Ct+5.0;
             PArray(ICell,0) = P; PArray(ICell,1) = P+1.0; PArray(ICell,2) = P+2.0;
          });

   /// Compute specific volume first
   TestEos->computeSpecVol(TArray, SArray, PArray);
   Array2DReal SpecVol = TestEos->SpecVol;

   /// Compute Brunt-Vaisala frequency
   TestEos->computeBruntVaisalaFreq(TArray, SArray, PArray, SpecVol, ZArray);
   Array2DReal BruntVaisalaFreq = TestEos->BruntVaisalaFreq; 

   /// Check all array values against expected value
   int numMismatches   = 0;
   if (!isApprox(BruntVaisalaFreq(1, 1), TeosBVFExpValue, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   //parallelReduce(
   //    "CheckSpecVolMatrix-Teos", {Mesh->NCellsAll, NVertLevels},
   //    KOKKOS_LAMBDA(int i, int j, int &localCount) {
   //       if (!isApprox(BruntVaisalaFreq(i, j), TeosBVFExpValue, RTol)) {
   //          localCount++;
   //       }
   //    },
   //    numMismatches);

   auto BruntVaisalaFreqH = createHostMirrorCopy(BruntVaisalaFreq);
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("EosTest: TEOS BruntVaisalaFreq isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                TeosBVFExpValue, BruntVaisalaFreqH(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("EosTest BruntVaisalaFreqCalc TEOS-10: PASS");
   }

   return Err;
}

/// Finalize and clean up all test infrastructure
void finalizeEosTest() {
   VertCoord::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

/// Test that the external GSW-C library returns the expected specific volume
int checkValueGswcSpecVol() {
   int Err         = 0;
   const Real RTol = 1e-10;

   /// Get specific volume from GSW-C library
   double SpecVol = gsw_specvol(Sa, Ct, P);
   /// Check the value against the expected TEOS-10 value
   bool Check = isApprox(SpecVol, TeosSVExpValue, RTol);
   if (!Check) {
      Err++;
      LOG_ERROR(
          "checkValueGswcSpecVol: SpecVol isApprox FAIL, expected {}, got {}",
          TeosSVExpValue, SpecVol);
   }
   if (Err == 0) {
      LOG_INFO("checkValueGswcSpecVol: PASS");
   }
   return Err;
}

/// Test that the external GSW-C library returns the expected N2
int checkValueGswcN2() {
   int Err         = 0;
   const Real RTol = 1e-10;

   // Number of intervals (nz)
   int nz = 2;

   // Input arrays: length nz+1
   double SA[4] = {Sa-1.0, Sa, Sa+1.0}; // Absolute Salinity (g/kg)
   double CT[4] = {Ct+15.0, Ct+10.0, Ct+5.0}; // Conservative Temperature (deg C)
   double p[4]  = {P, P+1.0, P+2.0};    // Pressure (dbar)

   // Latitude (degrees north)
   double latitude[4] = {0.0, 0.0, 0.0};

   // Output arrays: length nz
   double N2[nz];     // Brunt–Väisälä frequency squared
   double p_mid[nz];  // Midpoint pressure

   /// Get specific volume from GSW-C library
   gsw_nsquared(SA, CT, p, latitude, nz, N2, p_mid);

   /// Check the value against the expected TEOS-10 value
   bool Check = isApprox(N2[0], GswBVFExpValue, RTol);
   if (!Check) {
      Err++;
      LOG_ERROR(
          "checkValueGswcN2: N2 isApprox FAIL, expected {}, got {}",
          TeosBVFExpValue, N2[0]);
   }
   if (Err == 0) {
      LOG_INFO("checkValueGswcN2: PASS");
   }
   return Err;
}

// the main test (all in one to have the same log)
// Single value test:
// --> test calls the external GSW-C library
// and compares the specific volume to the published value
// Full array tests:
// --> one tests the value on a Eos with linear option
// --> next checks the value on a Eos with linear displaced option
// --> next checks the value on a Eos with TEOS-10 option
// --> next checks the value on a Eos with TEOS-10 displaced option
int eosTest(const std::string &MeshFile = "OmegaMesh.nc") {
   int Err = initEosTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("EosTest: Error initializing");
   }
   const auto &Mesh = HorzMesh::getDefault();

   LOG_INFO("Single value checks:");
   Err += checkValueGswcSpecVol();
   Err += checkValueGswcN2();

   LOG_INFO("Full array checks:");
   Err += testEosLinear();
   Err += testEosLinearDisplaced();
   Err += testBruntVaisalaFreqSimple();
   Err += testEosTeos10();
   Err += testEosTeos10Displaced();
   Err += testBruntVaisalaFreqTeos10();

   if (Err == 0) {
      LOG_INFO("EosTest: Successful completion");
   }
   finalizeEosTest();

   return Err;
}

// The test driver for Eos testing
int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += eosTest();

   Eos::destroyInstance();
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
