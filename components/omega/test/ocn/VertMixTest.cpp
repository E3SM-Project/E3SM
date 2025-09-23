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

using namespace OMEGA;

/// Test constants and expected values
constexpr int NVertLevels      = 60;

/// Published values (TEOS-10 and linear) to test against
const Real VertDiffExpValueN    = 1.1650682538996415;      // Expected value for diffusivity for positive BVF
const Real VertViscExpValueN    = 1.0515534893999001; // Expected value for viscosity for positive BVF
const Real VertDiffExpValueP    = 0.0010490674327851117; // Expected value for diffusivity for negative BVF
const Real VertViscExpValueP    = 0.0018542271307222758; // Expected value for viscosity for negative BVF
const Real VertDiffBackExp      = 1.0e-5;                   // Expected value for background diffusivity
const Real VertViscBackExp      = 1.0e-4;                   // Expected value for background viscosity
const Real VertDiffConvExp      = 1.0;                    // Expected value for convective diffusivity/viscosity
const Real VertDiffShearExp     = 0.1650582538996415;   // Expected value for shear diffusivity
const Real VertViscShearExp     = 0.0514534893999002;  // Expected value for shear viscosity

/// Test input values
double BVFP                 = 0.1;     // Positive Brunt-Vaisala frequency in s^-2
double BVFN                 = -0.1;    // Negative Brunt-Vaisala frequency in s^-2
double NV                   = 1.0;     // Normal velocity in m/s
double TV                   = 1.0;     // Tangential velocity in m/s
const Real RTol             = 1e-10;   // Relative tolerance for isApprox checks

/// The initialization routine for Eos testing. It calls various
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

   /// Initialize parallel IO
   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("VertMixTest: error initializing parallel IO");
   }

   /// Initialize decomposition
   Decomp::init(mesh);

   /// Initialize mesh
   HorzMesh::init();

   /// Create vertical dimension for test arrays
   const auto &Mesh = HorzMesh::getDefault();
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);

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
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   VertMix *TestVertMix = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NVArray  = Array2DReal("NVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal TVArray  = Array2DReal("TVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal BVFArray = Array2DReal("BVFArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal ZArray   = Array2DReal("ZArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NVArray, NV);
   deepCopy(TVArray, TV);
   deepCopy(BVFArray, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
          "populateArrays", {Mesh->NCellsAll, NVertLevels},
          KOKKOS_LAMBDA(I4 ICell, I4 K) {
             ZArray(ICell, K) = -K;
             NVArray(ICell, K) = NV + 0.5*K;
             TVArray(ICell, K) = TV + 0.5*K;
          });

   /// Compute only background vertical viscosity and diffusivity
   TestVertMix->BackDiff = 1.0e-5;
   TestVertMix->BackVisc = 1.0e-4;
   TestVertMix->ComputeVertMixConv.Enabled = false;
   TestVertMix->ComputeVertMixShear.Enabled = false;
   TestVertMix->computeVertMix(NVArray, TVArray, BVFArray, ZArray);
   Array2DReal BackVertVisc = TestVertMix->VertVisc;
   Array2DReal BackVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int numMismatches = 0;
   if (!isApprox(BackVertVisc(1, 1), VertViscBackExp, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixBack: VertVisc isApprox FAIL, "
                "expected {}, got {} with {} mismatches", 
                VertViscBackExp, BackVertVisc(1, 1), 
                numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixBack-Visc: PASS");
   }

   /// Check total Diff against linear addition of components
   numMismatches = 0;
   if (!isApprox(BackVertDiff(1, 1), VertDiffBackExp, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixBack: VertDiff isApprox FAIL, "
                "expected {}, got {} with {} mismatches", 
                VertDiffBackExp, BackVertDiff(1, 1), 
                numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixBack-Diff: PASS");
   }

   return Err;
}

int testPPConvVertMix() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   VertMix *TestVertMix = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NVArray  = Array2DReal("NVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal TVArray  = Array2DReal("TVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal BVFArray = Array2DReal("BVFArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal ZArray   = Array2DReal("ZArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NVArray, NV);
   deepCopy(TVArray, TV);
   deepCopy(BVFArray, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
          "populateArrays", {Mesh->NCellsAll, NVertLevels},
          KOKKOS_LAMBDA(I4 ICell, I4 K) {
             ZArray(ICell, K) = -K;
             NVArray(ICell, K) = NV + 0.5*K;
             TVArray(ICell, K) = TV + 0.5*K;
          });

   /// Compute only convective vertical viscosity and diffusivity
   TestVertMix->BackDiff = 0.0;
   TestVertMix->BackVisc = 0.0;
   TestVertMix->ComputeVertMixConv.Enabled = true;
   TestVertMix->ComputeVertMixShear.Enabled = false;
   TestVertMix->computeVertMix(NVArray, TVArray, BVFArray, ZArray);
   Array2DReal ConvVertVisc = TestVertMix->VertVisc;
   Array2DReal ConvVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int numMismatches = 0;
   if (!isApprox(ConvVertVisc(1, 1), VertDiffConvExp, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixConv: VertVisc isApprox FAIL, "
                "expected {}, got {} with {} mismatches", 
                VertDiffConvExp, ConvVertVisc(1, 1), 
                numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixConv-Visc: PASS");
   }

   /// Check total Diff against linear addition of components
   numMismatches = 0;
   if (!isApprox(ConvVertDiff(1, 1), VertDiffConvExp, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixConv: VertDiff isApprox FAIL, "
                "expected {}, got {} with {} mismatches", 
                VertDiffConvExp, ConvVertDiff(1, 1), 
                numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixConv-Diff: PASS");
   }

   return Err;
}

int testPPShearVertMix() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   VertMix *TestVertMix = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NVArray  = Array2DReal("NVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal TVArray  = Array2DReal("TVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal BVFArray = Array2DReal("BVFArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal ZArray   = Array2DReal("ZArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NVArray, NV);
   deepCopy(TVArray, TV);
   deepCopy(BVFArray, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
          "populateArrays", {Mesh->NCellsAll, NVertLevels},
          KOKKOS_LAMBDA(I4 ICell, I4 K) {
             ZArray(ICell, K) = -K;
             NVArray(ICell, K) = NV + 0.5*K;
             TVArray(ICell, K) = TV + 0.5*K;
          });

   /// Compute only shear vertical viscosity and diffusivity
   TestVertMix->BackDiff = 0.0;
   TestVertMix->BackVisc = 0.0;
   TestVertMix->ComputeVertMixConv.Enabled = false;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(NVArray, TVArray, BVFArray, ZArray);
   Array2DReal ShearVertVisc = TestVertMix->VertVisc;
   Array2DReal ShearVertDiff = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int numMismatches = 0;
   if (!isApprox(ShearVertVisc(1, 1), VertViscShearExp, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixShear: VertVisc isApprox FAIL, "
                "expected {}, got {} with {} mismatches", 
                VertViscShearExp, ShearVertVisc(1, 1), 
                numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixShear-Visc: PASS");
   }

   /// Check total Diff against linear addition of components
   numMismatches = 0;
   if (!isApprox(ShearVertDiff(1, 1), VertDiffShearExp, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixShear: VertDiff isApprox FAIL, "
                "expected {}, got {} with {} mismatches", 
                VertDiffShearExp, ShearVertDiff(1, 1), 
                numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixShear-Diff: PASS");
   }

   return Err;
}

/// Test vertical mixing coefficients calculation for all cells/levels
int testPPTotalVertMix() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   VertMix *TestVertMix = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal NVArray  = Array2DReal("NVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal TVArray  = Array2DReal("TVArray", Mesh->NEdgesAll, NVertLevels);
   Array2DReal BVFArray = Array2DReal("BVFArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal ZArray   = Array2DReal("ZArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(NVArray, NV);
   deepCopy(TVArray, TV);

   // Test with positive BVF first
   deepCopy(BVFArray, BVFP);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   parallelFor(
          "populateArrays", {Mesh->NCellsAll, NVertLevels},
          KOKKOS_LAMBDA(I4 ICell, I4 K) {
             ZArray(ICell, K) = -K;
             NVArray(ICell, K) = NV + 0.5*K;
             TVArray(ICell, K) = TV + 0.5*K;
          });

   /// Compute vertical viscosity and diffusivity
   TestVertMix->BackDiff = 1.0e-5;
   TestVertMix->BackVisc = 1.0e-4;
   TestVertMix->ComputeVertMixConv.Enabled = true;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(NVArray, TVArray, BVFArray, ZArray);

   /// Check all VertDiff array values against expected value
   int numMismatches   = 0;
   Array2DReal VertDiffP = TestVertMix->VertDiff;
   if (!isApprox(VertDiffP(1, 1), VertDiffExpValueP, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixTotal: VertDiffPositive isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertDiffExpValueP, VertDiffP(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMixTotal VertDiffPositive: PASS");
   }

   /// Check all VertVisc array values against expected value
   numMismatches = 0;
   Array2DReal VertViscP = TestVertMix->VertVisc;
   if (!isApprox(VertViscP(1, 1), VertViscExpValueP, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixTotal: VertViscPositive isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertViscExpValueP, VertViscP(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMixTotal VertViscPositive: PASS");
   }

   // Now test with negative BVF
   deepCopy(BVFArray, BVFN);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   /// Compute vertical viscosity and diffusivity
   TestVertMix->computeVertMix(NVArray, TVArray, BVFArray, ZArray);

   /// Check all VertDiff array values against expected value
   numMismatches = 0;
   Array2DReal VertDiffN = TestVertMix->VertDiff;
   if (!isApprox(VertDiffN(1, 1), VertDiffExpValueN, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMix: VertDiffNegative isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertDiffExpValueN, VertDiffN(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMix VertDiffNegative: PASS");
   }

   /// Check all VertVisc array values against expected value
   numMismatches = 0;
   Array2DReal VertViscN = TestVertMix->VertVisc;
   if (!isApprox(VertViscN(1, 1), VertViscExpValueN, RTol)) {
         numMismatches = 1;
   } else {
         numMismatches = 0;
   }
   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMix: VertViscNegative isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                VertViscExpValueN, VertViscN(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMix VertViscNegative : PASS");
   }

   return Err;
}

/// Finalize and clean up all test infrastructure
void finalizeVertMixTest() {
   HorzMesh::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

// the main test (all in one to have the same log)
// --> one tests the vertical diffusivity and viscosity
// with background, convective, and shear on
// --> one tests the linear superposition of the 
// background, convective, and shear contributions
// --> one tests the Brunt-Vaisala frequency calculation
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

// The test driver for Eos testing
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
