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
constexpr int NVertLevels    = 60;

/// Published values (TEOS-10 and linear) to test against
const Real VertDiffExpValue    = 0.0009732819628;      // Expected value for diffusivity
const Real VertViscExpValue    = 0.0009784735812133072; // Expected value for viscosity

/// Test input values
double BVF                  = 30.0;       // Absolute Salinity in g/kg
double U                    = 10.0;       // Conservative Temperature in degC
double V                    = 1000.0;     // Pressure in dbar
double SVol                 = 0.0;        // Specific volume
const Real RTol             = 1e-10;      // Relative tolerance for isApprox checks

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

/// Test vertical mixing coefficients calculation for all cells/levels
int testVertMix() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   VertMix *TestVertMix = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal UArray   = Array2DReal("UArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal VArray   = Array2DReal("VArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal BVFArray = Array2DReal("BVFArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(UArray, U);
   deepCopy(VArray, V);
   deepCopy(BVFArray, BVF);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   /// Compute vertical viscosity and diffusivity
   TestVertMix->computeVertMix(UArray, VArray, BVFArray);

   /// Check all VertDiff array values against expected value
   int numMismatches = 0;
   Array2DReal VertDiff = TestVertMix->VertDiff;
   parallelReduce("CheckVertDiffMatrix", {Mesh->NCellsAll, NVertLevels},
                  KOKKOS_LAMBDA(int i, int j, int &localCount) {
                     if (!isApprox(VertDiff(i, j), VertDiffExpValue, RTol)) {
                        localCount++;
                     }
                  },
                  numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMix: VertDiff isApprox FAIL, "
                "expected {}, got {} mismatches",
                VertDiffExpValue, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMix VertDiff: PASS");
   }

   /// Check all VertVisc array values against expected value
   numMismatches = 0;
   Array2DReal VertVisc = TestVertMix->VertVisc;
   parallelReduce("CheckVertViscMatrix", {Mesh->NCellsAll, NVertLevels},
                  KOKKOS_LAMBDA(int i, int j, int &localCount) {
                     if (!isApprox(VertVisc(i, j), VertViscExpValue, RTol)) {
                        localCount++;
                     }
                  },
                  numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMix: VertVisc isApprox FAIL, "
                "expected {}, got {} mismatches",
                VertViscExpValue, numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("TestVertMix VertVisc: PASS");
   }

   return Err;
}

/// Test linear superposition of vertical mixing coefficients
int testVertMixLinearity() {
   int Err          = 0;
   const auto *Mesh = HorzMesh::getDefault();
   /// Get Eos instance to test
   VertMix *TestVertMix = VertMix::getInstance();
   TestVertMix->VertMixChoice = VertMixType::PP;

   /// Create and fill ocean state arrays
   Array2DReal UArray   = Array2DReal("UArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal VArray   = Array2DReal("VArray", Mesh->NCellsAll, NVertLevels);
   Array2DReal BVFArray = Array2DReal("BVFArray", Mesh->NCellsAll, NVertLevels);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(UArray, U);
   deepCopy(VArray, V);
   deepCopy(BVFArray, BVF);
   deepCopy(TestVertMix->VertDiff, 0.0);
   deepCopy(TestVertMix->VertVisc, 0.0);

   /// Compute background vertical viscosity and diffusivity
   TestVertMix->BackDiff = 1.0e-5;
   TestVertMix->BackVisc = 1.0e-4;
   TestVertMix->ComputeVertMixConv.Enabled = false;
   TestVertMix->ComputeVertMixShear.Enabled = false;
   TestVertMix->computeVertMix(UArray, VArray, BVFArray);
   Array2DReal BackVertVisc = TestVertMix->VertVisc;
   Array2DReal BackVertDiff = TestVertMix->VertDiff;

   /// Compute convective vertical viscosity and diffusivity
   TestVertMix->BackDiff = 0.0;
   TestVertMix->BackVisc = 0.0;
   TestVertMix->ComputeVertMixConv.Enabled = true;
   TestVertMix->ComputeVertMixShear.Enabled = false;
   TestVertMix->computeVertMix(UArray, VArray, BVFArray);
   Array2DReal ConvVertVisc = TestVertMix->VertVisc;
   Array2DReal ConvVertDiff = TestVertMix->VertDiff;

   /// Compute shear vertical viscosity and diffusivity
   TestVertMix->BackDiff = 0.0;
   TestVertMix->BackVisc = 0.0;
   TestVertMix->ComputeVertMixConv.Enabled = false;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(UArray, VArray, BVFArray);
   Array2DReal ShearVertVisc = TestVertMix->VertVisc;
   Array2DReal ShearVertDiff = TestVertMix->VertDiff;

   /// Compute vertical viscosity and diffusivity
   TestVertMix->BackDiff = 1.0e-5;
   TestVertMix->BackVisc = 1.0e-4;
   TestVertMix->ComputeVertMixConv.Enabled = true;
   TestVertMix->ComputeVertMixShear.Enabled = true;
   TestVertMix->computeVertMix(UArray, VArray, BVFArray);
   Array2DReal VertViscTotal = TestVertMix->VertVisc;
   Array2DReal VertDiffTotal = TestVertMix->VertDiff;

   /// Check total Visc against linear addition of components
   int numMismatches = 0;
   parallelReduce("CheckVertViscMatrix-Linear", {Mesh->NCellsAll, NVertLevels},
                  KOKKOS_LAMBDA(int i, int j, int &localCount) {
                     if (!isApprox(VertViscTotal(i, j), BackVertVisc(i, j) + 
                        ConvVertVisc(i, j) + ShearVertVisc(i, j), RTol)) {
                        localCount++;
                     }
                  },
                  numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixLinearity: Linear addition of visc components isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                BackVertVisc(1, 1) + ConvVertVisc(1, 1) + ShearVertVisc(1, 1), 
                VertViscTotal(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixLinearity-Visc: PASS");
   }

   /// Check total Diff against linear addition of components
   numMismatches = 0;
   parallelReduce("CheckVertDiffMatrix-Linear", {Mesh->NCellsAll, NVertLevels},
                  KOKKOS_LAMBDA(int i, int j, int &localCount) {
                     if (!isApprox(VertDiffTotal(i, j), BackVertDiff(i, j) + 
                        ConvVertDiff(i, j) + ShearVertDiff(i, j), RTol)) {
                        localCount++;
                     }
                  },
                  numMismatches);

   if (numMismatches != 0) {
      Err++;
      LOG_ERROR("TestVertMixLinearity: Linear addition of diff components isApprox FAIL, "
                "expected {}, got {} with {} mismatches",
                BackVertDiff(1, 1) + ConvVertDiff(1, 1) + ShearVertDiff(1, 1), 
                VertDiffTotal(1, 1), numMismatches);
   }
   if (Err == 0) {
      LOG_INFO("VertMixTest TestVertMixLinearity-Diff: PASS");
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

   Err += testVertMix();
   Err += testVertMixLinearity();

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
