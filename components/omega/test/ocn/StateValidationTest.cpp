//===-- Test driver for OMEGA StateValidation -----------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for ocean state validation
///
/// Tests the validateOceanState function by verifying that it passes on valid
/// state data. Checks cover LayerThickness, KineticEnergyCell, Temperature,
/// and Salinity fields.
//
//===-----------------------------------------------------------------------===/

#include "StateValidation.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "mpi.h"

#include <cmath>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Initialize the Omega subsystems required for state validation testing

int initStateValidationTest(const std::string &MeshFile) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ StateValidation unit tests ------");

   Config("Omega");
   Config::readAll("omega.yml");

   TimeStepper::init1();

   IO::init(DefComm);
   Decomp::init(MeshFile);

   IOStream::init();

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("StateValidationTest: error initializing default halo");
   }

   HorzMesh::init();
   VertCoord::init();
   Tracers::init();

   int StateErr = OceanState::init();
   if (StateErr != 0) {
      Err++;
      LOG_ERROR("StateValidationTest: error initializing default state");
   }

   VertAdv::init();

   return Err;
}

//------------------------------------------------------------------------------
// Fill state and auxiliary/tracer arrays with known-valid values

static int fillValidState() {
   int Err = 0;

   auto *State      = OceanState::getDefault();
   const int NCells = State->NCellsAll;
   const int NVert  = State->NVertLayers;
   const int NEdges = State->NEdgesAll;

   // LayerThickness: fill with 100 m (valid range [1e-10, 1000])
   Array2DReal LayerThick = State->getLayerThickness(0);
   parallelFor(
       "FillLayerThick", {NCells, NVert},
       KOKKOS_LAMBDA(int ICell, int K) { LayerThick(ICell, K) = 100.0; });

   // NormalVelocity: fill with 0 (not checked, but needed for AuxState)
   Array2DReal NormalVel = State->getNormalVelocity(0);
   parallelFor(
       "FillNormalVel", {NEdges, NVert},
       KOKKOS_LAMBDA(int IEdge, int K) { NormalVel(IEdge, K) = 0.0; });

   // Exchange halos so auxiliary state computations are consistent
   State->exchangeHalo(0);

   // Tracers: fill Temperature with 10 C and Salinity with 35 g/kg
   // Use deepCopy with individual tracer subviews
   if (Tracers::getNumTracers() > 0) {
      // Temperature = 10.0 (valid: -10 to 50)
      if (Tracers::IndxTemp != Tracers::IndxInvalid) {
         Array2DReal TempArr = Tracers::getByIndex(0, Tracers::IndxTemp);
         deepCopy(TempArr, static_cast<Real>(10.0));
      }

      // Salinity = 35.0 (valid: -2 to 60)
      if (Tracers::IndxSalt != Tracers::IndxInvalid) {
         Array2DReal SaltArr = Tracers::getByIndex(0, Tracers::IndxSalt);
         deepCopy(SaltArr, static_cast<Real>(35.0));
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Run state validation tests

int testStateValidation() {
   int Err = 0;

   // Initialize the auxiliary state (needed for KineticEnergyCell)
   AuxiliaryState::init();
   auto *DefAuxState = AuxiliaryState::getDefault();

   if (!DefAuxState) {
      Err++;
      LOG_ERROR("StateValidationTest: Default AuxiliaryState not found");
      return Err;
   }

   auto *DefState = OceanState::getDefault();
   if (!DefState) {
      Err++;
      LOG_ERROR("StateValidationTest: Default OceanState not found");
      return Err;
   }

   // Fill state arrays with valid values
   Err += fillValidState();

   // Compute auxiliary variables so KineticEnergyCell is populated
   {
      Array3DReal AllTracers = Tracers::getAll(0);
      DefAuxState->computeAll(DefState, AllTracers, 0);
   }

   // Test: validation should pass on valid state (no abort expected)
   LOG_INFO("StateValidationTest: Testing validation on valid state");
   validateOceanState(DefState, DefAuxState, 0);
   LOG_INFO("StateValidationTest: Valid state validation PASS");

   AuxiliaryState::clear();
   return Err;
}

//------------------------------------------------------------------------------
// Finalize Omega objects

void finalizeStateValidationTest() {
   Tracers::clear();
   OceanState::clear();
   VertAdv::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

//------------------------------------------------------------------------------
// Main entry point

int main(int argc, char *argv[]) {
   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      int Err = initStateValidationTest("OmegaMesh.nc");
      if (Err != 0) {
         LOG_CRITICAL("StateValidationTest: Error during initialization");
      } else {
         RetVal += testStateValidation();
      }
      finalizeStateValidationTest();
   }

   if (RetVal == 0)
      LOG_INFO("------ StateValidation unit tests successful ------");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
