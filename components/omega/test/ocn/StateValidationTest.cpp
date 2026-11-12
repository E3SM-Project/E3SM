//===-- Test driver for OMEGA StateValidation -----------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for ocean state validation
///
/// Tests both the positive path (valid state passes) and the negative paths
/// (NaN and out-of-bounds values in each checked field are detected) of
/// validateOceanState / checkOceanState. Checked fields are:
///   - PseudoThickness, KineticEnergyCell, Temperature, Salinity.
//
//===-----------------------------------------------------------------------===/

#include "StateValidation.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
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
#include <limits>

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
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   Decomp::init(MeshFile);

   IOStream::init(ModelClock);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("StateValidationTest: error initializing default halo");
   }

   HorzMesh::init(ModelClock);
   VertCoord::init();
   Tracers::init();
   Eos::init();

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

   // PseudoThickness: fill with 100 m (valid range [1e-10, 1000])
   Array2DReal PseudoThick = State->getPseudoThickness(0);
   parallelFor(
       "FillPseudoThick", {NCells, NVert},
       KOKKOS_LAMBDA(int ICell, int K) { PseudoThick(ICell, K) = 100.0; });

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
      if (Tracers::IndxTemp != -1) {
         Array2DReal TempArr = Tracers::getByIndex(0, Tracers::IndxTemp);
         deepCopy(TempArr, static_cast<Real>(10.0));
      }

      // Salinity = 35.0 (valid: -2 to 60)
      if (Tracers::IndxSalt != -1) {
         Array2DReal SaltArr = Tracers::getByIndex(0, Tracers::IndxSalt);
         deepCopy(SaltArr, static_cast<Real>(35.0));
      }
   }

   return Err;
}

//------------------------------------------------------------------------------
// Restore the default state to valid values and recompute the auxiliary state

static void restoreValidState(OceanState *State, AuxiliaryState *AuxState) {
   fillValidState();
   Array3DReal AllTracers = Tracers::getAll(0);
   TimeInterval Interval(1., TimeUnits::Seconds);
   AuxState->computeAll(State, AllTracers, 0, Interval);
}

//------------------------------------------------------------------------------
// Positive test: validate a clean, valid state — expects 0 errors

static int testValidState(OceanState *State, AuxiliaryState *AuxState,
                          VertCoord *VCoord) {
   LOG_INFO("StateValidationTest: Testing validation on valid state");
   validateOceanState(State, AuxState, VCoord, 0);
   LOG_INFO("StateValidationTest: Valid state validation PASS");
   return 0;
}

//------------------------------------------------------------------------------
// Negative tests: inject an invalid value, verify checkOceanState returns > 0,
// then restore valid state.
//
// Returns 0 on success (i.e. the error was caught), 1 otherwise.

static int expectErrors(const char *TestName, OceanState *State,
                        AuxiliaryState *AuxState, VertCoord *VCoord) {
   auto [NaNs, OOBs] = checkOceanState(State, AuxState, VCoord, 0);
   auto Errs         = NaNs + OOBs;
   if (Errs == 0) {
      LOG_ERROR("StateValidationTest: {} - expected errors but got none",
                TestName);
      return 1;
   }
   LOG_INFO("StateValidationTest: {} PASS (caught {} error(s))", TestName,
            Errs);
   return 0;
}

// --- PseudoThickness ---

static int testNaNPseudoThickness(OceanState *State, AuxiliaryState *AuxState,
                                  VertCoord *VCoord) {
   restoreValidState(State, AuxState);
   const Real NaN   = std::numeric_limits<Real>::quiet_NaN();
   Array2DReal PT   = State->getPseudoThickness(0);
   const int NCells = State->NCellsAll;
   const int NVert  = State->NVertLayers;
   parallelFor(
       "InjectNaNPseudoThick", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { PT(I, K) = NaN; });
   return expectErrors("NaN in PseudoThickness", State, AuxState, VCoord);
}

static int testOOBHighPseudoThickness(OceanState *State,
                                      AuxiliaryState *AuxState,
                                      VertCoord *VCoord) {
   restoreValidState(State, AuxState);
   Array2DReal PT   = State->getPseudoThickness(0);
   const int NCells = State->NCellsAll;
   const int NVert  = State->NVertLayers;
   // 2000 m is above the valid max of 1000 m
   parallelFor(
       "InjectOOBHighPseudoThick", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { PT(I, K) = 2000.0; });
   return expectErrors("OOB-high in PseudoThickness", State, AuxState, VCoord);
}

static int testOOBLowPseudoThickness(OceanState *State,
                                     AuxiliaryState *AuxState,
                                     VertCoord *VCoord) {
   restoreValidState(State, AuxState);
   Array2DReal PT   = State->getPseudoThickness(0);
   const int NCells = State->NCellsAll;
   const int NVert  = State->NVertLayers;
   // -1.0 m is below the valid min of 1e-10 m
   parallelFor(
       "InjectOOBLowPseudoThick", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { PT(I, K) = -1.0; });
   return expectErrors("OOB-low in PseudoThickness", State, AuxState, VCoord);
}

// --- KineticEnergyCell ---

static int testNaNKineticEnergy(OceanState *State, AuxiliaryState *AuxState,
                                VertCoord *VCoord) {
   restoreValidState(State, AuxState);
   const Real NaN   = std::numeric_limits<Real>::quiet_NaN();
   Array2DReal KE   = AuxState->KineticAux.KineticEnergyCell;
   const int NCells = State->NCellsAll;
   const int NVert  = State->NVertLayers;
   parallelFor(
       "InjectNaNKE", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { KE(I, K) = NaN; });
   return expectErrors("NaN in KineticEnergyCell", State, AuxState, VCoord);
}

static int testOOBKineticEnergy(OceanState *State, AuxiliaryState *AuxState,
                                VertCoord *VCoord) {
   restoreValidState(State, AuxState);
   Array2DReal KE   = AuxState->KineticAux.KineticEnergyCell;
   const int NCells = State->NCellsAll;
   const int NVert  = State->NVertLayers;
   // 9999 J/kg is above the valid max of 10 J/kg
   parallelFor(
       "InjectOOBKE", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { KE(I, K) = 9999.0; });
   return expectErrors("OOB in KineticEnergyCell", State, AuxState, VCoord);
}

// --- Temperature tracer ---

static int testNaNTemperature(OceanState *State, AuxiliaryState *AuxState,
                              VertCoord *VCoord) {
   if (Tracers::IndxTemp == -1)
      return 0; // tracer not active; skip
   restoreValidState(State, AuxState);
   const Real NaN      = std::numeric_limits<Real>::quiet_NaN();
   Array2DReal TempArr = Tracers::getByIndex(0, Tracers::IndxTemp);
   const int NCells    = State->NCellsAll;
   const int NVert     = State->NVertLayers;
   parallelFor(
       "InjectNaNTemp", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { TempArr(I, K) = NaN; });
   return expectErrors("NaN in Temperature", State, AuxState, VCoord);
}

static int testOOBTemperature(OceanState *State, AuxiliaryState *AuxState,
                              VertCoord *VCoord) {
   if (Tracers::IndxTemp == -1)
      return 0; // tracer not active; skip
   restoreValidState(State, AuxState);
   Array2DReal TempArr = Tracers::getByIndex(0, Tracers::IndxTemp);
   const int NCells    = State->NCellsAll;
   const int NVert     = State->NVertLayers;
   // 9999 C is above the valid max of 50 C
   parallelFor(
       "InjectOOBTemp", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { TempArr(I, K) = 9999.0; });
   return expectErrors("OOB in Temperature", State, AuxState, VCoord);
}

// --- Salinity tracer ---

static int testNaNSalinity(OceanState *State, AuxiliaryState *AuxState,
                           VertCoord *VCoord) {
   if (Tracers::IndxSalt == -1)
      return 0; // tracer not active; skip
   restoreValidState(State, AuxState);
   const Real NaN      = std::numeric_limits<Real>::quiet_NaN();
   Array2DReal SaltArr = Tracers::getByIndex(0, Tracers::IndxSalt);
   const int NCells    = State->NCellsAll;
   const int NVert     = State->NVertLayers;
   parallelFor(
       "InjectNaNSalt", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { SaltArr(I, K) = NaN; });
   return expectErrors("NaN in Salinity", State, AuxState, VCoord);
}

static int testOOBSalinity(OceanState *State, AuxiliaryState *AuxState,
                           VertCoord *VCoord) {
   if (Tracers::IndxSalt == -1)
      return 0; // tracer not active; skip
   restoreValidState(State, AuxState);
   Array2DReal SaltArr = Tracers::getByIndex(0, Tracers::IndxSalt);
   const int NCells    = State->NCellsAll;
   const int NVert     = State->NVertLayers;
   // 9999 g/kg is above the valid max of 60 g/kg
   parallelFor(
       "InjectOOBSalt", {NCells, NVert},
       KOKKOS_LAMBDA(int I, int K) { SaltArr(I, K) = 9999.0; });
   return expectErrors("OOB in Salinity", State, AuxState, VCoord);
}

//------------------------------------------------------------------------------
// Run all state validation tests

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

   VertCoord *VCoord = VertCoord::getDefault();

   // Fill state arrays with valid values and compute auxiliary state
   fillValidState();
   {
      Array3DReal AllTracers = Tracers::getAll(0);
      TimeInterval Interval(1., TimeUnits::Seconds);
      DefAuxState->computeAll(DefState, AllTracers, 0, Interval);
   }

   // ---- Positive test ----
   Err += testValidState(DefState, DefAuxState, VCoord);

   // ---- Negative tests: each injects a bad value and verifies detection ----

   // PseudoThickness
   Err += testNaNPseudoThickness(DefState, DefAuxState, VCoord);
   Err += testOOBHighPseudoThickness(DefState, DefAuxState, VCoord);
   Err += testOOBLowPseudoThickness(DefState, DefAuxState, VCoord);

   // KineticEnergyCell
   Err += testNaNKineticEnergy(DefState, DefAuxState, VCoord);
   Err += testOOBKineticEnergy(DefState, DefAuxState, VCoord);

   // Temperature tracer
   Err += testNaNTemperature(DefState, DefAuxState, VCoord);
   Err += testOOBTemperature(DefState, DefAuxState, VCoord);

   // Salinity tracer
   Err += testNaNSalinity(DefState, DefAuxState, VCoord);
   Err += testOOBSalinity(DefState, DefAuxState, VCoord);

   AuxiliaryState::clear();
   return Err;
}

//------------------------------------------------------------------------------
// Finalize Omega objects

void finalizeStateValidationTest() {
   Tracers::clear();
   Eos::destroyInstance();
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
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
