//===-- Test driver for OMEGA State -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA state class
///
/// This driver tests that the OMEGA state class member variables are read in
/// correctly from a sample shperical mesh file. Also tests that the time level
/// update works as expected.
//
//===-----------------------------------------------------------------------===/

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
#include "VertCoord.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The initialization routine for State testing. It calls various
// init routines, including the creation of the default decomposition.

void initStateTest() {

   int Err = 0;
   Error Err1;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize Logging
   initLogging(DefEnv);

   // Read Omega config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize the default time stepper and retrieve the model clock
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Initialize IOStreams - this does not yet validate the contents
   // of each file, only creates streams from Config
   IOStream::init(ModelClock);

   // Initialize Field infrastructure
   Field::init(ModelClock);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0)
      ABORT_ERROR("State: error initializing default halo");

   // Initialize the default mesh
   HorzMesh::init();

   // Initialize the vertical coordinate
   VertCoord::init();

   // Initialize tracers
   Tracers::init();

   // Initialize Aux State variables
   AuxiliaryState::init();

   // Create tendencies
   Tendencies::init();

   // Finish time stepper initialization
   TimeStepper::init2();

   // Create a default ocean state
   Err = OceanState::init();
   if (Err != 0)
      ABORT_ERROR("ocnInit: Error initializing default state");

   // Now that all fields have been defined, validate all the streams
   // contents
   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid)
      ABORT_ERROR("ocnInit: Error validating IO Streams");

   // Read the state variables from the initial state stream
   Metadata ReqMeta; // no global metadata needed for init state read
   Err1 = IOStream::read("InitialState", ModelClock, ReqMeta);
   CHECK_ERROR_ABORT(Err1, "ocnInit: Error reading initial state from stream");

   // Finish initialization of initial state by filling halos and copying
   // to host. Current time level is zero.
   OceanState *DefState = OceanState::getDefault();
   DefState->exchangeHalo(0);
   DefState->copyToHost(0);

   return;
}

// Check for differences between layer thickness and normal velocity host arrays
int checkHost(OceanState *RefState, OceanState *TstState, int RefTimeLevel,
              int TstTimeLevel) {

   int count = 0;
   HostArray2DReal LayerThicknessHRef =
       RefState->getLayerThicknessH(RefTimeLevel);
   HostArray2DReal LayerThicknessHTst =
       TstState->getLayerThicknessH(TstTimeLevel);
   for (int Cell = 0; Cell < RefState->NCellsAll; Cell++) {
      for (int Layer = 0; Layer < RefState->NVertLayers; Layer++) {
         if (LayerThicknessHRef(Cell, Layer) !=
             LayerThicknessHTst(Cell, Layer)) {
            count++;
         }
      }
   }

   HostArray2DReal NormalVelocityHRef =
       RefState->getNormalVelocityH(RefTimeLevel);
   HostArray2DReal NormalVelocityHTst =
       TstState->getNormalVelocityH(TstTimeLevel);
   for (int Edge = 0; Edge < RefState->NEdgesAll; Edge++) {
      for (int Layer = 0; Layer < RefState->NVertLayers; Layer++) {
         if (NormalVelocityHRef(Edge, Layer) !=
             NormalVelocityHTst(Edge, Layer)) {
            count++;
         }
      }
   }

   return count;
}

// Check for differences between layer thickness and normal velocity device
// arrays
int checkDevice(OceanState *RefState, OceanState *TstState, int RefTimeLevel,
                int TstTimeLevel) {

   int count1;
   Array2DReal LayerThicknessRef = RefState->getLayerThickness(RefTimeLevel);
   Array2DReal LayerThicknessTst = TstState->getLayerThickness(TstTimeLevel);
   parallelReduce(
       "reduce", {RefState->NCellsAll, RefState->NVertLayers},
       KOKKOS_LAMBDA(int Cell, int Layer, int &Accum) {
          if (LayerThicknessRef(Cell, Layer) !=
              LayerThicknessTst(Cell, Layer)) {
             Accum++;
          }
       },
       count1);

   int count2;
   Array2DReal NormalVelocityRef = RefState->getNormalVelocity(RefTimeLevel);
   Array2DReal NormalVelocityTst = TstState->getNormalVelocity(TstTimeLevel);
   parallelReduce(
       "reduce", {RefState->NCellsAll, RefState->NVertLayers},
       KOKKOS_LAMBDA(int Cell, int Layer, int &Accum) {
          if (NormalVelocityRef(Cell, Layer) !=
              NormalVelocityTst(Cell, Layer)) {
             Accum++;
          }
       },
       count2);

   return count1 + count2;
}

//------------------------------------------------------------------------------
// The test driver for State -> This tests the time level update of state
// variables and verifies the state is read in correctly.
//
int main(int argc, char *argv[]) {

   int RetVal = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {

      // Call initialization routine to create default state and other
      // quantities
      initStateTest();

      // Get default mesh, halo and time data
      HorzMesh *DefHorzMesh   = HorzMesh::getDefault();
      Halo *DefHalo           = Halo::getDefault();
      TimeStepper *DefStepper = TimeStepper::getDefault();
      Clock *ModelClock       = DefStepper->getClock();
      int NTimeLevelsDef      = DefStepper->getNTimeLevels();

      // Test retrieval of the default state
      OceanState *DefState = OceanState::getDefault();
      int NVertLayers      = DefState->NVertLayers;
      int NCellsAll        = DefState->NCellsAll;
      int NEdgesAll        = DefState->NEdgesAll;
      int CurTime          = 0;
      int NewTime          = 1;
      if (DefState and NVertLayers > 0) {
         LOG_INFO("State: Default state retrieval PASS");
      } else {
         RetVal += 1;
         LOG_INFO("State: Default state retrieval FAIL");
      }

      // Get default state arrays and check for reasonable values
      HostArray2DReal LayerThickHDef = DefState->getLayerThicknessH(CurTime);
      HostArray2DReal NormalVelocityHDef =
          DefState->getNormalVelocityH(CurTime);
      Array2DReal LayerThickDef     = DefState->getLayerThickness(CurTime);
      Array2DReal NormalVelocityDef = DefState->getNormalVelocity(CurTime);

      int Count1 = 0;
      int Count2 = 0;
      for (int Cell = 0; Cell < NCellsAll; Cell++) {
         for (int Layer = 0; Layer < NVertLayers; Layer++) {
            R8 val = LayerThickHDef(Cell, Layer);
            if (val != 0.0)
               ++Count1;                     // check for all-zero array
            if (val < 0.0 and val > 300.0) { // out of range
               Count2++;
            }
         }
      }
      if (Count2 == 0 and Count1 > 0) {
         LOG_INFO("State: State read PASS");
      } else {
         RetVal += 1;
         LOG_INFO("State: State read FAIL");
         LOG_INFO("State: Out-of-range {} Non-zero: {}", Count2, Count1);
      }

      // Test time swapping with 2 and higher numbers of time levels
      for (int NTimeLevels = 2; NTimeLevels < 5; NTimeLevels++) {

         // Create new reference and test states with the number of time levels
         OceanState *RefState = OceanState::create(
             "Reference", DefHorzMesh, DefHalo, NVertLayers, NTimeLevels);
         OceanState *TstState = OceanState::create("Test", DefHorzMesh, DefHalo,
                                                   NVertLayers, NTimeLevels);

         // Test successful creation
         if (TstState and RefState) { // true if non-null ptr
            LOG_INFO("State: Test state creation (NTimeLevels={}) PASS",
                     NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO("State: Test state creation (NTimeLevels={}) FAIL",
                     NTimeLevels);
         }

         // Create test arrays and fill reference and test states at current
         // time with Default state
         HostArray2DReal LayerThickHRef = RefState->getLayerThicknessH(CurTime);
         HostArray2DReal LayerThickHTst = TstState->getLayerThicknessH(CurTime);
         HostArray2DReal NormalVelocityHRef =
             RefState->getNormalVelocityH(CurTime);
         HostArray2DReal NormalVelocityHTst =
             TstState->getNormalVelocityH(CurTime);

         for (int Cell = 0; Cell < NCellsAll; Cell++) {
            for (int Layer = 0; Layer < NVertLayers; Layer++) {
               LayerThickHRef(Cell, Layer) = LayerThickHDef(Cell, Layer);
               LayerThickHTst(Cell, Layer) = LayerThickHDef(Cell, Layer);
            }
         }
         for (int Edge = 0; Edge < NEdgesAll; Edge++) {
            for (int Layer = 0; Layer < NVertLayers; Layer++) {
               NormalVelocityHRef(Edge, Layer) =
                   NormalVelocityHDef(Edge, Layer);
               NormalVelocityHTst(Edge, Layer) =
                   NormalVelocityHDef(Edge, Layer);
            }
         }
         RefState->copyToDevice(CurTime);
         TstState->copyToDevice(CurTime);

         // Fill reference and test states at next time levels with the
         // Default state + 1
         LayerThickHRef     = RefState->getLayerThicknessH(NewTime);
         LayerThickHTst     = TstState->getLayerThicknessH(NewTime);
         NormalVelocityHRef = RefState->getNormalVelocityH(NewTime);
         NormalVelocityHTst = TstState->getNormalVelocityH(NewTime);
         for (int Cell = 0; Cell < NCellsAll; Cell++) {
            for (int Layer = 0; Layer < NVertLayers; Layer++) {
               LayerThickHRef(Cell, Layer) = LayerThickHDef(Cell, Layer) + 1;
               LayerThickHTst(Cell, Layer) = LayerThickHDef(Cell, Layer) + 1;
            }
         }

         for (int Edge = 0; Edge < NEdgesAll; Edge++) {
            for (int Layer = 0; Layer < NVertLayers; Layer++) {
               NormalVelocityHRef(Edge, Layer) =
                   NormalVelocityHDef(Edge, Layer) + 1;
               NormalVelocityHTst(Edge, Layer) =
                   NormalVelocityHDef(Edge, Layer) + 1;
            }
         }
         RefState->copyToDevice(NewTime);
         TstState->copyToDevice(NewTime);

         // Check initial values
         for (int N = 0; N <= 1; ++N) {
            Count1 = checkHost(RefState, TstState, N, N);
            Count2 = checkDevice(RefState, TstState, N, N);
            if (Count1 == 0 and Count2 == 0) {
               LOG_INFO(
                   "State: State compare (TimeLevel {}, NTimeLevels {}) PASS",
                   N, NTimeLevels);
            } else {
               RetVal += 1;
               LOG_INFO(
                   "State: State compare (TimeLevel {}, NTimeLevels {}) FAIL",
                   N, NTimeLevels);
            }
         }

         // Perform time level updates.
         for (int N = 1; N < NTimeLevels; ++N) {
            TstState->updateTimeLevels();

            // The time index represents the n + Ith level so the new
            // time n + 1 has index 1. Current time is 0. Previous
            // times (n-1, n-2, etc.) are represented by negative indices
            // (-1, -2, etc respectively) but cannot extend below
            // -(NTimeLevels-2). After an update the time indices shift to
            // one older level, but wraps around if they extend below the
            // lower limit.
            int NMin          = -(NTimeLevels - 2);
            int CurTimeUpdate = CurTime - N;
            int NewTimeUpdate = NewTime - N;
            if (CurTimeUpdate < NMin)
               CurTimeUpdate += NTimeLevels;
            if (NewTimeUpdate < NMin)
               NewTimeUpdate += NTimeLevels;

            // Check updated levels show up in the right place
            Count1 = checkHost(RefState, TstState, CurTime, CurTimeUpdate);
            Count2 = checkDevice(RefState, TstState, CurTime, CurTimeUpdate);

            if (Count1 == 0 and Count2 == 0) {
               LOG_INFO("State: NTimeLevels={} After update {} "
                        "Current level: PASS",
                        NTimeLevels, N);
            } else {
               RetVal += 1;
               LOG_INFO("State: NTimeLevels={} After update {} "
                        "Current level: FAIL",
                        NTimeLevels, N);
            }

            Count1 = checkHost(RefState, TstState, NewTime, NewTimeUpdate);
            Count2 = checkDevice(RefState, TstState, NewTime, NewTimeUpdate);

            if (Count1 == 0 and Count2 == 0) {
               LOG_INFO("State: NTimeLevels={} After update {} "
                        "New time level: PASS",
                        NTimeLevels, N);
            } else {
               RetVal += 1;
               LOG_INFO("State: NTimeLevels={} After update {} "
                        "New time level: FAIL",
                        NTimeLevels, N);
            }
         }

         // Erase the reference and test cases to prep for next cycle
         OceanState::erase("Reference");
         OceanState::erase("Test");
      }

      // Finalize Omega objects
      OceanState::clear();
      Tracers::clear();
      AuxiliaryState::clear();
      Tendencies::clear();
      TimeStepper::clear();
      HorzMesh::clear();
      VertCoord::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
      FieldGroup::clear();
      Field::clear();
      Dimension::clear();

      if (RetVal == 0)
         LOG_INFO("State: Successful completion");
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
