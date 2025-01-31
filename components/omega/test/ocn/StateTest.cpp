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
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TimeStepper.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The initialization routine for State testing. It calls various
// init routines, including the creation of the default decomposition.

int initStateTest() {

   int Err = 0;

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
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("State: Error reading config file");
      return Err;
   }

   // Initialize the default time stepper and retrieve the model clock
   Err = TimeStepper::init1();
   if (Err != 0)
      LOG_ERROR("State: error initializing default time stepper");

   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   Err = IO::init(DefComm);
   if (Err != 0)
      LOG_ERROR("State: error initializing parallel IO");

   // Initialize IOStreams - this does not yet validate the contents
   // of each file, only creates streams from Config
   Err = IOStream::init(ModelClock);
   if (Err != 0) {
      LOG_CRITICAL("State: Error initializing IOStreams");
      return Err;
   }

   // Initialize Field infrastructure
   Err = Field::init(ModelClock);
   if (Err != 0) {
      LOG_CRITICAL("State: Error initializing Fields");
      return Err;
   }

   // Create the default decomposition (initializes the decomposition)
   Err = Decomp::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default decomposition");

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default halo");

   // Initialize the default mesh
   Err = HorzMesh::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default mesh");

   // Get the number of vertical levels from Config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config DimConfig("Dimension");
   Err = OmegaConfig->get(DimConfig);
   if (Err != 0) {
      LOG_CRITICAL("State: Dimension group not found in Config");
      return Err;
   }
   I4 NVertLevels;
   Err = DimConfig.get("NVertLevels", NVertLevels);
   if (Err != 0) {
      LOG_CRITICAL("State: NVertLevels not found in Dimension Config");
      return Err;
   }
   // Create vertical dimension
   auto VertDim = Dimension::create("NVertLevels", NVertLevels);

   // Initialize tracers
   Err = Tracers::init();
   if (Err != 0) {
      LOG_CRITICAL("State: Error initializing tracers infrastructure");
      return Err;
   }

   // Initialize Aux State variables
   Err = AuxiliaryState::init();
   if (Err != 0) {
      LOG_CRITICAL("State: Error initializing default aux state");
      return Err;
   }

   // Create tendencies
   Err = Tendencies::init();
   if (Err != 0) {
      LOG_CRITICAL("State: Error initializing default tendencies");
      return Err;
   }

   // Finish time stepper initialization
   Err = TimeStepper::init2();
   if (Err != 0) {
      LOG_CRITICAL("State: Error phase 2 initializing default time stepper");
      return Err;
   }

   // Create a default ocean state
   Err = OceanState::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default state");
      return Err;
   }

   // Now that all fields have been defined, validate all the streams
   // contents
   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid) {
      LOG_CRITICAL("ocnInit: Error validating IO Streams");
      return Err;
   }

   // Read the state variables from the initial state stream
   Metadata ReqMeta; // no global metadata needed for init state read
   Err = IOStream::read("InitialState", ModelClock, ReqMeta);
   if (Err != IOStream::Success) {
      LOG_CRITICAL("ocnInit: Error reading initial state from stream");
      return Err;
   }

   // Finish initialization of initial state by filling halos and copying
   // to host. Current time level is zero.
   OceanState *DefState = OceanState::getDefault();
   DefState->exchangeHalo(0);
   DefState->copyToHost(0);

   return Err;
}

// Check for differences between layer thickness and normal velocity host arrays
int checkHost(OceanState *RefState, OceanState *TstState, int RefTimeLevel,
              int TstTimeLevel) {

   int count = 0;
   HostArray2DReal LayerThicknessHRef;
   HostArray2DReal LayerThicknessHTst;
   RefState->getLayerThicknessH(LayerThicknessHRef, RefTimeLevel);
   TstState->getLayerThicknessH(LayerThicknessHTst, TstTimeLevel);
   for (int Cell = 0; Cell < RefState->NCellsAll; Cell++) {
      for (int Level = 0; Level < RefState->NVertLevels; Level++) {
         if (LayerThicknessHRef(Cell, Level) !=
             LayerThicknessHTst(Cell, Level)) {
            count++;
         }
      }
   }

   HostArray2DReal NormalVelocityHRef;
   HostArray2DReal NormalVelocityHTst;
   RefState->getNormalVelocityH(NormalVelocityHRef, RefTimeLevel);
   TstState->getNormalVelocityH(NormalVelocityHTst, TstTimeLevel);
   for (int Edge = 0; Edge < RefState->NEdgesAll; Edge++) {
      for (int Level = 0; Level < RefState->NVertLevels; Level++) {
         if (NormalVelocityHRef(Edge, Level) !=
             NormalVelocityHTst(Edge, Level)) {
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
   Array2DReal LayerThicknessRef;
   Array2DReal LayerThicknessTst;
   RefState->getLayerThickness(LayerThicknessRef, RefTimeLevel);
   TstState->getLayerThickness(LayerThicknessTst, TstTimeLevel);
   parallelReduce(
       "reduce", {RefState->NCellsAll, RefState->NVertLevels},
       KOKKOS_LAMBDA(int Cell, int Level, int &Accum) {
          if (LayerThicknessRef(Cell, Level) !=
              LayerThicknessTst(Cell, Level)) {
             Accum++;
          }
       },
       count1);

   int count2;
   Array2DReal NormalVelocityRef;
   Array2DReal NormalVelocityTst;
   RefState->getNormalVelocity(NormalVelocityRef, RefTimeLevel);
   TstState->getNormalVelocity(NormalVelocityTst, TstTimeLevel);
   parallelReduce(
       "reduce", {RefState->NCellsAll, RefState->NVertLevels},
       KOKKOS_LAMBDA(int Cell, int Level, int &Accum) {
          if (NormalVelocityRef(Cell, Level) !=
              NormalVelocityTst(Cell, Level)) {
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
   {

      // Call initialization routine to create default state and other
      // quantities
      int Err = initStateTest();
      if (Err != 0)
         LOG_CRITICAL("State: Error initializing");

      // Get default mesh, halo and time data
      HorzMesh *DefHorzMesh   = HorzMesh::getDefault();
      Halo *DefHalo           = Halo::getDefault();
      TimeStepper *DefStepper = TimeStepper::getDefault();
      Clock *ModelClock       = DefStepper->getClock();
      int NTimeLevelsDef      = DefStepper->getNTimeLevels();

      // Test retrieval of the default state
      OceanState *DefState = OceanState::getDefault();
      int NVertLevels      = DefState->NVertLevels;
      int NCellsAll        = DefState->NCellsAll;
      int NEdgesAll        = DefState->NEdgesAll;
      int CurTime          = 0;
      int NewTime          = 1;
      if (DefState and NVertLevels > 0) {
         LOG_INFO("State: Default state retrieval PASS");
      } else {
         RetVal += 1;
         LOG_INFO("State: Default state retrieval FAIL");
      }

      // Get default state arrays and check for reasonable values
      HostArray2DReal LayerThickHDef;
      HostArray2DReal NormalVelocityHDef;
      Array2DReal LayerThickDef;
      Array2DReal NormalVelocityDef;
      DefState->getLayerThicknessH(LayerThickHDef, CurTime);
      DefState->getLayerThickness(LayerThickDef, CurTime);
      DefState->getNormalVelocityH(NormalVelocityHDef, CurTime);
      DefState->getNormalVelocity(NormalVelocityDef, CurTime);

      int Count1 = 0;
      int Count2 = 0;
      for (int Cell = 0; Cell < NCellsAll; Cell++) {
         for (int Level = 0; Level < NVertLevels; Level++) {
            R8 val = LayerThickHDef(Cell, Level);
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
             "Reference", DefHorzMesh, DefHalo, NVertLevels, NTimeLevels);
         OceanState *TstState = OceanState::create("Test", DefHorzMesh, DefHalo,
                                                   NVertLevels, NTimeLevels);

         // Test successful creation
         if (TstState and RefState) { // true if non-null ptr
            LOG_INFO("State: Test state creation (NTimeLevels={}) PASS",
                     NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO("State: Test state creation (NTimeLevels={}) FAIL",
                     NTimeLevels);
         }

         // Create test arrays
         HostArray2DReal LayerThickHRef;
         HostArray2DReal LayerThickHTst;
         HostArray2DReal NormalVelocityHRef;
         HostArray2DReal NormalVelocityHTst;

         // Fill reference and test states at current time with Default state
         RefState->getLayerThicknessH(LayerThickHRef, CurTime);
         TstState->getLayerThicknessH(LayerThickHTst, CurTime);
         RefState->getNormalVelocityH(NormalVelocityHRef, CurTime);
         TstState->getNormalVelocityH(NormalVelocityHTst, CurTime);
         for (int Cell = 0; Cell < NCellsAll; Cell++) {
            for (int Level = 0; Level < NVertLevels; Level++) {
               LayerThickHRef(Cell, Level) = LayerThickHDef(Cell, Level);
               LayerThickHTst(Cell, Level) = LayerThickHDef(Cell, Level);
            }
         }
         for (int Edge = 0; Edge < NEdgesAll; Edge++) {
            for (int Level = 0; Level < NVertLevels; Level++) {
               NormalVelocityHRef(Edge, Level) =
                   NormalVelocityHDef(Edge, Level);
               NormalVelocityHTst(Edge, Level) =
                   NormalVelocityHDef(Edge, Level);
            }
         }
         RefState->copyToDevice(CurTime);
         TstState->copyToDevice(CurTime);

         // Fill reference and test states at next time levels with the
         // Default state + 1
         RefState->getLayerThicknessH(LayerThickHRef, NewTime);
         TstState->getLayerThicknessH(LayerThickHTst, NewTime);
         RefState->getNormalVelocityH(NormalVelocityHRef, NewTime);
         TstState->getNormalVelocityH(NormalVelocityHTst, NewTime);
         for (int Cell = 0; Cell < NCellsAll; Cell++) {
            for (int Level = 0; Level < NVertLevels; Level++) {
               LayerThickHRef(Cell, Level) = LayerThickHDef(Cell, Level) + 1;
               LayerThickHTst(Cell, Level) = LayerThickHDef(Cell, Level) + 1;
            }
         }

         for (int Edge = 0; Edge < NEdgesAll; Edge++) {
            for (int Level = 0; Level < NVertLevels; Level++) {
               NormalVelocityHRef(Edge, Level) =
                   NormalVelocityHDef(Edge, Level) + 1;
               NormalVelocityHTst(Edge, Level) =
                   NormalVelocityHDef(Edge, Level) + 1;
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
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
      FieldGroup::clear();
      Field::clear();
      Dimension::clear();

      if (RetVal == 0)
         LOG_INFO("State: Successful completion");
   }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
