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
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TimeStepper.h"
#include "mpi.h"

#include <iostream>

//------------------------------------------------------------------------------
// The initialization routine for State testing. It calls various
// init routines, including the creation of the default decomposition.

int initStateTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
   MPI_Comm DefComm       = DefEnv->getComm();

   OMEGA::initLogging(DefEnv);

   // Open config file
   OMEGA::Config("Omega");
   Err = OMEGA::Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("State: Error reading config file");
      return Err;
   }

   // Initialize the default time stepper
   Err = OMEGA::TimeStepper::init1();
   if (Err != 0)
      LOG_ERROR("State: error initializing default time stepper");

   // Initialize the IO system
   Err = OMEGA::IO::init(DefComm);
   if (Err != 0)
      LOG_ERROR("State: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Err = OMEGA::Decomp::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default decomposition");

   // Initialize the default halo
   Err = OMEGA::Halo::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default halo");

   // Initialize the default mesh
   Err = OMEGA::HorzMesh::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default mesh");

   return Err;
}

// Check for differences between layer thickness and normal velocity host arrays
int checkHost(OMEGA::OceanState *DefState, OMEGA::OceanState *TestState,
              int DefLevel, int TestLevel) {

   int count = 0;
   OMEGA::HostArray2DReal LayerThicknessH_def;
   OMEGA::HostArray2DReal LayerThicknessH_test;
   DefState->getLayerThicknessH(LayerThicknessH_def, DefLevel);
   TestState->getLayerThicknessH(LayerThicknessH_test, TestLevel);
   for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
      for (int Level = 0; Level < DefState->NVertLevels; Level++) {
         if (LayerThicknessH_def(Cell, Level) !=
             LayerThicknessH_test(Cell, Level)) {
            count++;
         }
      }
   }

   OMEGA::HostArray2DReal NormalVelocityH_def;
   OMEGA::HostArray2DReal NormalVelocityH_test;
   DefState->getNormalVelocityH(NormalVelocityH_def, DefLevel);
   TestState->getNormalVelocityH(NormalVelocityH_test, TestLevel);
   for (int Edge = 0; Edge < DefState->NEdgesAll; Edge++) {
      for (int Level = 0; Level < DefState->NVertLevels; Level++) {
         if (NormalVelocityH_def(Edge, Level) !=
             NormalVelocityH_test(Edge, Level)) {
            count++;
         }
      }
   }

   return count;
}

// Check for differences between layer thickness and normal velocity device
// arrays
int checkDevice(OMEGA::OceanState *DefState, OMEGA::OceanState *TestState,
                int DefLevel, int TestLevel) {

   int count1;
   OMEGA::Array2DReal LayerThickness_def;
   OMEGA::Array2DReal LayerThickness_test;
   DefState->getLayerThickness(LayerThickness_def, DefLevel);
   TestState->getLayerThickness(LayerThickness_test, TestLevel);
   OMEGA::parallelReduce(
       "reduce", {DefState->NCellsAll, DefState->NVertLevels},
       KOKKOS_LAMBDA(int Cell, int Level, int &Accum) {
          if (LayerThickness_def(Cell, Level) !=
              LayerThickness_test(Cell, Level)) {
             Accum++;
          }
       },
       count1);

   int count2;
   OMEGA::Array2DReal NormalVelocity_def;
   OMEGA::Array2DReal NormalVelocity_test;
   DefState->getNormalVelocity(NormalVelocity_def, DefLevel);
   TestState->getNormalVelocity(NormalVelocity_test, TestLevel);
   OMEGA::parallelReduce(
       "reduce", {DefState->NCellsAll, DefState->NVertLevels},
       KOKKOS_LAMBDA(int Cell, int Level, int &Accum) {
          if (NormalVelocity_def(Cell, Level) !=
              NormalVelocity_test(Cell, Level)) {
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

      // Call initialization routine to create the default decomposition
      int Err = initStateTest();
      if (Err != 0)
         LOG_CRITICAL("State: Error initializing");

      // Get MPI vars if needed
      OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
      MPI_Comm Comm          = DefEnv->getComm();
      OMEGA::I4 MyTask       = DefEnv->getMyTask();
      OMEGA::I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster          = DefEnv->isMasterTask();

      OMEGA::HorzMesh *DefHorzMesh = OMEGA::HorzMesh::getDefault();
      OMEGA::Decomp *DefDecomp     = OMEGA::Decomp::getDefault();
      OMEGA::Halo *DefHalo         = OMEGA::Halo::getDefault();

      // These hard-wired variables need to be upated
      // with retrivals/config options
      int NVertLevels = 60;
      int NTimeLevels = 2;

      // Create dimensions (Horz dims computed in Mesh init)
      auto VertDim = OMEGA::Dimension::create("NVertLevels", NVertLevels);

      for (int NTimeLevels = 2; NTimeLevels < 5; NTimeLevels++) {

         int NewLevel = 1;
         int CurLevel = 0;

         // Create "default" state
         if (NTimeLevels == 2) {

            OMEGA::OceanState::init();
            OMEGA::OceanState *DefOceanState = OMEGA::OceanState::getDefault();

         } else {
            OMEGA::OceanState *DefState = OMEGA::OceanState::create(
                "Default", DefHorzMesh, DefHalo, NVertLevels, NTimeLevels);
            DefState->loadStateFromFile(DefHorzMesh->MeshFileName, DefDecomp);
         }

         // Test retrieval of the default state
         OMEGA::OceanState *DefState = OMEGA::OceanState::get("Default");
         if (DefState) { // true if non-null ptr
            LOG_INFO("State: Default state retrieval (NTimeLevels={}) PASS",
                     NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO("State: Default state retrieval (NTimeLevels={}) FAIL",
                     NTimeLevels);
         }

         // Create "test" state
         OMEGA::OceanState::create("Test", DefHorzMesh, DefHalo, NVertLevels,
                                   NTimeLevels);

         OMEGA::OceanState *TestState = OMEGA::OceanState::get("Test");

         if (TestState) { // true if non-null ptr
            LOG_INFO("State: Test state retrieval (NTimeLevels={}) PASS",
                     NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO("State: Test state retrieval (NTimeLevels={}) FAIL",
                     NTimeLevels);
         }

         // Initially fill test state with the same values as the default state
         TestState->loadStateFromFile(DefHorzMesh->MeshFileName, DefDecomp);

         // Test that reasonable values have been read in for LayerThickness
         OMEGA::HostArray2DReal LayerThickH;
         DefState->getLayerThicknessH(LayerThickH, CurLevel);
         int count = 0;
         for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
            int colCount = 0;
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               OMEGA::R8 val = LayerThickH(Cell, Level);
               if (val > 0.0 && val < 300.0) {
                  colCount++;
               }
            }
            if (colCount < 2) {
               count++;
            }
         }

         if (count == 0) {
            LOG_INFO("State: State read (NTimeLevels={}) PASS", NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO("State: State read (NTimeLevels={}) FAIL", NTimeLevels);
         }

         // Initialize NormalVelocity values
         OMEGA::HostArray2DReal NormalVelocityHDef;
         OMEGA::HostArray2DReal NormalVelocityHTest;
         DefState->getNormalVelocityH(NormalVelocityHDef, CurLevel);
         TestState->getNormalVelocityH(NormalVelocityHTest, CurLevel);
         for (int Edge = 0; Edge < DefState->NEdgesAll; Edge++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               NormalVelocityHDef(Edge, Level)  = Edge;
               NormalVelocityHTest(Edge, Level) = Edge;
            }
         }
         DefState->exchangeHalo(CurLevel);
         TestState->exchangeHalo(CurLevel);

         // Test that initally the 0 time levels of the
         // Def and Test state arrays match
         int count1 = checkHost(DefState, TestState, CurLevel, CurLevel);
         DefState->copyToDevice(CurLevel);
         TestState->copyToDevice(CurLevel);
         int count2 = checkDevice(DefState, TestState, CurLevel, CurLevel);

         if (count1 + count2 == 0) {
            LOG_INFO(
                "State: Default test state comparison (NTimeLevels={}) PASS",
                NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO(
                "State: Default test state comparison (NTimeLevels={}) FAIL",
                NTimeLevels);
         }

         // Perform time level update.
         DefState->updateTimeLevels();

         // Test that the time level update is correct.
         // Time levels should be different after one update
         count1 = checkHost(DefState, TestState, CurLevel, CurLevel);
         count2 = checkDevice(DefState, TestState, CurLevel, CurLevel);

         if (count1 + count2 != 0) {
            LOG_INFO("State: time levels different after single update "
                     "(NTimeLevels={}) PASS",
                     NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO("State: time levels different after single update "
                     "(NTimeLevels={}) FAIL",
                     NTimeLevels);
         }

         // Perform time level updates to cycle back to inital index
         for (int i = 0; i < NTimeLevels - 1; i++) {
            DefState->updateTimeLevels();
         }

         // Test that the time level update is correct.
         // Time levels should be the same again
         count1 = checkHost(DefState, TestState, CurLevel, CurLevel);
         count2 = checkDevice(DefState, TestState, CurLevel, CurLevel);

         if (count1 + count2 == 0) {
            LOG_INFO("State: time level update (NTimeLevels={}) PASS",
                     NTimeLevels);
         } else {
            RetVal += 1;
            LOG_INFO("State: time level update (NTimeLevels={}) FAIL",
                     NTimeLevels);
         }

         OMEGA::OceanState::clear();
      }

      // Finalize Omega objects
      OMEGA::TimeStepper::clear();
      OMEGA::HorzMesh::clear();
      OMEGA::Decomp::clear();
      OMEGA::MachEnv::removeAll();
      OMEGA::FieldGroup::clear();
      OMEGA::Field::clear();
      OMEGA::Dimension::clear();

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
