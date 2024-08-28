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

   // Initialize the default time stepper
   Err = OMEGA::TimeStepper::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default time stepper");

   return Err;
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

      for (int NTimeLevels = 2; NTimeLevels < 4; NTimeLevels++) {

         int CurLevel = NTimeLevels - 2;
         int NewLevel = CurLevel - 1;
         if (NewLevel < 0) {
            NewLevel = NTimeLevels - 1;
         }

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
            LOG_INFO("State: Default state retrieval PASS");
         } else {
            RetVal += 1;
            LOG_INFO("State: Default state retrieval FAIL");
         }

         // Create "test" state
         OMEGA::OceanState::create("Test", DefHorzMesh, DefHalo, NVertLevels,
                                   NTimeLevels);

         OMEGA::OceanState *TestState = OMEGA::OceanState::get("Test");

         if (TestState) { // true if non-null ptr
            LOG_INFO("State: Test state retrieval PASS");
         } else {
            RetVal += 1;
            LOG_INFO("State: Test state retrieval FAIL");
         }

         // Initially fill test state with the same values as the default state
         TestState->loadStateFromFile(DefHorzMesh->MeshFileName, DefDecomp);

         // Test that reasonable values have been read in for LayerThickness
         int count = 0;
         for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
            int colCount = 0;
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               OMEGA::R8 val = DefState->LayerThicknessH[CurLevel](Cell, Level);
               if (val > 0.0 && val < 300.0) {
                  colCount++;
               }
            }
            if (colCount < 2) {
               count++;
            }
         }

         if (count == 0) {
            LOG_INFO("State: State read PASS");
         } else {
            RetVal += 1;
            LOG_INFO("State: State read FAIL");
         }

         // Initialize NormalVelocity values
         for (int Edge = 0; Edge < DefState->NEdgesAll; Edge++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               DefState->NormalVelocityH[CurLevel](Edge, Level)  = Edge;
               TestState->NormalVelocityH[CurLevel](Edge, Level) = Edge;
            }
         }

         // Test that initally the 0 time levels of the
         // Def and Test state arrays match
         count                     = 0;
         auto LayerThicknessH_def  = DefState->LayerThicknessH[CurLevel];
         auto LayerThicknessH_test = TestState->LayerThicknessH[CurLevel];
         for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               if (LayerThicknessH_def(Cell, Level) !=
                   LayerThicknessH_test(Cell, Level)) {
                  count++;
               }
            }
         }

         auto NormalVelocityH_def  = DefState->NormalVelocityH[CurLevel];
         auto NormalVelocityH_test = TestState->NormalVelocityH[CurLevel];
         for (int Edge = 0; Edge < DefState->NEdgesAll; Edge++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               if (NormalVelocityH_def(Edge, Level) !=
                   NormalVelocityH_test(Edge, Level)) {
                  count++;
               }
            }
         }

         if (count == 0) {
            LOG_INFO("State: Default test state comparison PASS");
         } else {
            RetVal += 1;
            LOG_INFO("State: Default test state comparison FAIL");
         }

         // Perform time level update.
         DefState->updateTimeLevels();

         // Test that the time level update is correct.
         count                = 0;
         LayerThicknessH_def  = DefState->LayerThicknessH[NewLevel];
         LayerThicknessH_test = TestState->LayerThicknessH[CurLevel];
         for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               if (LayerThicknessH_def(Cell, Level) !=
                   LayerThicknessH_test(Cell, Level)) {
                  count++;
               }
            }
         }

         LayerThicknessH_def  = DefState->LayerThicknessH[CurLevel];
         LayerThicknessH_test = TestState->LayerThicknessH[NewLevel];
         for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               if (LayerThicknessH_def(Cell, Level) !=
                   LayerThicknessH_test(Cell, Level)) {
                  count++;
               }
            }
         }

         NormalVelocityH_def  = DefState->NormalVelocityH[NewLevel];
         NormalVelocityH_test = TestState->NormalVelocityH[CurLevel];
         for (int Edge = 0; Edge < DefState->NEdgesAll; Edge++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               if (NormalVelocityH_def(Edge, Level) !=
                   NormalVelocityH_test(Edge, Level)) {
                  count++;
               }
            }
         }

         NormalVelocityH_def  = DefState->NormalVelocityH[CurLevel];
         NormalVelocityH_test = TestState->NormalVelocityH[NewLevel];
         for (int Edge = 0; Edge < DefState->NEdgesAll; Edge++) {
            for (int Level = 0; Level < DefState->NVertLevels; Level++) {
               if (NormalVelocityH_def(Edge, Level) !=
                   NormalVelocityH_test(Edge, Level)) {
                  count++;
               }
            }
         }

         if (count == 0) {
            LOG_INFO("State: time level update PASS");
         } else {
            RetVal += 1;
            LOG_INFO("State: time level update FAIL");
         }

         // Test time level update on device
         int count1;
         auto LayerThickness_def  = DefState->LayerThickness[NewLevel];
         auto LayerThickness_test = TestState->LayerThickness[CurLevel];
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
         LayerThickness_def  = DefState->LayerThickness[CurLevel];
         LayerThickness_test = TestState->LayerThickness[NewLevel];
         OMEGA::parallelReduce(
             "reduce", {DefState->NCellsAll, DefState->NVertLevels},
             KOKKOS_LAMBDA(int Cell, int Level, int &Accum) {
                if (LayerThickness_def(Cell, Level) !=
                    LayerThickness_test(Cell, Level)) {
                   Accum++;
                }
             },
             count2);

         int count3;
         auto NormalVelocity_def  = DefState->NormalVelocity[CurLevel];
         auto NormalVelocity_test = TestState->NormalVelocity[NewLevel];
         OMEGA::parallelReduce(
             "reduce", {DefState->NEdgesAll, DefState->NVertLevels},
             KOKKOS_LAMBDA(int Edge, int Level, int &Accum) {
                if (NormalVelocity_def(Edge, Level) !=
                    NormalVelocity_test(Edge, Level)) {
                   Accum++;
                }
             },
             count3);

         int count4;
         NormalVelocity_def  = DefState->NormalVelocity[NewLevel];
         NormalVelocity_test = TestState->NormalVelocity[CurLevel];
         OMEGA::parallelReduce(
             "reduce", {DefState->NEdgesAll, DefState->NVertLevels},
             KOKKOS_LAMBDA(int Edge, int Level, int &Accum) {
                if (NormalVelocity_def(Edge, Level) !=
                    NormalVelocity_test(Edge, Level)) {
                   Accum++;
                }
             },
             count4);

         if (count1 == 0 && count2 == 0 && count3 == 0 && count4 == 0) {
            LOG_INFO("State: time level update (GPU) PASS");
         } else {
            RetVal += 1;
            LOG_INFO("State: time level update (GPU) FAIL");
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
