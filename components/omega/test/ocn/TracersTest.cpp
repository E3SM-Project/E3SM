//===-- Test driver for OMEGA Tracers -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA tracers class
///
/// This driver tests that the OMEGA tracers module.
//
//===-----------------------------------------------------------------------===/

#include "Tracers.h"

#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeStepper.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

// misc. variables for testing
const Real RefReal = 3.0;

//------------------------------------------------------------------------------
// The initialization routine for Tracers testing. It calls various
// init routines, including the creation of the default decomposition.

I4 initTracersTest() {

   I4 Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_ERROR("Tracers: Error reading config file");
      return Err;
   }

   // Initialize the default time stepper
   Err = TimeStepper::init1();
   if (Err != 0) {
      LOG_ERROR("Tracers: error initializing default time stepper");
      return Err;
   }

   // Initialize the IO system
   Err = IO::init(DefComm);
   if (Err != 0) {
      LOG_ERROR("Tracers: error initializing parallel IO");
      return Err;
   }

   // Create the default decomposition (initializes the decomposition)
   Err = Decomp::init();
   if (Err != 0) {
      LOG_ERROR("Tracers: error initializing default decomposition");
      return Err;
   }

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0) {
      LOG_ERROR("Tracers: error initializing default halo");
      return Err;
   }

   // Initialize the default mesh
   Err = HorzMesh::init();
   if (Err != 0) {
      LOG_ERROR("Tracers: error initializing default mesh");
      return Err;
   }

   return 0;
}

//------------------------------------------------------------------------------
// The test driver for Tracers infrastructure
//
int main(int argc, char *argv[]) {

   int RetVal = 0;
   I4 Err;
   I4 count = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      // Call initialization routine
      Err = initTracersTest();
      if (Err != 0)
         LOG_ERROR("Tracers: Error initializing");

      // Get MPI vars if needed
      MachEnv *DefEnv = MachEnv::getDefault();
      MPI_Comm Comm   = DefEnv->getComm();
      I4 MyTask       = DefEnv->getMyTask();
      I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster   = DefEnv->isMasterTask();

      HorzMesh *DefHorzMesh       = HorzMesh::getDefault();
      Decomp *DefDecomp           = Decomp::getDefault();
      Halo *DefHalo               = Halo::getDefault();
      TimeStepper *DefTimeStepper = TimeStepper::getDefault();

      // initialize Tracers infrastructure
      Err = Tracers::init();
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("Tracers: initialzation FAIL");
      }

      I4 NTracers    = Tracers::getNumTracers();
      I4 NCellsOwned = DefHorzMesh->NCellsOwned;
      I4 NCellsSize  = DefHorzMesh->NCellsSize;
      I4 NVertLevels = DefHorzMesh->NVertLevels;
      I4 NTimeLevels = DefTimeStepper->getNTimeLevels();

      // Get group names
      std::vector<std::string> GroupNames = Tracers::getGroupNames();

      // Check if "Base" group exists
      if (std::find(GroupNames.begin(), GroupNames.end(), "Base") !=
          GroupNames.end()) {
         LOG_INFO("Tracers: Group, 'Base', exists PASS");
      } else {
         RetVal += 1;
         LOG_ERROR("Tracers: Group, 'Base', does not exist FAIL");
      }

      // Check if "Debug" group exists
      if (std::find(GroupNames.begin(), GroupNames.end(), "Debug") !=
          GroupNames.end()) {
         LOG_INFO("Tracers: Group, 'Debug', exists PASS");
      } else {
         RetVal += 1;
         LOG_ERROR("Tracers: Group, 'Debug', does not exist FAIL");
      }

      I4 TotalLength = 0;

      for (std::string GroupName : GroupNames) {
         std::pair<I4, I4> GroupRange;
         Err = Tracers::getGroupRange(GroupRange, GroupName);

         if (Err != 0) {
            LOG_ERROR("Tracers: getGroupRange returns {} FAIL", Err);
            RetVal += 1;
         }

         auto [StartIndex, GroupLength] = GroupRange;

         TotalLength += GroupLength;

         // Check if a group contains more than one tracers
         if (GroupLength > 0) {
            LOG_INFO("Tracers: {} tracers retrieval PASS", GroupName);
         } else {
            RetVal += 1;
            LOG_ERROR("Tracers: {} tracers retrieval FAIL", GroupName);
         }

         // Check if tracer index is a member of the Group
         for (I4 TracerIndex = StartIndex;
              TracerIndex < StartIndex + GroupLength; ++TracerIndex) {
            if (Tracers::isGroupMemberByIndex(TracerIndex, GroupName)) {
               LOG_INFO("Tracers: {} group has the tracer index, {} PASS",
                        GroupName, TracerIndex);
            } else {
               RetVal += 1;
               LOG_ERROR(
                   "Tracers: {} group does not have the tracer index, {} FAIL",
                   GroupName, TracerIndex);
            }
         }

         // Check if tracer index:name mapping is correct
         for (I4 TracerIndex = StartIndex;
              TracerIndex < StartIndex + GroupLength; ++TracerIndex) {
            std::string TracerName;
            Err = Tracers::getName(TracerName, TracerIndex);
            if (Err != 0) {
               LOG_ERROR("Tracers: getName returns {} FAIL", Err);
               RetVal += 1;
            }

            I4 RetTracerIndex;
            Err = Tracers::getIndex(RetTracerIndex, TracerName);
            if (Err != 0) {
               LOG_ERROR("Tracers: getIndex returns {} FAIL", Err);
               RetVal += 1;
            }

            if (TracerIndex == RetTracerIndex) {
               LOG_INFO("Tracers: {} group tracer:name mapping for {} is "
                        "correct PASS",
                        GroupName, TracerName);
            } else {
               RetVal += 1;
               LOG_ERROR("Tracers: {} group tracer:name mapping for {} is not "
                         "correct FAIL",
                         GroupName, TracerName);
            }
         }

         // check if tracer has a corresponding field
         for (I4 TracerIndex = StartIndex;
              TracerIndex < StartIndex + GroupLength; ++TracerIndex) {
            auto TracerField = Tracers::getFieldByIndex(TracerIndex);

            if (TracerField) {
               LOG_INFO("Tracers: getFieldByIndex returns a field PASS");
            } else {
               RetVal += 1;
               LOG_ERROR("Tracers: getFieldByIndex returns nullptr FAIL");
            }
         }
      }

      // Check if total number of tracers is correct
      if (TotalLength == NTracers) {
         LOG_INFO("Tracers: getNumTracers() returns correct tracer size PASS");
      } else {
         RetVal += 1;
         LOG_ERROR(
             "Tracers: getNumTracers() returns incorrect tracer size FAIL");
      }

      // Reference host array of current time level
      HostArray3DReal RefHostArray =
          HostArray3DReal("RefHostArray", NTracers, NCellsSize, NVertLevels);

      // intialize tracer elements of all time levels
      for (I4 TimeLevel = 1; TimeLevel + NTimeLevels > 1; --TimeLevel) {
         HostArray3DReal TempHostArray;
         Err = Tracers::getAllHost(TempHostArray, TimeLevel);
         if (Err != 0) {
            LOG_ERROR("getAllHost(TempHostArray, TimeLevel) returns non-zero "
                      "code: {}",
                      Err);
            RetVal += 1;
         }

         for (I4 Tracer = 0; Tracer < NTracers; ++Tracer) {
            for (I4 Cell = 0; Cell < NCellsSize; Cell++) {
               for (I4 Vert = 0; Vert < NVertLevels; Vert++) {
                  TempHostArray(Tracer, Cell, Vert) =
                      RefReal + Tracer + Cell + Vert + TimeLevel;
                  if (TimeLevel == 1)
                     RefHostArray(Tracer, Cell, Vert) =
                         TempHostArray(Tracer, Cell, Vert);
               }
            }
         }
         Tracers::copyToDevice(TimeLevel);
      }

      // Reference device array of new time level
      Array3DReal RefArray =
          Array3DReal("RefArray", NTracers, NCellsSize, NVertLevels);

      Err = Tracers::getAll(RefArray, 1);
      if (Err != 0) {
         LOG_ERROR("getAll(RefArray, 1) returns non-zero code: {}", Err);
         RetVal += 1;
      }

      // deepCopy(RefArray, RefArray); TODO: remove this

      // Reference field data of all tracers
      std::vector<Array2DReal> RefFieldDataArray;

      // get field references of all tracers
      for (I4 Tracer = 0; Tracer < NTracers; ++Tracer) {
         auto TracerField = Tracers::getFieldByIndex(Tracer);
         RefFieldDataArray.push_back(TracerField->getDataArray<Array2DReal>());
      }

      // update time levels
      Tracers::updateTimeLevels();

      // getAll of current time level(0) should return the same to RefArray
      Array3DReal CurArray;
      Err = Tracers::getAll(CurArray, 0);
      if (Err != 0) {
         LOG_ERROR("getAll(CurArray, 0) returns non-zero code: {}", Err);
         RetVal += 1;
      }

      count = -1;

      // check if time level shift works
      parallelReduce(
          "reduce1", {NTracers, NCellsOwned, NVertLevels},
          KOKKOS_LAMBDA(I4 Tracer, I4 Cell, I4 Vert, I4 & Accum) {
             if (std::abs(CurArray(Tracer, Cell, Vert) -
                          RefArray(Tracer, Cell, Vert)) > 1e-9) {
                Accum++;
             }
          },
          count);

      if (count == 0) {
         LOG_INFO("Tracers: Tracer data match after updateTimeLevels() PASS");
      } else {
         RetVal += 1;
         LOG_ERROR("Tracers: Not all tracer data match after "
                   "updateTimeLevels():{} FAIL",
                   count);
      }

      // test getByName and getByIndex
      for (I4 Tracer = 0; Tracer < NTracers; ++Tracer) {
         std::string TracerName;
         Tracers::getName(TracerName, Tracer);

         Array2DReal CurTracer;
         Err = Tracers::getByName(CurTracer, 0, TracerName);
         if (Err != 0) {
            LOG_ERROR("getByName(CurTracer, 0, TracerName) returns non-zero "
                      "code: {}",
                      Err);
            RetVal += 1;
         }

         count = -1;

         parallelReduce(
             "reduce2", {NCellsOwned, NVertLevels},
             KOKKOS_LAMBDA(I4 Cell, I4 Vert, I4 & Accum) {
                if (std::abs(CurTracer(Cell, Vert) -
                             (RefReal + Tracer + Cell + Vert + 1)) > 1e-9) {
                   Accum++;
                }
             },
             count);

         if (count == 0) {
            LOG_INFO("Tracers: Tracer data from getByName match after "
                     "updateTimeLevels() PASS");
         } else {
            RetVal += 1;
            LOG_ERROR("Tracers: Not all tracer data from getByName match after "
                      "updateTimeLevels():{} FAIL",
                      count);
         }
      }

      // test field data - this test should generate positive count value.
      // Referece and test field data are different due to updateTimeLevels()
      for (I4 Tracer = 0; Tracer < NTracers; ++Tracer) {

         auto TracerField          = Tracers::getFieldByIndex(Tracer);
         Array2DReal TestFieldData = TracerField->getDataArray<Array2DReal>();
         Array2DReal RefFieldData  = RefFieldDataArray[Tracer];

         count = -1;

         parallelReduce(
             "reduce3", {NCellsOwned, NVertLevels},
             KOKKOS_LAMBDA(I4 Cell, I4 Vert, I4 & Accum) {
                if (std::abs(RefFieldData(Cell, Vert) -
                             TestFieldData(Cell, Vert)) > 1e-9) {
                   Accum++;
                }
             },
             count);

         if (count > 0) {
            LOG_INFO("Tracers: Tracer field data correctly catch the "
                     "difference after updateTimeLevels() PASS");
         } else {
            RetVal += 1;
            LOG_ERROR("Tracers: Tracer field data should not match after "
                      "updateTimeLevels() FAIL");
         }
      }

      // update time levels to cycle back to original index
      for (I4 TimeLevel = 0; TimeLevel + NTimeLevels > 1; --TimeLevel) {
         // update time levels
         Tracers::updateTimeLevels();
      }

      // test field data - this test should generate zero count value
      // Referece and test field data are the same because updateTimeLevels() is
      // called NTimeLevels
      for (I4 Tracer = 0; Tracer < NTracers; ++Tracer) {
         auto TracerField          = Tracers::getFieldByIndex(Tracer);
         Array2DReal TestFieldData = TracerField->getDataArray<Array2DReal>();
         Array2DReal RefFieldData  = RefFieldDataArray[Tracer];

         count = -1;

         parallelReduce(
             "reduce4", {NCellsOwned, NVertLevels},
             KOKKOS_LAMBDA(I4 Cell, I4 Vert, I4 & Accum) {
                if (std::abs(RefFieldData(Cell, Vert) -
                             TestFieldData(Cell, Vert)) > 1e-9) {
                   Accum++;
                }
             },
             count);

         if (count == 0) {
            LOG_INFO("Tracers: Tracer field data correctly match after "
                     "updateTimeLevels() back to original index PASS");
         } else {
            RetVal += 1;
            LOG_ERROR("Tracers: Not all tracer data match after "
                      "updateTimeLevels() back to original index FAIL");
         }
      }

      count = 0;

      Array2DReal SaltTracerByName;
      Err = Tracers::getByName(SaltTracerByName, 1, "Salinity");

      Array2DReal SaltTracerByIndexVar;
      Err = Tracers::getByIndex(SaltTracerByIndexVar, 1, Tracers::IndxSalt);

      count = -1;

      parallelReduce(
          "reduce5", {NCellsOwned, NVertLevels},
          KOKKOS_LAMBDA(I4 Cell, I4 Vert, I4 & Accum) {
             if (std::abs(SaltTracerByName(Cell, Vert) -
                          SaltTracerByIndexVar(Cell, Vert)) > 1e-9) {
                Accum++;
             }
          },
          count);

      if (count == 0) {
         LOG_INFO("Tracers: Index variable IndxSalt correctly retrieves Tracer "
                  "data PASS");
      } else {
         RetVal += 1;
         LOG_ERROR("Tracers: Index variable IndxSalt incorrectly retrieves "
                   "Tracer data FAIL");
      }

      count = 0;

      // Finally, check if getHostByName returns the correct host array
      for (I4 Tracer = 0; Tracer < NTracers; ++Tracer) {
         std::string TracerName;
         Tracers::getName(TracerName, Tracer);

         HostArray2DReal TestHostArray;
         Err = Tracers::getHostByName(TestHostArray, 1, TracerName);
         if (Err != 0) {
            LOG_ERROR("getHostByName(TestHostArray, 1, TracerName) returns "
                      "non-zero code: {}",
                      Err);
            RetVal += 1;
         }

         for (I4 Cell = 0; Cell < NCellsOwned; Cell++) {
            for (I4 Vert = 0; Vert < NVertLevels; Vert++) {
               if (std::abs(RefHostArray(Tracer, Cell, Vert) -
                            TestHostArray(Cell, Vert)) > 1e-9)
                  ++count;
            }
         }
      }

      if (count == 0) {
         LOG_INFO("Tracers: Tracer getHostByName correctly retreive tracer "
                  "data PASS");
      } else {
         RetVal += 1;
         LOG_ERROR("Tracers: Tracer getHostByName retreives incorrect tracer "
                   "data FAIL");
      }

      Tracers::clear();
      TimeStepper::clear();
      HorzMesh::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
      FieldGroup::clear();
      Field::clear();
      Dimension::clear();

      if (RetVal == 0)
         LOG_INFO("Tracers: Successful completion");
   }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
