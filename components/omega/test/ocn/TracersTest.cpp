//===-- Test driver for OMEGA Tracers -----------------------------*- C++-*-===/
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
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "VertCoord.h"
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
   Config::readAll("omega.yml");

   // Initialize the default time stepper and model clock
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0) {
      LOG_ERROR("Tracers: error initializing default halo");
      return Err;
   }

   // Read in the default mesh
   Field::init(ModelClock);
   IOStream::init(ModelClock);
   HorzMesh::init(ModelClock);

   // Initialize the vertical coordinate
   VertCoord::init(false);

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
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
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
      Tracers::init();

      I4 NTracers    = Tracers::getNumTracers();
      I4 NCellsOwned = DefHorzMesh->NCellsOwned;
      I4 NCellsSize  = DefHorzMesh->NCellsSize;
      I4 NVertLayers = VertCoord::getDefault()->NVertLayers;
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
          HostArray3DReal("RefHostArray", NTracers, NCellsSize, NVertLayers);

      // intialize tracer elements of all time levels
      for (I4 TimeLevel = 1; TimeLevel + NTimeLevels > 1; --TimeLevel) {
         HostArray3DReal TempHostArray = Tracers::getAllHost(TimeLevel);

         for (I4 Tracer = 0; Tracer < NTracers; ++Tracer) {
            for (I4 Cell = 0; Cell < NCellsSize; Cell++) {
               for (I4 Vert = 0; Vert < NVertLayers; Vert++) {
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
          Array3DReal("RefArray", NTracers, NCellsSize, NVertLayers);

      RefArray = Tracers::getAll(1);

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
      Array3DReal CurArray = Tracers::getAll(0);

      count = -1;

      // check if time level shift works
      parallelReduce(
          "reduce1", {NTracers, NCellsOwned, NVertLayers},
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

         Array2DReal CurTracer = Tracers::getByName(0, TracerName);

         count = -1;

         parallelReduce(
             "reduce2", {NCellsOwned, NVertLayers},
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
             "reduce3", {NCellsOwned, NVertLayers},
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
             "reduce4", {NCellsOwned, NVertLayers},
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

      Array2DReal SaltTracerByName = Tracers::getByName(1, "Salinity");

      Array2DReal SaltTracerByIndexVar =
          Tracers::getByIndex(1, Tracers::IndxSalt);

      count = -1;

      parallelReduce(
          "reduce5", {NCellsOwned, NVertLayers},
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

         HostArray2DReal TestHostArray = Tracers::getHostByName(1, TracerName);

         for (I4 Cell = 0; Cell < NCellsOwned; Cell++) {
            for (I4 Vert = 0; Vert < NVertLayers; Vert++) {
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
      VertCoord::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
      FieldGroup::clear();
      Field::clear();
      Dimension::clear();

      if (RetVal == 0)
         LOG_INFO("Tracers: Successful completion");
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
