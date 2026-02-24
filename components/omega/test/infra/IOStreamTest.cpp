//===-- Test driver for OMEGA IO Streams -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA IO Streams
///
/// This driver tests the ability to write to/from files within Omega. An
/// IOStream is defined for each unique combination of read/write frequency,
/// contents of the file and other properties.
//
//===-----------------------------------------------------------------------===/

#include "IOStream.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "Pacer.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"
#include <chrono>
#include <thread>
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
// A simple test evaluation function
template <typename T>
void TestEval(const std::string &TestName, T TestVal, T ExpectVal, Error &Err) {
   if (TestVal != ExpectVal) {
      std::string ErrMsg = TestName + ": FAIL";
      Err += Error(ErrorCode::Fail, ErrMsg);
   }
}
//------------------------------------------------------------------------------
// Initialization routine to create reference Fields
void initIOStreamTest(Clock *&ModelClock // Model clock
) {

   Error Err;
   int TmpErr = 0;

   // Initialize machine environment and logging
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();
   initLogging(DefEnv);
   LOG_INFO("------ IOStream Unit Tests ------");

   // Read the model configuration
   Config("Omega");
   Config::readAll("omega.yml");
   Config *OmegaConfig = Config::getOmegaConfig();

   // Initialize the default time stepper (phase 1) that includes the
   // num time levels, calendar, model clock and start/stop times and alarms
   TimeStepper::init1();

   // Use alternative model clock rather than the input config
   // for testing
   TimeInstant SimStartTime(0001, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(2, TimeUnits::Hours);
   ModelClock = new Clock(SimStartTime, TimeStep);

   // Initialize base-level IO
   IO::init(DefComm);

   // Initialize decomposition
   Decomp::init();
   Decomp *DefDecomp = Decomp::getDefault();

   // Initialize Halo updates
   Halo::init();
   OMEGA::Halo *DefHalo = OMEGA::Halo::getDefault();

   // Initialize Field
   Field::init(ModelClock);

   // Initialize IOStreams
   IOStream::init(ModelClock);

   // Initialize HorzMesh - this should read Mesh stream
   HorzMesh::init();
   HorzMesh *DefMesh = HorzMesh::getDefault();

   // Initialize the vertical coordinate
   VertCoord::init();

   // Initialize State
   TmpErr = OceanState::init();
   if (TmpErr != 0)
      ABORT_ERROR("IOStreamTest: Error initializing OceanState");

   PressureGrad::init();

   // Intialize Tendencies
   Tendencies::init();

   // Initialize Tracers
   Tracers::init();

   // Initialize Aux State
   AuxiliaryState::init();

   // Add some global (Model and Simulation) metadata
   std::shared_ptr<Field> CodeField = Field::get(CodeMeta);
   std::shared_ptr<Field> SimField  = Field::get(SimMeta);

   CodeField->addMetadata("CodeIntTest", 3);
   CodeField->addMetadata("CodeRealTest", 4.567);
   CodeField->addMetadata("CodeBoolTest", true);
   std::string CodeStrVal = "ASampleString";
   CodeField->addMetadata("CodeStrTest", CodeStrVal);
   CodeField->addMetadata("CodeVersion", "V0.0");
   SimField->addMetadata("ExpName", "IOStreamsTest");
   std::string StartTimeStr = SimStartTime.getString(4, 2, "_");
   SimField->addMetadata("SimStartTime", StartTimeStr);

   // Validate all streams (Mesh stream already validated in HorzMesh?)
   bool AllValidated = IOStream::validateAll();
   TestEval("IOStream Validation", AllValidated, true, Err);

   // End of init
   CHECK_ERROR_ABORT(Err, "IOStreamTest: Error in initialization");
   return;

} // End initialization IOStream test

//------------------------------------------------------------------------------
// We will test the IOStream interfaces by defining Fields and reading stream
// configurations during init. We will test periodic writing of data and
// reading/writing restart and initial data.

int main(int argc, char **argv) {

   Error Err; // internal error code
   int Err1   = 0;
   int ErrRef = 0;

   // Initialize the global MPI and Kokkos environments
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {

      Clock *ModelClock = nullptr;

      // Call initialization to create reference IO field
      initIOStreamTest(ModelClock);

      // Retrieve dimension lengths and some mesh info
      I4 NCellsSize     = Dimension::getDimLengthLocal("NCells");
      I4 NVertLayers    = Dimension::getDimLengthLocal("NVertLayers");
      Decomp *DefDecomp = Decomp::getDefault();
      I4 NCellsOwned    = DefDecomp->NCellsOwned;
      Array1DI4 CellID  = DefDecomp->CellID;

      // Read restart file for initial temperature and salinity data
      Metadata ReqMetadata; // leave empty for now - no required metadata
      Err = IOStream::read("InitialState", ModelClock, ReqMetadata);
      CHECK_ERROR_ABORT(Err, "IOStreamTest: Error reading initial state");

      // Overwrite salinity array with values associated with global cell
      // ID to test proper indexing of IO
      Array2DReal Test("Test", NCellsSize, NVertLayers);
      Array2DReal Salt = Tracers::getByIndex(0, Tracers::IndxSalt);

      parallelFor(
          {NCellsSize, NVertLayers}, KOKKOS_LAMBDA(int Cell, int K) {
             Salt(Cell, K) = 0.0001_Real * (CellID(Cell) + K);
             Test(Cell, K) = Salt(Cell, K);
          });

      // Create a stop alarm at 1 year for time stepping
      TimeInstant StopTime(0002, 1, 1, 0, 0, 0.0);
      Alarm StopAlarm("Stop Time", StopTime);
      ModelClock->attachAlarm(&StopAlarm);

      // Overwrite
      // Step forward in time and write files if it is time
      while (!StopAlarm.isRinging()) {
         ModelClock->advance();
         TimeInstant CurTime    = ModelClock->getCurrentTime();
         std::string CurTimeStr = CurTime.getString(4, 2, " ");

         IOStream::writeAll(ModelClock);
      }

      // Force read the latest restart and check the results
      // A delay is needed here to make sure the restart file is completely
      // written before we read.
      std::this_thread::sleep_for(std::chrono::seconds(5));
      bool ForceRead = true;
      Err = IOStream::read("RestartRead", ModelClock, ReqMetadata, ForceRead);
      CHECK_ERROR_ABORT(Err, "IOStreamTest: Error reading restart");

      Err1             = 0;
      auto DataReducer = Kokkos::Sum<I4>(Err1);

      parallelReduce(
          {NCellsOwned, NVertLayers},
          KOKKOS_LAMBDA(int Cell, int K, I4 &Err1) {
             if (Salt(Cell, K) != Test(Cell, K))
                ++Err1;
          },
          DataReducer);
      TestEval("Check Salt array ", Err1, ErrRef, Err);

      // Write final output and remove all streams
      IOStream::finalize(ModelClock);
   }

   // Clean up environments
   TimeStepper::clear();
   PressureGrad::clear();
   Tendencies::clear();
   Tracers::clear();
   OceanState::clear();
   AuxiliaryState::clear();
   HorzMesh::clear();
   VertCoord::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   CHECK_ERROR_ABORT(Err, "IOStream Unit Tests: FAIL");
   LOG_INFO("------ IOStream Unit Tests successful ------");

   // End of testing
   return 0;
}
//===--- End test driver for IO Streams ------------------------------------===/
