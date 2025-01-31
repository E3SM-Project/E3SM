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
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "mpi.h"
#include <chrono>
#include <thread>
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Set some constant reference values for simplicity
const I4 RefI4           = 3;
const I8 RefI8           = 400000000;
const R4 RefR4           = 5.1;
const R8 RefR8           = 6.123456789;
const std::string RefStr = "Reference String";

//------------------------------------------------------------------------------
// A simple test evaluation function
template <typename T>
void TestEval(const std::string &TestName, T TestVal, T ExpectVal, int &Error) {

   if (TestVal == ExpectVal) {
      LOG_INFO("{}: PASS", TestName);
   } else {
      LOG_ERROR("{}: FAIL", TestName);
      ++Error;
   }
}
//------------------------------------------------------------------------------
// Initialization routine to create reference Fields
int initIOStreamTest(Clock *&ModelClock // Model clock
) {

   int Err    = 0;
   int Err1   = 0;
   int ErrRef = 0;

   // Initialize machine environment and logging
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();
   initLogging(DefEnv);

   // Read the model configuration
   Config("Omega");
   Err1 = Config::readAll("omega.yml");
   TestEval("Config read all", Err1, ErrRef, Err);
   Config *OmegaConfig = Config::getOmegaConfig();

   // Initialize the default time stepper (phase 1) that includes the
   // num time levels, calendar, model clock and start/stop times and alarms
   Err = TimeStepper::init1();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error phase 1 initializing default time stepper");
      return Err;
   }

   // Use alternative model clock rather than the input config
   // for testing
   TimeInstant SimStartTime(0001, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(2, TimeUnits::Hours);
   ModelClock = new Clock(SimStartTime, TimeStep);

   // Initialize base-level IO
   Err1 = IO::init(DefComm);
   TestEval("IO Initialization", Err1, ErrRef, Err);

   // Initialize decomposition
   Decomp::init();
   Decomp *DefDecomp = Decomp::getDefault();

   // Initialize Halo updates
   Halo::init();
   OMEGA::Halo *DefHalo = OMEGA::Halo::getDefault();

   // Initialize Field
   Err1 = Field::init(ModelClock);
   TestEval("IO Field initialization", Err1, ErrRef, Err);

   // Initialize IOStreams
   Err1 = IOStream::init(ModelClock);
   TestEval("IOStream Initialization", Err1, ErrRef, Err);

   // Initialize HorzMesh - this should read Mesh stream
   Err1 = HorzMesh::init();
   TestEval("Horizontal mesh initialization", Err1, ErrRef, Err);
   HorzMesh *DefMesh = HorzMesh::getDefault();
   I4 NCellsSize     = DefMesh->NCellsSize;

   // Set vertical levels and time levels
   I4 NVertLevels = 60;
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);

   // Initialize State
   Err1 = OceanState::init();
   TestEval("Ocean state initialization", Err1, ErrRef, Err);

   // Initialize Aux State
   Err1 = AuxiliaryState::init();
   TestEval("Ocean aux state initialization", Err1, ErrRef, Err);

   // Initialize Tracers
   Err1 = Tracers::init();
   TestEval("Ocean tracer initialization", Err1, ErrRef, Err);

   // Add some global (Model and Simulation) metadata
   std::shared_ptr<Field> CodeField = Field::get(CodeMeta);
   std::shared_ptr<Field> SimField  = Field::get(SimMeta);

   Err1 = CodeField->addMetadata("CodeIntTest", 3);
   TestEval("Add code metadata int", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeRealTest", 4.567);
   TestEval("Add code metadata real", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeBoolTest", true);
   TestEval("Add code metadata bool", Err1, ErrRef, Err);
   std::string CodeStrVal = "ASampleString";
   Err1                   = CodeField->addMetadata("CodeStrTest", CodeStrVal);
   TestEval("Add code metadata str", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeVersion", "V0.0");
   TestEval("Add code metadata str literal", Err1, ErrRef, Err);
   Err1 = SimField->addMetadata("ExpName", "IOStreamsTest");
   TestEval("Add ExpName metadata", Err1, ErrRef, Err);
   std::string StartTimeStr = SimStartTime.getString(4, 2, "_");
   Err1 = SimField->addMetadata("SimStartTime", StartTimeStr);
   TestEval("Add SimStartTime metadata", Err1, ErrRef, Err);

   // Validate all streams (Mesh stream already validated in HorzMesh?)
   bool AllValidated = IOStream::validateAll();
   TestEval("IOStream Validation", AllValidated, true, Err);

   // End of init
   return Err;

} // End initialization IOStream test

//------------------------------------------------------------------------------
// We will test the IOStream interfaces by defining Fields and reading stream
// configurations during init. We will test periodic writing of data and
// reading/writing restart and initial data.

int main(int argc, char **argv) {

   int Err    = 0;
   int Err1   = 0;
   int ErrRef = 0;

   // Initialize the global MPI and Kokkos environments
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      Clock *ModelClock = nullptr;

      // Call initialization to create reference IO field
      Err1 = initIOStreamTest(ModelClock);
      TestEval("Initialize IOStream test", Err1, ErrRef, Err);

      // Retrieve dimension lengths and some mesh info
      I4 NCellsSize     = Dimension::getDimLengthLocal("NCells");
      I4 NVertLevels    = Dimension::getDimLengthLocal("NVertLevels");
      Decomp *DefDecomp = Decomp::getDefault();
      I4 NCellsOwned    = DefDecomp->NCellsOwned;
      Array1DI4 CellID  = DefDecomp->CellID;

      // Read restart file for initial temperature and salinity data
      Metadata ReqMetadata; // leave empty for now - no required metadata
      Err1 = IOStream::read("InitialState", ModelClock, ReqMetadata);
      TestEval("Read restart file", Err1, IOStream::Success, Err);

      // Overwrite salinity array with values associated with global cell
      // ID to test proper indexing of IO
      Array2DReal Test("Test", NCellsSize, NVertLevels);
      Array2DReal Salt;
      Err1 = Tracers::getByIndex(Salt, 0, Tracers::IndxSalt);
      TestEval("Retrieve Salinity", Err1, ErrRef, Err);

      parallelFor(
          {NCellsSize, NVertLevels}, KOKKOS_LAMBDA(int Cell, int K) {
             Salt(Cell, K) = 0.0001_Real * (CellID(Cell) + K);
             Test(Cell, K) = Salt(Cell, K);
          });

      // Create a stop alarm at 1 year for time stepping
      TimeInstant StopTime(0002, 1, 1, 0, 0, 0.0);
      Alarm StopAlarm("Stop Time", StopTime);
      Err1 = ModelClock->attachAlarm(&StopAlarm);

      // Overwrite
      // Step forward in time and write files if it is time
      while (!StopAlarm.isRinging()) {
         ModelClock->advance();
         TimeInstant CurTime    = ModelClock->getCurrentTime();
         std::string CurTimeStr = CurTime.getString(4, 2, " ");

         Err1 = IOStream::writeAll(ModelClock);
         if (Err1 != 0) // to prevent too much output in log
            TestEval("Write all streams " + CurTimeStr, Err1, ErrRef, Err);
      }

      // Force read the latest restart and check the results
      // A delay is needed here to make sure the restart file is completely
      // written before we read.
      std::this_thread::sleep_for(std::chrono::seconds(5));
      bool ForceRead = true;
      Err1 = IOStream::read("RestartRead", ModelClock, ReqMetadata, ForceRead);
      TestEval("Restart force read", Err1, IOStream::Success, Err);

      Err1             = 0;
      auto DataReducer = Kokkos::Sum<I4>(Err1);

      parallelReduce(
          {NCellsOwned, NVertLevels},
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
   Tracers::clear();
   OceanState::clear();
   AuxiliaryState::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   Kokkos::finalize();
   MPI_Finalize();

   if (Err >= 256)
      Err = 255;

   // End of testing
   return Err;
}
//===--- End test driver for IO Streams ------------------------------------===/
