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
int initIOStreamTest(std::shared_ptr<Clock> &ModelClock // Model clock
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
   Err1 = Field::init();
   TestEval("IO Field initialization", Err1, ErrRef, Err);

   // Create the model clock and time step
   // Get Calendar from time management config group
   Config TimeMgmtConfig("TimeManagement");
   Err = OmegaConfig->get(TimeMgmtConfig);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: TimeManagement group not found in Config");
      return Err;
   }
   std::string ConfigCalStr;
   Err = TimeMgmtConfig.get("CalendarType", ConfigCalStr);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: CalendarType not found in TimeMgmtConfig");
      return Err;
   }
   Calendar::init(ConfigCalStr);

   // Use internal start time and time step rather than Config
   TimeInstant SimStartTime(0001, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(2, TimeUnits::Hours);
   ModelClock = std::make_shared<Clock>(SimStartTime, TimeStep);

   // Initialize IOStreams
   Err1 = IOStream::init(*ModelClock);
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

   // Initialize time stepper needed before ocean state (for time levels)
   Err1 = TimeStepper::init1();
   TestEval("Ocean time step initialization", Err1, ErrRef, Err);

   // Initialize State
   Err1 = OceanState::init();
   TestEval("Ocean state initialization", Err1, ErrRef, Err);

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

   // Define temperature and salinity tracer fields and create a tracer
   // group

   std::vector<std::string> DimNames(2);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";

   // 2D Fields on device

   DimNames[0]    = "NCells";
   DimNames[1]    = "NVertLevels";
   Real FillValue = -1.2345e-30;
   auto TempField = Field::create(
       "Temperature", "Potential temperature at cell centers", "deg C",
       "sea_water_pot_tem", -3.0, 100.0, FillValue, 2, DimNames);
   auto SaltField =
       Field::create("Salinity", "Salinity at cell centers", "",
                     "sea_water_salinity", 0.0, 100.0, FillValue, 2, DimNames);

   // Create Tracer group
   auto TracerGroup = FieldGroup::create("Tracers");

   // Add fields to tracer group
   Err1 = TracerGroup->addField("Temperature");
   TestEval("Add Temperature to tracer group", Err1, ErrRef, Err);
   Err1 = TracerGroup->addField("Salinity");
   TestEval("Add Salinity to tracer group", Err1, ErrRef, Err);

   // Also create a Restart group
   auto RestartGroup = FieldGroup::create("Restart");

   // Add fields to restart group
   Err1 = RestartGroup->addField("Temperature");
   TestEval("Add Temperature to restart group", Err1, ErrRef, Err);
   Err1 = RestartGroup->addField("Salinity");
   TestEval("Add Salinity to restart group", Err1, ErrRef, Err);

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

      std::shared_ptr<Clock> ModelClock = nullptr;

      // Call initialization to create reference IO field
      Err1 = initIOStreamTest(ModelClock);
      TestEval("Initialize IOStream test", Err1, ErrRef, Err);

      // Retrieve dimension lengths and some mesh info
      I4 NCellsSize     = Dimension::getDimLengthLocal("NCells");
      I4 NVertLevels    = Dimension::getDimLengthLocal("NVertLevels");
      Decomp *DefDecomp = Decomp::getDefault();
      I4 NCellsOwned    = DefDecomp->NCellsOwned;
      Array1DI4 CellID  = DefDecomp->CellID;

      // Create data arrays

      Array2DR8 Temp("Temp", NCellsSize, NVertLevels);
      Array2DR8 Salt("Salt", NCellsSize, NVertLevels);
      Array2DR8 Test("Test", NCellsSize, NVertLevels);

      // Attach data arrays to fields

      Err1 = Field::attachFieldData<Array2DR8>("Temperature", Temp);
      TestEval("Attach temperature data to field", Err1, ErrRef, Err);
      Err1 = Field::attachFieldData<Array2DR8>("Salinity", Salt);
      TestEval("Attach salinity data to field", Err1, ErrRef, Err);

      // Validate all streams (Mesh stream already validated in HorzMesh?)
      bool AllValidated = IOStream::validateAll();
      TestEval("IOStream Validation", AllValidated, true, Err);

      // Read restart file for initial temperature and salinity data
      Metadata ReqMetadata; // leave empty for now - no required metadata
      Err1 = IOStream::read("InitialState", *ModelClock, ReqMetadata);
      TestEval("Read restart file", Err1, ErrRef, Err);

      // Overwrite salinity array with values associated with global cell
      // ID to test proper indexing of IO
      parallelFor(
          {NCellsSize, NVertLevels}, KOKKOS_LAMBDA(int Cell, int K) {
             Salt(Cell, K) = 0.0001_Real * (CellID(Cell) + K);
             Test(Cell, K) = Salt(Cell, K);
          });

      // Create a stop alarm at 1 year for time stepping
      TimeInstant StopTime(0002, 1, 1, 0, 0, 0.0);
      Alarm StopAlarm("Stop Time", StopTime);
      Err1 = ModelClock->attachAlarm(&StopAlarm);
      TestEval("Attach stop alarm", Err1, ErrRef, Err);

      // Overwrite
      // Step forward in time and write files if it is time
      while (!StopAlarm.isRinging()) {
         ModelClock->advance();
         TimeInstant CurTime    = ModelClock->getCurrentTime();
         std::string CurTimeStr = CurTime.getString(4, 2, " ");

         Err1 = IOStream::writeAll(*ModelClock);
         if (Err1 != 0) // to prevent too much output in log
            TestEval("Write all streams " + CurTimeStr, Err1, ErrRef, Err);
      }

      // Force read the latest restart and check the results
      // A delay is needed here to make sure the restart file is completely
      // written before we read.
      std::this_thread::sleep_for(std::chrono::seconds(5));
      bool ForceRead = true;
      Err1 = IOStream::read("RestartRead", *ModelClock, ReqMetadata, ForceRead);
      TestEval("Restart force read", Err1, ErrRef, Err);

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
      IOStream::finalize(*ModelClock);
   }

   // Clean up environments
   TimeStepper::clear();
   OceanState::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Decomp::clear();
   Kokkos::finalize();
   MPI_Finalize();

   if (Err >= 256)
      Err = 255;

   // End of testing
   return Err;
}
//===--- End test driver for IO Streams ------------------------------------===/
