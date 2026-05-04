//===-- ocn/OceanInit.cpp - Ocean Initialization ----------------*- C++ -*-===//
//
// This file contians ocnInit and associated methods which initialize Omega.
// The ocnInit process reads the config file and uses the config options to
// initialize time management and call all the individual initialization
// routines for each module in Omega.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "PGrad.h"
#include "Pacer.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"

#include "mpi.h"

namespace OMEGA {

namespace Timing {
// Flag to determine if timing info should be printed from all ranks
// Set by ocnInit. Access outside of this file is provided by
// the printAllRanks() function below
static bool PrintAllRanks = false;
} // namespace Timing

// Accessor function for the Timing::PrintAllRanks flag
bool printTimingAllRanks() { return Timing::PrintAllRanks; }

// Read timing configuration and set Pacer options
static void readTimingConfig() {
   Error Err;

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimingConfig("Timing");
   Err += OmegaConfig->get(TimingConfig);
   CHECK_ERROR_ABORT(Err, "Timing: Timing group not found in Config");

   int TimingLevel;
   Err += TimingConfig.get("Level", TimingLevel);
   CHECK_ERROR_ABORT(Err, "Timing: Level not found in TimingConfig");
   OMEGA_REQUIRE(TimingLevel >= 0, "Invalid timing level {} < 0", TimingLevel);
   Pacer::setTimingLevel(TimingLevel);

   bool AutoFence;
   Err += TimingConfig.get("AutoFence", AutoFence);
   CHECK_ERROR_ABORT(Err, "Timing: AutoFence not found in TimingConfig");
   if (AutoFence) {
      Pacer::enableAutoFence();
   }

   bool TimingBarriers;
   Err += TimingConfig.get("TimingBarriers", TimingBarriers);
   CHECK_ERROR_ABORT(Err, "Timing: TimingBarriers not found in TimingConfig");
   if (TimingBarriers) {
      Pacer::enableTimingBarriers();
   }

   Err += TimingConfig.get("PrintAllRanks", Timing::PrintAllRanks);
   CHECK_ERROR_ABORT(Err, "Timing: PrintAllRanks not found in TimingConfig");
}

int ocnInit(MPI_Comm Comm ///< [in] ocean MPI communicator
) {

   I4 Err = 0; // return error code

   // Init the default machine environment based on input MPI communicator
   MachEnv::init(Comm);
   MachEnv *DefEnv = MachEnv::getDefault();

   // Initialize Omega logging
   initLogging(DefEnv);

   // Read config file into Config object
   Config("Omega");
   Config::readAll("omega.yml");
   Config *OmegaConfig = Config::getOmegaConfig();

   readTimingConfig();

   // initialize remaining Omega modules
   Err = initOmegaModules(Comm);
   if (Err != 0)
      ABORT_ERROR("ocnInit: Error initializing Omega modules");

   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Now that all fields have been defined, validate all the streams
   // contents
   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid) {
      ABORT_ERROR("ocnInit: Error validating IO Streams");
   }

   // Initialize data from Restart or InitialState files
   std::string SimTimeStr          = " "; // create SimulationTime metadata
   std::shared_ptr<Field> SimField = Field::get(SimMeta);
   SimField->addMetadata("SimulationTime", SimTimeStr);
   Error Err1;
   Error Err2;

   // read from initial state if this is starting a new simulation
   Metadata ReqMeta; // no requested metadata for initial state
   Err1 = IOStream::read("InitialState", ModelClock, ReqMeta);

   // read restart if starting from restart
   SimTimeStr                = " ";
   ReqMeta["SimulationTime"] = SimTimeStr;
   Err2 = IOStream::read("RestartRead", ModelClock, ReqMeta);

   // One of the above two streams must be successful to initialize the
   // state and other fields used in the model
   if (Err1.isFail() and Err2.isFail()) {
      CHECK_ERROR(Err1, "Errors encountered reading InitialState");
      CHECK_ERROR(Err2, "Errors encountered reading RestartRead");
      ABORT_ERROR("Error initializing ocean variables from input streams");
   }

   // If reading from restart, reset the current time to the input time
   if (SimTimeStr != " ") {
      TimeInstant NewCurrentTime(SimTimeStr);
      ModelClock->setCurrentTime(NewCurrentTime);
   }

   // Update Halos and Host arrays with new state, auxiliary state, and tracer
   // fields

   OceanState *DefState = OceanState::getDefault();
   I4 CurTimeLevel      = 0;
   DefState->exchangeHalo(CurTimeLevel);
   DefState->copyToHost(CurTimeLevel);

   AuxiliaryState *DefAuxState = AuxiliaryState::getDefault();
   DefAuxState->exchangeHalo();

   // Now update tracers - assume using same time level index
   Err = Tracers::exchangeHalo(CurTimeLevel);
   if (Err != 0) {
      ABORT_ERROR("Error updating tracer halo after restart");
   }
   Tracers::copyToHost(CurTimeLevel);

   return Err;

} // end ocnInit

// Call init routines for remaining Omega modules
// Internal helper — all module init after TimeStepper::init1 is called.
// Called by both initOmegaModules overloads.
static int initOmegaModulesImpl(MPI_Comm Comm) {

   // error and return codes
   int Err = 0;

   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize IOStreams - this does not yet validate the contents
   // of each file, only creates streams from Config
   IOStream::init(ModelClock);

   IO::init(Comm);
   Field::init(ModelClock);
   Decomp::init();

   Err = Halo::init();
   if (Err != 0) {
      ABORT_ERROR("ocnInit: Error initializing default halo");
   }

   HorzMesh::init();
   VertCoord::init();
   Tracers::init();
   VertAdv::init();
   AuxiliaryState::init();
   Eos::init();
   PressureGrad::init();
   Tendencies::init();

   // Validate SurfaceTracerRestoring configuration
   Tendencies *DefTend = Tendencies::getDefault();
   if (DefTend->SurfaceTracerRestoring.Enabled &&
       DefTend->SurfaceTracerRestoring.NTracersToRestore == 0) {
      ABORT_ERROR("OceanInit: SurfaceTracerRestoring is enabled but "
                  "TracersToRestore is empty");
   }

   TimeStepper::init2();

   Err = OceanState::init();
   if (Err != 0) {
      ABORT_ERROR("ocnInit: Error initializing default state");
   }

   return Err;

} // end initOmegaModulesImpl

int initOmegaModules(MPI_Comm Comm) {
   // Initialize the default time stepper (phase 1) that includes the
   // calendar, model clock and start/stop times and alarms with all options
   // read from the config file
   TimeStepper::init1();
   return initOmegaModulesImpl(Comm);
}

int initOmegaModules(MPI_Comm Comm, const TimeInitParams &TParams) {
   // Initialize time stepper (phase 1) using coupler provided time parameters
   // Calendar should have already been initalized
   TimeStepper::init1(TParams);
   return initOmegaModulesImpl(Comm);
}

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
