//===-- ocn/OceanInit.cpp - Ocean Initialization ----------------*- C++ -*-===//
//
// This file contians ocnInit and associated methods which initialize Omega.
// The ocnInit process reads the config file and uses the config options to
// initialize time management and call all the individual initialization
// routines for each module in Omega.
//
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "Forcing.h"
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
#include "SfcCoupling.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "VertMix.h"

#include "mpi.h"

namespace OMEGA {

// Convienvence converter of an int to a StartType enum, with error checking
StartType safeIntToStartType(int val) {
   switch (val) {
   case 0:
      return StartType::StartUp;
   case 1:
      return StartType::Continue;
   case 2:
      return StartType::Branch;
   default:
      ABORT_ERROR("Invalid start type value: {}", val);
   }
}

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

   // Update Halo/Host arrays with new state, auxiliary state, and tracer fields
   Err = initUpdateHaloAndHostArrays();

   return Err;
} // end ocnInit

int ocnInit1(MPI_Comm Comm,                 ///< [in] ocean MPI communicator
             const int OcnId,               ///< [in] mct comp id for ocean
             const std::string &ConfigFile, ///< [in] path to yaml config file
             const std::string &LogFile,    ///< [in] path to log file
             const StartType StartType,     ///< [in] simulation start type
             const TimeInitParams &TimeParams, ///< [in] simulation start time
             const CouplingInitParams &CouplingParams ///< [in] coupler info
) {

   I4 Err = 0; // return error code

   MachEnv *DefEnv = MachEnv::getDefault();

   // Read config file into Config object
   Config("Omega");
   Config::readAll(ConfigFile);
   Config *OmegaConfig = Config::getOmegaConfig();

   readTimingConfig();

   // initialize remaining Omega modules
   Err = initOmegaModules(Comm, TimeParams, CouplingParams);
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

   Metadata ReqMeta;
   if (StartType == StartType::StartUp) {
      // read from initial state if this is starting a new simulation
      Error IOError = IOStream::read("InitialState", ModelClock, ReqMeta);
      if (IOError.isFail()) {
         ABORT_ERROR("Errors encountered reading InitialState");
      }
   } else if (StartType == StartType::Continue ||
              StartType == StartType::Branch) {
      // read restart if starting from restart
      ReqMeta["SimulationTime"] = " ";
      Error IOError = IOStream::read("RestartRead", ModelClock, ReqMeta);
      if (IOError.isFail()) {
         ABORT_ERROR("Errors encountered reading RestartRead");
      }
      // TODO: check restat time matches coupler time.
   };

   // Advance clock one coupling interval, to be in sync with couplers clock
   if (StartType == StartType::StartUp) {
      SfcCoupling *DefCoupling = SfcCoupling::getDefault();
      while (!DefCoupling->CouplingAlarm.isRinging()) {
         ModelClock->advance();
      }
   }

   return Err;
} // end ocnInit1

// Coupled init phase 2: attach the coupler's MCT buffers and exchange the
// initial coupled state; split from ocnInit1 since these buffers don't exist
// until the coupler has sized/allocated them using Omega's decomposition
int ocnInit2(const Real *CplToOcnData, Real *OcnToCplData) {
   SfcCoupling *DefCoupling = SfcCoupling::getDefault();
   DefCoupling->attachData(CplToOcnData, OcnToCplData);

   DefCoupling->exportToCoupler();
   DefCoupling->importFromCoupler();
   DefCoupling->applyImportFields(Forcing::getDefault());

   return initUpdateHaloAndHostArrays();
} // end ocnInit2

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

   HorzMesh::init(ModelClock);
   VertCoord::init();
   Tracers::init();
   VertAdv::init();
   Forcing::init();
   AuxiliaryState::init();
   Eos::init();
   PressureGrad::init();
   VertMix::init();
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

   Analysis::init();

   return Err;

} // end initOmegaModulesImpl

int initOmegaModules(MPI_Comm Comm) {
   // Initialize the default time stepper (phase 1) that includes the
   // calendar, model clock and start/stop times and alarms with all options
   // read from the config file
   TimeStepper::init1();
   return initOmegaModulesImpl(Comm);
}

int initOmegaModules(MPI_Comm Comm, const TimeInitParams &TParams,
                     const CouplingInitParams &CParams) {
   int Err = 0;
   // Initialize time stepper (phase 1) using coupler provided time parameters
   // Calendar should have already been initalized
   TimeStepper::init1(TParams);
   Err = initOmegaModulesImpl(Comm);
   SfcCoupling::init(CParams);

   return Err;
}

int initUpdateHaloAndHostArrays() {
   // Update Halo/Host arrays with new state, auxiliary state, and tracer fields
   int Err = 0;

   OceanState *DefState = OceanState::getDefault();
   I4 CurTimeLevel      = 0;
   DefState->exchangeHalo(CurTimeLevel);

   // Enforce layer masks on state and tracer variables: fully-inactive layers
   // get FillValueReal, boundary layers get 0, active layers keep their
   // IC/restart value.
   DefState->applyLayerMasks(CurTimeLevel);

   DefState->copyToHost(CurTimeLevel);
   VertCoord::getDefault()->initSurfacePressure(Halo::getDefault());

   Forcing *DefForcing = Forcing::getDefault();
   DefForcing->exchangeHalo();

   AuxiliaryState *DefAuxState = AuxiliaryState::getDefault();
   DefAuxState->exchangeHalo();

   // Now update tracers - assume using same time level index
   Err = Tracers::exchangeHalo(CurTimeLevel);
   if (Err != 0) {
      ABORT_ERROR("Error updating tracer halo");
   }
   Tracers::copyToHost(CurTimeLevel);

   return Err;
} // end initUpdateHaloAndHostArrays

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
