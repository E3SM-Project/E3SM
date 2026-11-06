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
#include "Frazil.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "KPPMix.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "PGrad.h"
#include "Pacer.h"
#include "SfcCoupling.h"
#include "SubmesoEddies.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "VertMix.h"

#include "mpi.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Timing initialization routines

namespace Timing {
// Flag to determine if timing info should be printed from all ranks
// Set by ocnInit. Access outside of this file is provided by
// the printAllRanks() function below
static bool PrintAllRanks = false;
} // namespace Timing

// Accessor function for the Timing::PrintAllRanks flag
bool printTimingAllRanks() { return Timing::PrintAllRanks; }

// Read timing configuration and set Pacer options
static void readTimingConfig(Config *OmegaConfig) {
   Error Err;

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

// Records whether the coupled init read the state from a restart file. Set in
// ocnInit1, where the start type is known, and used in ocnInit2, which is where
// the state is far enough along for initStateForTimeStepper
static bool CoupledReadRestart = false;

//------------------------------------------------------------------------------
// Perform any time-stepper specific state initialization that can only be done
// once the input has been read and the halos have been exchanged. For the
// split-explicit stepper this separates the input velocity into its barotropic
// and baroclinic parts. Must be called after initUpdateHaloAndHostArrays.
static void initStateForTimeStepper(
    bool ReadRestart ///< [in] true if restart input initialized the state
) {
   // Both split-explicit variants need their velocity split established
   // before the first step.
   TimeStepper *DefStepper = TimeStepper::getDefault();
   OceanState *DefState    = OceanState::getDefault();
   DefStepper->initializeStateFromInput(DefState, ReadRestart);
}

//------------------------------------------------------------------------------
// Ocean initialization - standalone case
int ocnInit(MPI_Comm Comm ///< [in] ocean MPI communicator
) {

   I4 Err = 0; // return error code

   // Init the default machine environment based on input MPI communicator
   MachEnv::init(Comm);
   MachEnv *DefEnv = MachEnv::getDefault();
   OMEGA_REQUIRE(DefEnv, "Null default MachEnv pointer in ocnInit");

   // Initialize Omega logging
   initLogging(DefEnv);

   // Read config file into Config object
   Config("Omega");
   Config::readAll("omega.yml");
   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in ocnInit");

   readTimingConfig(OmegaConfig);

   // initialize remaining Omega modules
   Err = initOmegaModules(Comm);
   if (Err != 0)
      ABORT_ERROR("ocnInit: Error initializing Omega modules");

   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Now that all fields have been defined, validate all the streams
   // contents
   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid)
      ABORT_ERROR("ocnInit: Error validating IO Streams");

   if (Tendencies::getDefault()->PenetratingShortwave.Enabled) {
      Metadata ShortwaveExtinctionMeta;
      Error ShortwaveExtinctionError = IOStream::read(
          "ShortwaveExtinctionIn", ModelClock, ShortwaveExtinctionMeta);
      if (ShortwaveExtinctionError.isFail()) {
         ABORT_ERROR("Errors encountered reading ShortwaveExtinctionIn");
      }
   }

   // Initialize data from Restart or InitialState files
   std::shared_ptr<Field> SimField = Field::get(SimMeta);
   std::string SimTimeStr          = " ";
   SimField->addMetadata("SimulationTime", SimTimeStr);
   Error Err1;

   Metadata ReqMeta; // empty requested metadata from file
   TimeStepperStartType StartType = DefStepper->getStartType();
   bool ReadRestart               = false;
   // ContinueRestart is True only for true temporal continuations
   // (Continue/Branch), which resume at the restart time.
   bool ContinueRestart = false;

   // Read from either initial state stream or restart stream based
   // on the start option
   switch (StartType) {

   // Starting from scratch using an initial state
   case TimeStepperStartType::StartUp:
      Err1 = IOStream::read("InitialState", ModelClock, ReqMeta);
      CHECK_ERROR_ABORT(Err1, "Error reading InitialState file");
      break;

   // Continue simulation from a restart file and reset current time
   // to the restart time read from restart metadata
   case TimeStepperStartType::Branch:
   case TimeStepperStartType::Continue: {
      ReqMeta["SimulationTime"] = SimTimeStr; // request current sim time
      Err1            = IOStream::read("RestartRead", ModelClock, ReqMeta);
      ReadRestart     = true;
      ContinueRestart = true;

      // Reset the current time to the input time from restart file and
      // update the end alarm and stop time.
      SimTimeStr = std::any_cast<std::string>(ReqMeta["SimulationTime"]);
      if (SimTimeStr == " ")
         ABORT_ERROR("Error reading current time from restart file");
      TimeInstant NewCurrentTime(SimTimeStr);
      ModelClock->setCurrentTime(NewCurrentTime);
      DefStepper->resetEndAlarm();
   } break;

   // Start a simulation from a previous restart file but keep the
   // simulation StartTime rather than the restart time. Because the clock is
   // re-dated to StartTime, analysis accumulators must begin a fresh period, so
   // ContinueRestart is intentionally left false (no applyRestartState below).
   case TimeStepperStartType::Hybrid:
      Err1 = IOStream::read("RestartRead", ModelClock, ReqMeta);
      CHECK_ERROR_ABORT(Err1, "Error reading restart file for hybrid run");
      ReadRestart = true;
      break;

   // If the type is coupler, we have reached here in error. For coupled
   // simulations the input coupled start type is translated into one of the
   // above, so we only reached here if it's coupled but in standalone mode
   case TimeStepperStartType::Coupled:
      ABORT_ERROR("Attempted use of Coupled StartType for Standalone run");

   default:
      ABORT_ERROR("Unknown StartType in OcnInit");

   } // end switch StartType

   // Update Halo/Host arrays with new state, auxiliary state, and tracer fields
   Err = initUpdateHaloAndHostArrays();

   // Restore the sample count and period flag for each restart-persisting
   // analysis operator from the checkpoint Fields just read. (The partial sum
   // is restored directly by the restart read.) Only done for true temporal
   // continuations (Continue/Branch); Hybrid keeps the config StartTime and so
   // begins a fresh averaging period, with the first sample overwriting any
   // partial sum left in the output array by the restart read.
   if (ContinueRestart) {
      Analysis::getDefault()->applyRestartState();
   }

   // Finish any time-stepper specific state initialization
   initStateForTimeStepper(ReadRestart);

   return Err;
} // end ocnInit

//------------------------------------------------------------------------------
// Ocean initialization - coupling case
int ocnInit1(MPI_Comm Comm,                 ///< [in] ocean MPI communicator
             const int OcnId,               ///< [in] mct comp id for ocean
             const std::string &ConfigFile, ///< [in] path to yaml config file
             const std::string &LogFile,    ///< [in] path to log file
             const TimeStepperStartType StartType, ///< [in] sim start type
             const TimeInstant &StartTime, ///< [in] simulation start time
             const CouplingInitParams &CouplingParams, ///< [in] coupler info
             const int IOSysID,       ///< [in] driver-owned IO system id
             const int IOType,        ///< [in] driver-owned IO type
             const int IONetCDFFormat ///< [in] driver-owned netCDF format
) {

   I4 Err = 0; // return error code

   MachEnv *DefEnv = MachEnv::getDefault();
   OMEGA_REQUIRE(DefEnv, "Null default MachEnv pointer in ocnInit1");

   // Read config file into Config object
   // Note that most info comes from the coupler in this case, so only timing
   // configuration is read from the actual config file
   Config("Omega");
   Config::readAll(ConfigFile);
   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in ocnInit1");

   readTimingConfig(OmegaConfig);

   // initialize remaining Omega modules
   Err = initOmegaModules(Comm, StartType, StartTime, CouplingParams, IOSysID,
                          IOType, IONetCDFFormat);
   if (Err != 0)
      ABORT_ERROR("ocnInit: Error initializing Omega modules");

   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Now that all fields are defined, validate all stream contents
   bool StreamsValid = IOStream::validateAll();

   if (!StreamsValid) {
      ABORT_ERROR("ocnInit: Error validating IO Streams");
   }

   // Initialize data from Restart or InitialState files
   std::string SimTimeStr          = " "; // create SimulationTime metadata
   std::shared_ptr<Field> SimField = Field::get(SimMeta);
   SimField->addMetadata("SimulationTime", SimTimeStr);
   Error Err1;

   Metadata ReqMeta; // empty requested metadata from file

   Metadata ShortwaveExtinctionMeta;
   if (Tendencies::getDefault()->PenetratingShortwave.Enabled) {
      Metadata ShortwaveExtinctionMeta;
      Error ShortwaveExtinctionError = IOStream::read(
          "ShortwaveExtinctionIn", ModelClock, ShortwaveExtinctionMeta);
      if (ShortwaveExtinctionError.isFail()) {
         ABORT_ERROR("Errors encountered reading ShortwaveExtinctionIn");
      }
   }

   // Read from either initial state stream or restart stream based
   // on the start option
   switch (StartType) {

   // Starting from scratch using an initial state
   case TimeStepperStartType::StartUp:
      Err1 = IOStream::read("InitialState", ModelClock, ReqMeta);
      CHECK_ERROR_ABORT(Err1, "Error reading InitialState file");
      CoupledReadRestart = false;
      break;

   // Continue simulation from a restart file and reset current time
   // to the restart time read from restart metadata
   case TimeStepperStartType::Branch:
   case TimeStepperStartType::Continue: {
      ReqMeta["SimulationTime"] = SimTimeStr; // request current sim time
      Err1               = IOStream::read("RestartRead", ModelClock, ReqMeta);
      CoupledReadRestart = true;

      // Reset the current time to the input time from restart file and
      // update the end alarm and stop time.
      SimTimeStr = std::any_cast<std::string>(ReqMeta["SimulationTime"]);
      if (SimTimeStr == " ")
         ABORT_ERROR("Error reading current time from restart file");
      TimeInstant NewCurrentTime(SimTimeStr);
      ModelClock->setCurrentTime(NewCurrentTime);
      DefStepper->resetEndAlarm();
   } break;

   // For coupled simulations, these two options are translated elsewhere
   // into one of the above options, so if we have reached here, it is an error
   case TimeStepperStartType::Hybrid:
   case TimeStepperStartType::Coupled:
      ABORT_ERROR("Incorrect translation of StartType by coupled OcnInit");

   default:
      ABORT_ERROR("Unknown StartType in OcnInit");

   } // end switch StartType

   // Advance clock one coupling interval, to be in sync with couplers clock
   if (StartType == TimeStepperStartType::StartUp) {
      SfcCoupling *DefCoupling = SfcCoupling::getDefault();
      while (!DefCoupling->getCouplingAlarm()->isRinging()) {
         ModelClock->advance();
      }
   }

   Err = initUpdateHaloAndHostArrays();

   // Restore the sample count and period flag for each restart-persisting
   // analysis operator from the checkpoint Fields just read. (The partial sum
   // is restored directly by the restart read.)
   if (CoupledReadRestart) {
      Analysis::getDefault()->applyRestartState();
   }

   // Finish any time-stepper specific state initialization
   initStateForTimeStepper(CoupledReadRestart);

   return Err;
} // end ocnInit1

//------------------------------------------------------------------------------
// Coupled init phase 2: attach the coupler's MCT buffers and exchange the
// initial coupled state; split from ocnInit1 since these buffers don't exist
// until the coupler has sized/allocated them using Omega's decomposition
void ocnInit2(const Real *CplToOcnData, Real *OcnToCplData) {

   SfcCoupling *DefSfcCoupling = SfcCoupling::getDefault();
   OceanState *DefOceanState   = OceanState::getDefault();

   DefSfcCoupling->attachData(CplToOcnData, OcnToCplData);

   // Populate export fields with initial SST, SSS, velocities, and SSH
   DefSfcCoupling->updateExportFields(DefOceanState, Tracers::getAll(0));

   DefSfcCoupling->exportToCoupler();

   // No importFromCoupler/applyImportFields here. The coupler has not filled
   // x2o yet, so importing would overwrite the initial-state SurfacePressure
   // with uninitialized memory. ocnRun calls both at the start of every
   // coupling interval, including the first, so nothing is missed.
} // end ocnInit2

//------------------------------------------------------------------------------
// Call init routines for remaining Omega modules
// Internal helper — all module init after TimeStepper::init1 is called.
// Called by both initOmegaModules overloads.
static int initOmegaModulesImpl() {

   // error and return codes
   int Err = 0;

   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize IOStreams - this does not yet validate the contents
   // of each file, only creates streams from Config
   IOStream::init(ModelClock);

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
   Frazil::init();
   PressureGrad::init();
   VertMix::init();
   KPPMix::init();
   SubmesoEddies::init();
   Tendencies::init();

   // Validate SurfaceTracerRestoring configuration
   Tendencies *DefTend = Tendencies::getDefault();
   if (DefTend->SurfaceTracerRestoring.Enabled &&
       DefTend->SurfaceTracerRestoring.NTracersToRestore == 0) {
      ABORT_ERROR("OceanInit: SurfaceTracerRestoring is enabled but "
                  "TracersToRestore is empty");
   }

   // Add fields to time stepper
   TimeStepper::init2();

   Err = OceanState::init();
   if (Err != 0) {
      ABORT_ERROR("ocnInit: Error initializing default state");
   }

   Analysis::init();

   return Err;

} // end initOmegaModulesImpl

//------------------------------------------------------------------------------
// Initializes time stepper and other modules using Config inputs (standalone)
int initOmegaModules(MPI_Comm Comm) {
   // Initialize the default time stepper (phase 1) that includes the
   // calendar, model clock and start/stop times and alarms with all options
   // read from the config file
   TimeStepper::init1();
   IO::init(Comm);
   return initOmegaModulesImpl();
}

//------------------------------------------------------------------------------
// Initializes time stepper and other modules using coupler provided inputs
int initOmegaModules(MPI_Comm Comm, TimeStepperStartType StartType,
                     const TimeInstant &StartTime,
                     const CouplingInitParams &CParams, int IOSysID, int IOType,
                     int IONetCDFFormat) {
   int Err = 0;
   // Initialize time stepper (phase 1) using coupler provided time parameters
   // Calendar should have already been initalized
   TimeStepper::init1(StartType, StartTime);
   IO::init(IOSysID, IOType, IONetCDFFormat);
   Err = initOmegaModulesImpl();
   SfcCoupling::init(CParams);

   return Err;
}

//------------------------------------------------------------------------------
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
