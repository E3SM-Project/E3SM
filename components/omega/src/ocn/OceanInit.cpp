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
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"
#include "TimeStepper.h"

#include "mpi.h"

namespace OMEGA {

int ocnInit(MPI_Comm Comm,          ///< [in] ocean MPI communicator
            Calendar &OmegaCal,     ///< [out] sim calendar
            TimeInstant &StartTime, ///< [out] sim start time
            Alarm &EndAlarm         ///< [out] alarm to end simulation
) {

   // Error codes
   I4 RetErr = 0;
   I4 Err    = 0;

   // Init the default machine environment based on input MPI communicator
   MachEnv::init(Comm);
   MachEnv *DefEnv = MachEnv::getDefault();

   // Initialize Omega logging
   initLogging(DefEnv);

   // Read config file into Config object
   Config("omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_ERROR("ocnInit: Error reading config file");
      ++RetErr;
   }
   Config *OmegaConfig = Config::getOmegaConfig();

   // read and save time management options from Config
   Err = initTimeManagement(OmegaCal, StartTime, EndAlarm, OmegaConfig);
   if (Err != 0) {
      LOG_ERROR("ocnInit: Error initializing time management");
      ++RetErr;
   }

   // initialize remaining Omega modules
   Err = initOmegaModules(Comm);
   if (Err != 0) {
      LOG_ERROR("ocnInit: Error initializing Omega modules");
      ++RetErr;
   }

   return RetErr;
} // end ocnInit

// Read time management options from config
int initTimeManagement(Calendar &OmegaCal, TimeInstant &StartTime,
                       Alarm &EndAlarm, Config *OmegaConfig) {

   // error code
   I4 RetErr = 0;

   // create RunInterval and zero length interval for comparison
   TimeInterval RunInterval, ZeroInterval;

   // extract variables for time management group
   Config TimeMgmtConfig("TimeManagement");
   if (OmegaConfig->existsGroup("TimeManagement")) {
      int Err1 = OmegaConfig->get(TimeMgmtConfig);
      // check requested calendar is a valid option, return error if not found
      if (TimeMgmtConfig.existsVar("CalendarType")) {
         std::string ConfigCalStr;
         CalendarKind ConfigCalKind = CalendarUnknown;
         I4 Err1     = TimeMgmtConfig.get("CalendarType", ConfigCalStr);
         I4 ICalType = 9;
         for (I4 I = 0; I < NUM_SUPPORTED_CALENDARS; ++I) {
            if (ConfigCalStr == CalendarKindName[I]) {
               ICalType      = I;
               ConfigCalKind = (CalendarKind)(ICalType + 1);
               break;
            }
         }
         if (ICalType == 9) {
            LOG_ERROR("ocnInit: Requested Calendar type not found");
            ++RetErr;
            return RetErr;
         }
         // destroy default Calendar to keep static NumCalendars member
         // accurate, then construct requested Calendar
         OmegaCal.~Calendar();
         OmegaCal = Calendar(ConfigCalStr, ConfigCalKind);
      }

      // check for start time and set if found
      if (TimeMgmtConfig.existsVar("StartTime")) {
         std::string StartTimeStr;
         I4 Err1 = TimeMgmtConfig.get("StartTime", StartTimeStr);

         StartTime = TimeInstant(&OmegaCal, StartTimeStr);
      }

      std::string NoneStr("none");
      // set RunInterval by checking for StopTime and RunDuration in Config,
      // if both are present and inconsistent, use RunDuration
      if (TimeMgmtConfig.existsVar("StopTime")) {
         std::string StopTimeStr;
         I4 Err1 = TimeMgmtConfig.get("StopTime", StopTimeStr);
         if (StopTimeStr != NoneStr) {
            TimeInstant StopTime(&OmegaCal, StopTimeStr);
            RunInterval = StopTime - StartTime;
         }
      }
      if (TimeMgmtConfig.existsVar("RunDuration")) {
         std::string RunDurationStr;
         I4 Err1 = TimeMgmtConfig.get("RunDuration", RunDurationStr);
         if (RunDurationStr != NoneStr) {
            TimeInterval IntervalFromStr(RunDurationStr);
            if (IntervalFromStr != RunInterval) {
               LOG_WARN("ocnInit: RunDuration and StopTime are inconsistent: "
                        "using RunDuration");
               RunInterval = IntervalFromStr;
            }
         }
      }
   }

   // return error if RunInterval set to zero
   if (RunInterval == ZeroInterval) {
      LOG_ERROR("ocnInit: Simulation run duration set to zero");
      ++RetErr;
   }

   // set EndAlarm based on length of RunInterval
   TimeInstant EndTime = StartTime + RunInterval;
   EndAlarm            = Alarm("End Alarm", EndTime);

   return RetErr;
} // end initTimeManagement

// Call init routines for remaining Omega modules
int initOmegaModules(MPI_Comm Comm) {

   // error codes
   I4 RetErr = 0;
   I4 Err    = 0;

   // initialize all necessary Omega modules
   Err = IO::init(Comm);
   if (Err != 0) {
      LOG_ERROR("ocnInit: Error initializing parallel IO");
      ++RetErr;
   }

   Err = Decomp::init();
   if (Err != 0) {
      LOG_ERROR("ocnInit: Error initializing default decomposition");
      ++RetErr;
   }

   Err = Halo::init();
   if (Err != 0) {
      LOG_ERROR("ocnInit: Error initializing default halo");
      ++RetErr;
   }

   Err = HorzMesh::init();
   if (Err != 0) {
      ++RetErr;
      LOG_ERROR("ocnInit: Error initializing default mesh");
   }

   Err = AuxiliaryState::init();
   if (Err != 0) {
      ++RetErr;
      LOG_ERROR("ocnInit: Error initializing default aux state");
   }

   Err = Tendencies::init();
   if (Err != 0) {
      ++RetErr;
      LOG_ERROR("ocnInit: Error initializing default tendencies");
   }

   Err = TimeStepper::init();
   if (Err != 0) {
      ++RetErr;
      LOG_ERROR("ocnInit: Error initializing default time stepper");
   }

   Err = OceanState::init();
   if (Err != 0) {
      ++RetErr;
      LOG_ERROR("ocnInit: Error initializing default state");
   }

   return RetErr;
} // end initOmegaModules

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
