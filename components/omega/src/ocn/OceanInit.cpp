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

   I4 Err = 0; // Error code

   // Init the default machine environment based on input MPI communicator
   MachEnv::init(Comm);
   MachEnv *DefEnv = MachEnv::getDefault();

   // Initialize Omega logging
   initLogging(DefEnv);

   // Read config file into Config object
   Config("omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error reading config file");
      return Err;
   }
   Config *OmegaConfig = Config::getOmegaConfig();

   // read and save time management options from Config
   Err = initTimeManagement(OmegaCal, StartTime, EndAlarm, OmegaConfig);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing time management");
      return Err;
   }

   // initialize remaining Omega modules
   Err = initOmegaModules(Comm);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing Omega modules");
      return Err;
   }

   return Err;
} // end ocnInit

// Read time management options from config
int initTimeManagement(Calendar &OmegaCal, TimeInstant &StartTime,
                       Alarm &EndAlarm, Config *OmegaConfig) {

   // error code
   I4 Err = 0;

   // create RunInterval and zero length interval for comparison
   TimeInterval RunInterval, ZeroInterval;

   // extract variables for time management group
   Config TimeMgmtConfig("TimeManagement");
   Err = OmegaConfig->get(TimeMgmtConfig);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: TimeManagement group not found in Config");
      return Err;
   }

   // check requested calendar is a valid option, return error if not found
   std::string ConfigCalStr;
   Err = TimeMgmtConfig.get("CalendarType", ConfigCalStr);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: CalendarType not found in TimeMgmtConfig");
      return Err;
   }
   CalendarKind ConfigCalKind = CalendarUnknown;
   I4 ICalType                = CalendarUnknown;
   for (I4 I = 0; I < NUM_SUPPORTED_CALENDARS; ++I) {
      if (ConfigCalStr == CalendarKindName[I]) {
         ICalType      = I;
         ConfigCalKind = (CalendarKind)(ICalType + 1);
         break;
      }
   }
   if (ICalType == CalendarUnknown) {
      LOG_CRITICAL("ocnInit: Requested Calendar type not found");
      Err = -1;
      return Err;
   }
   // destroy default Calendar to keep static NumCalendars member
   // accurate, then construct requested Calendar
   OmegaCal.~Calendar();
   OmegaCal = Calendar(ConfigCalStr, ConfigCalKind);

   // retrieve start time from config
   std::string StartTimeStr;
   Err = TimeMgmtConfig.get("StartTime", StartTimeStr);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: StartTime not found in TimeMgmtConfig");
      return Err;
   }
   StartTime = TimeInstant(&OmegaCal, StartTimeStr);

   std::string NoneStr("none");

   // set RunInterval by checking for StopTime and RunDuration in Config,
   // if both are present and inconsistent, use RunDuration
   std::string StopTimeStr;
   I4 Err1 = TimeMgmtConfig.get("StopTime", StopTimeStr);
   if (Err1 != 0) {
      LOG_WARN("ocnInit: StopTime not found in TimeMgmtConfig");
   } else if (StopTimeStr != NoneStr) {
      TimeInstant StopTime(&OmegaCal, StopTimeStr);
      RunInterval = StopTime - StartTime;
   }
   std::string RunDurationStr;
   I4 Err2 = TimeMgmtConfig.get("RunDuration", RunDurationStr);
   if (Err2 != 0) {
      LOG_WARN("ocnInit: RunDuration not found in TimeMgmtConfig");
   } else if (RunDurationStr != NoneStr) {
      TimeInterval IntervalFromStr(RunDurationStr);
      if (IntervalFromStr != RunInterval) {
         LOG_WARN("ocnInit: RunDuration and StopTime are inconsistent: "
                  "using RunDuration");
         RunInterval = IntervalFromStr;
      }
   }

   // return error if RunInterval set to zero
   if (RunInterval == ZeroInterval) {
      LOG_CRITICAL("ocnInit: Simulation run duration set to zero");
      Err = -1;
      return Err;
   }

   // set EndAlarm based on length of RunInterval
   TimeInstant EndTime = StartTime + RunInterval;
   EndAlarm            = Alarm("End Alarm", EndTime);

   return Err;
} // end initTimeManagement

// Call init routines for remaining Omega modules
int initOmegaModules(MPI_Comm Comm) {

   // error code
   I4 Err = 0;

   // initialize all necessary Omega modules
   Err = IO::init(Comm);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing parallel IO");
      return Err;
   }

   Err = Field::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing Fields");
      return Err;
   }

   Err = Decomp::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default decomposition");
      return Err;
   }

   Err = Halo::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default halo");
      return Err;
   }

   Err = HorzMesh::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default mesh");
      return Err;
   }

   Err = AuxiliaryState::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default aux state");
      return Err;
   }

   Err = Tendencies::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default tendencies");
      return Err;
   }

   Err = TimeStepper::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default time stepper");
      return Err;
   }

   Err = OceanState::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default state");
      return Err;
   }

   return Err;
} // end initOmegaModules

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
