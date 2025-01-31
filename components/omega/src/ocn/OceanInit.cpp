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
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"

#include "mpi.h"

namespace OMEGA {

int ocnInit(MPI_Comm Comm ///< [in] ocean MPI communicator
) {

   I4 Err = 0; // Error code

   // Init the default machine environment based on input MPI communicator
   MachEnv::init(Comm);
   MachEnv *DefEnv = MachEnv::getDefault();

   // Initialize Omega logging
   initLogging(DefEnv);

   // Read config file into Config object
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error reading config file");
      return Err;
   }
   Config *OmegaConfig = Config::getOmegaConfig();

   // initialize remaining Omega modules
   Err = initOmegaModules(Comm);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing Omega modules");
      return Err;
   }

   return Err;
} // end ocnInit

// Call init routines for remaining Omega modules
int initOmegaModules(MPI_Comm Comm) {

   // error code
   I4 Err = 0;

   // Initialize the default time stepper (phase 1) that includes the
   // calendar, model clock and start/stop times and alarms
   Err = TimeStepper::init1();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error phase 1 initializing default time stepper");
      return Err;
   }

   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize IOStreams - this does not yet validate the contents
   // of each file, only creates streams from Config
   Err = IOStream::init(ModelClock);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing IOStreams");
      return Err;
   }

   Err = IO::init(Comm);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing parallel IO");
      return Err;
   }

   Err = Field::init(ModelClock);
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

   // Create the vertical dimension - this will eventually move to
   // a vertical mesh later
   Config *OmegaConfig = Config::getOmegaConfig();
   Config DimConfig("Dimension");
   Err = OmegaConfig->get(DimConfig);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Dimension group not found in Config");
      return Err;
   }
   I4 NVertLevels;
   Err = DimConfig.get("NVertLevels", NVertLevels);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: NVertLevels not found in Dimension Config");
      return Err;
   }
   auto VertDim = OMEGA::Dimension::create("NVertLevels", NVertLevels);

   Err = Tracers::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing tracers infrastructure");
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

   Err = TimeStepper::init2();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error phase 2 initializing default time stepper");
      return Err;
   }

   Err = OceanState::init();
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: Error initializing default state");
      return Err;
   }

   // Now that all fields have been defined, validate all the streams
   // contents
   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid) {
      LOG_CRITICAL("ocnInit: Error validating IO Streams");
      return Err;
   }

   // Initialize data from Restart or InitialState files
   std::string SimTimeStr          = " "; // create SimulationTime metadata
   std::shared_ptr<Field> SimField = Field::get(SimMeta);
   SimField->addMetadata("SimulationTime", SimTimeStr);
   int Err1 = IOStream::Success;
   int Err2 = IOStream::Success;

   // read from initial state if this is starting a new simulation
   Metadata ReqMeta; // no requested metadata for initial state
   Err1 = IOStream::read("InitialState", ModelClock, ReqMeta);

   // read restart if starting from restart
   SimTimeStr                = " ";
   ReqMeta["SimulationTime"] = SimTimeStr;
   Err2 = IOStream::read("RestartRead", ModelClock, ReqMeta);

   // One of the above two streams must be successful to initialize the
   // state and other fields used in the model
   if (Err1 != IOStream::Success and Err2 != IOStream::Success) {
      LOG_CRITICAL("Error initializing ocean variables from input streams");
      return Err1 + Err2;
   }

   // If reading from restart, reset the current time to the input time
   if (SimTimeStr != " ") {
      TimeInstant NewCurrentTime(SimTimeStr);
      Err = ModelClock->setCurrentTime(NewCurrentTime);
      if (Err != 0) {
         LOG_CRITICAL("Error resetting the simulation time from restart");
         return Err;
      }
   }

   // Update Halos and Host arrays with new state and tracer fields

   OceanState *DefState = OceanState::getDefault();
   I4 CurTimeLevel      = 0;
   DefState->exchangeHalo(CurTimeLevel);
   DefState->copyToHost(CurTimeLevel);

   // Now update tracers - assume using same time level index
   Err = Tracers::exchangeHalo(CurTimeLevel);
   if (Err != 0) {
      LOG_CRITICAL("Error updating tracer halo after restart");
      return Err;
   }
   Err = Tracers::copyToHost(CurTimeLevel);
   if (Err != 0) {
      LOG_CRITICAL("Error updating tracer device arrays after restart");
      return Err;
   }

   return Err;

} // end initOmegaModules

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
