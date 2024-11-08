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

   return Err;
} // end initOmegaModules

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
