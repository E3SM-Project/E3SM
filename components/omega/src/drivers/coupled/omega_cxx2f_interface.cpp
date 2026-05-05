//===-- drivers/standalone/OceanDriver.cpp - Coupled driver --*- C++ -*-===//
//
//
//===----------------------------------------------------------------------===//
#include "DataTypes.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include <cstdio>
#include <mpi.h>

extern "C" {

void omega_ocn_init(
    const MPI_Fint FComm,       // [in] MPI communicator from Fortran
    const int OcnID,            // [in] mct comp id for ocn mode
    const char *YamlConfigFile, // [in] yaml file name for ocean model
    const char *OcnLogFile,     // [in] log file name for ocean model
    const char *CalendarName,   // [in] CIME calendar name
    const int RunStartYMD,      // [in] run start date in YYYYMMDD
    const int RunStartTOD       // [in] run start time in seconds of day
) {

   // Create the C MPI_Comm from the Fortran one
   MPI_Comm Comm = MPI_Comm_f2c(FComm);

   // initialize Kokkos
   Kokkos::initialize();

   // initialize Pacer timing in coupled mode
   Pacer::initialize(Comm, Pacer::PACER_INTEGRATED);

   // Convert calendar from char to string
   std::string CalendarKindStr = CalendarName;
   std::string LogFileStr      = OcnLogFile;

   // Recall that e3sm uses the int YYYYMMDD to store a date
   OMEGA::I8 Year   = (RunStartYMD / 100) / 100;
   OMEGA::I8 Month  = (RunStartYMD / 100) % 100;
   OMEGA::I8 Day    = RunStartYMD % 100;
   OMEGA::I8 Hour   = (RunStartTOD / 60) / 60;
   OMEGA::I8 Minute = (RunStartTOD / 60) % 60;
   OMEGA::R8 Second = RunStartTOD % 60;

   // Initialize machine environment and logging before initializing
   // the calendar to ensure all output goes to the correct log file
   OMEGA::MachEnv::init(Comm);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();

   initLogging(DefEnv, LogFileStr);

   // Initialize the Calendar prior to creating the TimeInstant
   OMEGA::Calendar::init(CalendarKindStr);

   OMEGA::TimeInstant StartTime(Year, Month, Day, Hour, Minute, Second);

   // Pacer::start("Init", 0);

   OMEGA::ocnInit(Comm, OcnID, YamlConfigFile, OcnLogFile, StartTime);

   // Pacer::stop("Init", 0);

   LOG_INFO("ocnInit: Finished initializing ocean model");
   int ErrAll;
}

} // extern "C"
