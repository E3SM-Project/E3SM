#ifndef OMEGA_DRIVER_H
#define OMEGA_DRIVER_H
//===-- ocn/OceanDriver.h ---------------------------------------*- C++ -*-===//
//
/// \file
/// \brief Defines ocean driver methods
///
/// This Header defines methods to drive Omega. These methods are designed to
/// run Omega as either a standalone ocean model or as a component of E3SM.
/// This process is divided into three phases: init, run, and finalize.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "TimeMgr.h"

#include "mpi.h"

namespace OMEGA {

/// Read the config file and call all the inidividual initialization routines
/// for each Omega module
int ocnInit(MPI_Comm Comm, Calendar &OmegaCal, TimeInstant &StartTime,
            Alarm &EndAlarm);

/// Advance the model from starting from CurrTime until EndAlarm rings
int ocnRun(TimeInstant &CurrTime, Alarm &EndAlarm);

/// Clean up all Omega objects
int ocnFinalize(const TimeInstant &CurrTime);

/// Extract time management options from Config
int initTimeManagement(Calendar &OmegaCal, TimeInstant &StartTime,
                       Alarm &EndAlarm, Config *OmegaConfig);

/// Initialize Omega modules needed to run ocean model
int initOmegaModules(MPI_Comm Comm);

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
