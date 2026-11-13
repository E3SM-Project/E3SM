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
#include "IO.h"
#include "SfcCoupling.h"
#include "TimeMgr.h"
#include "TimeStepper.h"

#include "mpi.h"

namespace OMEGA {

/// Should timing info be printed from all ranks
bool printTimingAllRanks();

/// Read the config file and call all the inidividual initialization routines
/// for each Omega module
int ocnInit(MPI_Comm Comm);

/// Coupled init phase 1: everything up through mesh/decomp/state, before the
/// coupler-owned MCT buffers exist for attachData
int ocnInit1(
    MPI_Comm Comm,                            ///< [in] ocean MPI communicator
    const int OcnId,                          ///< [in] mct comp id for ocean
    const std::string &ConfigFile,            ///< [in] path to yaml config file
    const std::string &LogFile,               ///< [in] path to log file
    const TimeStepperStartType StartType,     ///< [in] simulation start type
    const TimeInstant &StartTime,             ///< [in] simulation start time
    const CouplingInitParams &CouplingParams, ///< [in] coupling parameters
    const int IOSysID,       ///< [in] driver-owned IO system id
    const int IOType,        ///< [in] driver-owned IO type
    const int IONetCDFFormat ///< [in] driver-owned netCDF file format
);

/// Coupled init phase 2: runs once the coupler has allocated its MCT buffers;
/// attaches them and exchanges the initial coupled state
void ocnInit2(const Real *CplToOcnData, ///< [in] coupler import data pointer
              Real *OcnToCplData        ///< [out] coupler export data pointer
);

/// Advance the model from starting from CurrTime until EndAlarm rings
int ocnRun(TimeInstant &CurrTime);

/// Advance the model from CurrTime until the CouplingAlarm rings
int ocnRun(TimeInstant &CurrTime, bool WriteRestart);

/// Clean up all Omega objects
int ocnFinalize(const TimeInstant &CurrTime);

/// Initialize Omega modules needed to run ocean model
int initOmegaModules(MPI_Comm Comm);

/// Initialize Omega modules with coupler-provided time parameters
int initOmegaModules(
    MPI_Comm Comm,                     ///< [in] MPI communicator for ocn
    TimeStepperStartType StartType,    ///< [in] option for starting this leg
    const TimeInstant &StartTime,      ///< [in] start time for full simulation
    const CouplingInitParams &CParams, ///< [in] struct with coupling params
    int IOSysID,                       ///< [in] driver-owned IO system id
    int IOType,                        ///< [in] driver-owned IO type
    int IONetCDFFormat                 ///< [in] driver-owned netCDF file format
);

/// Update Halo/Host arrays with new state, auxiliary state, and tracer fields
int initUpdateHaloAndHostArrays();

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
