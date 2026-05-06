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

int omega_get_ncells_local() {
   OMEGA::Decomp *OcnDecomp = OMEGA::Decomp::getDefault();

   return static_cast<int>(OcnDecomp->NCellsOwned);
}

int omega_get_ncells_global() {
   OMEGA::Decomp *OcnDecomp = OMEGA::Decomp::getDefault();

   return static_cast<int>(OcnDecomp->NCellsGlobal);
}

void omega_get_index_to_cell_id(int *CellID) {
   OMEGA::Decomp *OcnDecomp = OMEGA::Decomp::getDefault();

   // Sync device array back to host mirror
   OMEGA::deepCopy(OcnDecomp->CellIDH, OcnDecomp->CellID);

   for (int Cell = 0; Cell < OcnDecomp->NCellsOwned; ++Cell) {
      CellID[Cell] = static_cast<int>(OcnDecomp->CellIDH[Cell]);
   }
}

void omega_get_lonlat_cell(double *LonCell, double *LatCell) {
   OMEGA::HorzMesh *HMesh = OMEGA::HorzMesh::getDefault();

   for (int Cell = 0; Cell < HMesh->NCellsOwned; ++Cell) {
      LonCell[Cell] =
          static_cast<double>(HMesh->LonCellH[Cell] * OMEGA::Rad2Deg);
      LatCell[Cell] =
          static_cast<double>(HMesh->LatCellH[Cell] * OMEGA::Rad2Deg);
   }
}

void omega_get_area_cell(double *AreaCell) {
   OMEGA::HorzMesh *HMesh = OMEGA::HorzMesh::getDefault();

   // Sync device array back to host mirror
   OMEGA::deepCopy(HMesh->AreaCellH, HMesh->AreaCell);

   for (int Cell = 0; Cell < HMesh->NCellsOwned; ++Cell) {
      // TODO: Use the HMesh->SphereRadius attribute once PR#382 is merged
      AreaCell[Cell] = static_cast<double>(HMesh->AreaCellH[Cell] /
                                           (OMEGA::REarth * OMEGA::REarth));
   }
}
} // extern "C"
