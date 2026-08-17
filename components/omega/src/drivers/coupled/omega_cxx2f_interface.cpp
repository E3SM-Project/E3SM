//===-- drivers/standalone/OceanDriver.cpp - Coupled driver --*- C++ -*-===//
//
//
//===----------------------------------------------------------------------===//
#include "DataTypes.h"
#include "Decomp.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include <mpi.h>
#include <vector>

#ifdef HAVE_MOAB
#include "MoabInterface.h"
#endif

// helper C++ functions
namespace {

// coupler provides (year, month, day) and (time of day) as ints
OMEGA::TimeInstant buildStartTime(int RunStartYMD, int RunStartTOD) {
   OMEGA::I8 Year   = (RunStartYMD / 100) / 100;
   OMEGA::I8 Month  = (RunStartYMD / 100) % 100;
   OMEGA::I8 Day    = RunStartYMD % 100;
   OMEGA::I8 Hour   = (RunStartTOD / 60) / 60;
   OMEGA::I8 Minute = (RunStartTOD / 60) % 60;
   OMEGA::R8 Second = RunStartTOD % 60;
   return OMEGA::TimeInstant(Year, Month, Day, Hour, Minute, Second);
}

std::map<std::string, int> buildFieldIndexMap(const char *FieldNames,
                                              const int *FieldIndices,
                                              const int NFields) {
   std::map<std::string, int> FieldIdx;
   for (int I = 0; I < NFields; ++I) {
      // null-terminator stops string, so only start index is needed
      std::string FieldName(FieldNames + I * 32);
      // convert from 1-based index (Fortran) to 0-based index (C++)
      FieldIdx[FieldName] = FieldIndices[I] - 1;
   }
   return FieldIdx;
}

#ifdef HAVE_MOAB
// Coupling buffers for the MOAB path: attached once (by address) to
// SfcCoupling in omega_ocn_init2, refilled from/drained to MOAB tag storage
// around each omega_ocn_run call.
int MoabNCouplerImports = 0;
int MoabNCouplerExports = 0;
std::vector<OMEGA::Real> MoabCplToOcn;
std::vector<OMEGA::Real> MoabOcnToCpl;

// The ocean's own MOAB application id, returned by moabInit. The Fortran
// cap must copy this into seq_comm_mct's mpoid module variable, since that
// is what the coupler-side migration/mapping code reads to find the
// ocean's mesh.
int MoabPid = -1;
#endif

} // namespace

extern "C" {

void omega_ocn_init1(
    const MPI_Fint FComm,          // [in] MPI communicator from Fortran
    const int OcnID,               // [in] mct comp id for ocn mode
    const char *YamlConfigFile,    // [in] yaml file name for ocean model
    const char *OcnLogFile,        // [in] log file name for ocean model
    const int StartType,           // [in] 0=startup, 1=continue, 2=branch
    const char *CalendarName,      // [in] CIME calendar name
    const int RunStartYMD,         // [in] run start date in YYYYMMDD
    const int RunStartTOD,         // [in] run start time in seconds of day
    const int CouplingTimeStep,    // [in] coupling time step in seconds
    const int NCouplerImports,     // [in] number of coupler import fields
    const int NCouplerExports,     // [in] number of coupler export fields
    const int NOmegaImports,       // [in] number of omega import fields
    const int NOmegaExports,       // [in] number of omega export fields
    const char *ImportFieldNames,  // [in] array of import field names
    const char *ExportFieldNames,  // [in] array of export field names
    const int *ImportFieldIndices, // [in] array of import field indices
    const int *ExportFieldIndices, // [in] array of export field indices
    const char *Cpl2OcnFieldNames, // [in] full CIME x2o field list (MOAB)
    const char *Ocn2CplFieldNames  // [in] full CIME o2x field list (MOAB)
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

   // Initialize machine environment and logging before initializing
   // the calendar to ensure all output goes to the correct log file
   OMEGA::MachEnv::init(Comm);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();

   initLogging(DefEnv, LogFileStr);

   // Init Calendar prior to creating the TimeInstant/TimeInterval objects
   OMEGA::Calendar::init(CalendarKindStr);

   OMEGA::TimeInstant StartTime = buildStartTime(RunStartYMD, RunStartTOD);
   OMEGA::TimeInterval CouplingInterval(CouplingTimeStep,
                                        OMEGA::TimeUnits::Seconds);

   std::map<std::string, int> ImportIdxMap =
       buildFieldIndexMap(ImportFieldNames, ImportFieldIndices, NOmegaImports);
   std::map<std::string, int> ExportIdxMap =
       buildFieldIndexMap(ExportFieldNames, ExportFieldIndices, NOmegaExports);

   Pacer::start("Init1", 0);

   OMEGA::StartType StartTypeEnum = OMEGA::safeIntToStartType(StartType);
   OMEGA::TimeInitParams TimeParams{StartTime, std::nullopt};
   OMEGA::CouplingInitParams CouplingParams{
       NCouplerImports, NCouplerExports,  ImportIdxMap,
       ExportIdxMap,    CouplingInterval,
#ifdef HAVE_MOAB
       OMEGA::CouplingLayout::MOAB};
#else
       OMEGA::CouplingLayout::MCT};
#endif

   OMEGA::ocnInit1(Comm, OcnID, YamlConfigFile, OcnLogFile, StartTypeEnum,
                   TimeParams, CouplingParams);

   Pacer::stop("Init1", 0);

   LOG_INFO("ocnInit: Finished initializing ocean model");
   int ErrAll;

#ifdef HAVE_MOAB
   // Decomp/HorzMesh exist by now (built inside OMEGA::ocnInit1 above), so
   // the MOAB mesh can be constructed from them.
   MoabNCouplerImports = NCouplerImports;
   MoabNCouplerExports = NCouplerExports;
   MoabPid              = OMEGA::moabInit(Comm, OcnID);
   OMEGA::moabDefineTagStorage(MoabPid, Cpl2OcnFieldNames, Ocn2CplFieldNames);
#endif
}

#ifdef HAVE_MOAB
int omega_get_moab_pid() { return MoabPid; }
#endif

void omega_ocn_init2(const double *cpl_to_ocn_data, double *ocn_to_cpl_data) {
   Pacer::start("Init2", 0);
#ifdef HAVE_MOAB
   // The MCT attribute-vector pointers above are meaningless under MOAB;
   // attach Omega's own buffers instead, filled from MOAB tag storage.
   // These are refilled/drained in place (same memory address) by
   // omega_ocn_run on every subsequent coupling interval.
   const int NCellsOwned = static_cast<int>(OMEGA::Decomp::getDefault()->NCellsOwned);
   MoabCplToOcn.assign(static_cast<size_t>(MoabNCouplerImports) * NCellsOwned,
                       0);
   MoabOcnToCpl.assign(static_cast<size_t>(MoabNCouplerExports) * NCellsOwned,
                       0);
   OMEGA::moabImportTagStorage(MoabCplToOcn.data(),
                               static_cast<int>(MoabCplToOcn.size()));
   OMEGA::ocnInit2(MoabCplToOcn.data(), MoabOcnToCpl.data());
#else
   OMEGA::ocnInit2(cpl_to_ocn_data, ocn_to_cpl_data);
#endif
   Pacer::stop("Init2", 0);
}

int omega_ocn_run(bool WriteRestart) {

   int ErrRun;

   OMEGA::TimeStepper *DefStepper = OMEGA::TimeStepper::getDefault();
   OMEGA::Clock *ModelClock       = DefStepper->getClock();
   OMEGA::TimeInstant CurrTime    = ModelClock->getCurrentTime();

   Pacer::start("Run", 0);
#ifdef HAVE_MOAB
   OMEGA::moabImportTagStorage(MoabCplToOcn.data(),
                               static_cast<int>(MoabCplToOcn.size()));
#endif
   ErrRun = OMEGA::ocnRun(CurrTime, WriteRestart);
#ifdef HAVE_MOAB
   OMEGA::moabExportTagStorage(MoabOcnToCpl.data(),
                               static_cast<int>(MoabOcnToCpl.size()));
#endif
   Pacer::stop("Run", 0);

   return ErrRun;
}

int omega_ocn_finalize() {

   int ErrFinalize;

   OMEGA::TimeStepper *DefStepper = OMEGA::TimeStepper::getDefault();
   OMEGA::Clock *ModelClock       = DefStepper->getClock();
   OMEGA::TimeInstant CurrTime    = ModelClock->getCurrentTime();

   Pacer::start("Finalize", 0);
   OMEGA::ocnFinalize(CurrTime);
   if (ErrFinalize != 0) {
      LOG_ERROR("Error finalizing OMEGA");
   } else {
      LOG_INFO("OMEGA successfully completed");
   }
   Pacer::stop("Finalize", 0);

   // no Pacer::print or Pacer::finalize in coupled mode; cpl will handle it

   // finalize Kokkos
   Kokkos::finalize();

   return ErrFinalize;
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

   OMEGA::R8 SphereRadius = HMesh->SphereRadius;

   // Sync device array back to host mirror
   OMEGA::deepCopy(HMesh->AreaCellH, HMesh->AreaCell);

   for (int Cell = 0; Cell < HMesh->NCellsOwned; ++Cell) {
      AreaCell[Cell] = static_cast<double>(HMesh->AreaCellH[Cell] /
                                           (SphereRadius * SphereRadius));
   }
}

} // extern "C"
