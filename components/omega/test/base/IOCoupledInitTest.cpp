//===-- Test driver for OMEGA coupled IO initialization ----------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA coupled IO initialization
///
/// This driver tests the IO::init overload that takes a SCORPIO iosystem id,
/// iotype and netCDF format. In a coupled run the iosystem is created by the
/// driver (CIME/shr_pio), so Omega must reuse it instead of creating its own.
/// This test creates an iosystem the way a driver would, using IO settings
/// that differ from the Omega config, and checks that IO::init adopts that
/// iosystem, its iotype and its netCDF format. It then writes and reads a file
/// with each driver iotype and netCDF format and checks the netCDF file
/// version in the file header.
///
//
//===-----------------------------------------------------------------------===/

#include "Config.h"
#include "DataTypes.h"
#include "Error.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "mpi.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Returns the netCDF file version (1=classic, 2=64-bit offset, 5=64-bit data)
// from the header of a classic-family netCDF file, or -1 if not netCDF.
static int readNetCDFVersion(const std::string &Filename) {
   std::ifstream File(Filename, std::ios::binary);
   char Magic[4] = {0, 0, 0, 0};
   File.read(Magic, 4);
   if (!File || Magic[0] != 'C' || Magic[1] != 'D' || Magic[2] != 'F')
      return -1;
   return static_cast<int>(Magic[3]);
}

//------------------------------------------------------------------------------
// Initializes IO with the given driver settings, then writes and reads back a
// distributed array using the default file format and checks that the file
// was written with the expected netCDF version.
static void testWriteRead(int DriverSysID,     // driver iosystem id
                          int DriverIOType,    // driver iotype
                          int DriverNCFormat,  // driver netCDF format flag
                          int ExpectedVersion, // expected netCDF file version
                          MPI_Comm Comm,       // communicator of the iosystem
                          I4 MyTask,           // this task's rank
                          I4 NumTasks          // number of tasks
) {

   IO::init(DriverSysID, DriverIOType, DriverNCFormat);

   const std::string Filename = "IOCoupledInitTest_" +
                                std::to_string(DriverIOType) + "_" +
                                std::to_string(DriverNCFormat) + ".nc";

   // One element per task so every task owns a distinct global index
   std::vector<int> DimLengths = {NumTasks};
   std::vector<int> GlobalIndx = {MyTask};
   int DecompID = IO::createDecomp(IO::IOTypeI4, 1, DimLengths, 1, GlobalIndx,
                                   IO::DefaultRearr);

   // Default file format, so the driver iotype and netCDF format are used
   int FileID;
   bool NewFile;
   IO::openFileWrite(FileID, Filename, NewFile, IO::IfExists::Replace);
   int DimID    = IO::defineDim(FileID, "NTasks", NumTasks);
   int DimIDs[] = {DimID};
   int VarID    = IO::defineVar(FileID, "Data", IO::IOTypeI4, 1, DimIDs);
   IO::endDefinePhase(FileID);

   I4 Fill     = -999;
   I4 WriteVal = 10 * MyTask + 1;
   IO::writeArray(&WriteVal, 1, &Fill, FileID, DecompID, VarID);
   IO::closeFile(FileID);

   int ReadFileID;
   IO::openFileRead(ReadFileID, Filename);
   I4 ReadVal = -1;
   int ReadVarID;
   Error Err =
       IO::readArray(&ReadVal, 1, "Data", ReadFileID, DecompID, ReadVarID);
   CHECK_ERROR_ABORT(Err, "IOCoupledInitTest: FAIL reading array from {}",
                     Filename);
   IO::closeFile(ReadFileID);

   if (ReadVal != WriteVal)
      ABORT_ERROR("IOCoupledInitTest: FAIL read {} but wrote {} in {}", ReadVal,
                  WriteVal, Filename);

   // Only the master task inspects the file; share the result
   MPI_Barrier(Comm);
   int Version = -1;
   if (MyTask == 0)
      Version = readNetCDFVersion(Filename);
   MPI_Bcast(&Version, 1, MPI_INT, 0, Comm);

   if (Version != ExpectedVersion)
      ABORT_ERROR("IOCoupledInitTest: FAIL {} is netCDF version {} but the "
                  "driver iotype {} and netCDF format {} should give {}",
                  Filename, Version, DriverIOType, DriverNCFormat,
                  ExpectedVersion);

   IO::destroyDecomp(DecompID);

   MPI_Barrier(Comm);
   if (MyTask == 0)
      std::filesystem::remove(Filename);
}

//------------------------------------------------------------------------------
// The test driver for a driver-owned SCORPIO iosystem.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      // Initialize the Machine Environment class and retrieve the default
      // environment and communicator
      MachEnv::init(MPI_COMM_WORLD);
      MachEnv *DefEnv  = MachEnv::getDefault();
      MPI_Comm DefComm = DefEnv->getComm();
      I4 MyTask        = DefEnv->getMyTask();
      I4 NumTasks      = DefEnv->getNumTasks();

      // Initialize the Logging system
      initLogging(DefEnv);
      LOG_INFO("----- Driver-owned IO System Unit Testing -----");

      // Needs more than one task so two IO tasks can be created
      if (NumTasks < 2)
         ABORT_ERROR("IOCoupledInitTest: FAIL test requires at least 2 tasks, "
                     "got {}",
                     NumTasks);

      // Create an iosystem as the driver would. The layout is deliberately
      // different from the Omega default of a single IO task, so adopting the
      // driver iosystem can be told apart from creating a new one with
      // config values.
      const I4 DriverNumIOTasks   = 2;
      const I4 DriverStride       = NumTasks / DriverNumIOTasks;
      const I4 DriverBaseTask     = 0;
      const I4 DriverRearr        = PIO_REARR_SUBSET;
      const I4 DriverIOType       = PIO_IOTYPE_NETCDF; // not IO::FmtDefault
      const I4 DriverNetCDFFormat = PIO_64BIT_OFFSET;  // not IO default

      int DriverSysID = -1;
      int PIOErr =
          PIOc_Init_Intracomm(DefComm, DriverNumIOTasks, DriverStride,
                              DriverBaseTask, DriverRearr, &DriverSysID);
      if (PIOErr != PIO_NOERR)
         ABORT_ERROR("IOCoupledInitTest: FAIL could not create the driver "
                     "iosystem");

      // Initialize IO with the driver-owned iosystem
      IO::init(DriverSysID, DriverIOType, DriverNetCDFFormat);

      if (IO::SysID != DriverSysID)
         ABORT_ERROR("IOCoupledInitTest: FAIL SysID is {} but the driver "
                     "supplied {}",
                     IO::SysID, DriverSysID);

      if (IO::DefaultFileFmt != IO::FmtNetCDF3)
         ABORT_ERROR("IOCoupledInitTest: FAIL default file format is {} but "
                     "the driver supplied iotype {}",
                     static_cast<int>(IO::DefaultFileFmt), DriverIOType);

      if (IO::DefaultNetCDFMode != DriverNetCDFFormat)
         ABORT_ERROR("IOCoupledInitTest: FAIL netCDF mode is {} but the "
                     "driver supplied {}",
                     IO::DefaultNetCDFMode, DriverNetCDFFormat);

      // A rearranger of 0 defers to the driver iosystem's rearranger
      if (static_cast<int>(IO::DefaultRearr) != 0)
         ABORT_ERROR("IOCoupledInitTest: FAIL default rearranger is {} but "
                     "expected 0 (use driver iosystem default)",
                     static_cast<int>(IO::DefaultRearr));

      // The iosystem in use must have the driver's IO task layout
      int NumIOTasks = 0;
      PIOErr         = PIOc_get_numiotasks(IO::SysID, &NumIOTasks);
      if (PIOErr != PIO_NOERR)
         ABORT_ERROR("IOCoupledInitTest: FAIL could not query SCORPIO IO "
                     "tasks");

      if (NumIOTasks != DriverNumIOTasks)
         ABORT_ERROR("IOCoupledInitTest: FAIL iosystem has {} IO tasks but the "
                     "driver created {}",
                     NumIOTasks, DriverNumIOTasks);

      // Decompositions created with the default rearranger must work
      std::vector<int> DimLengths = {NumTasks};
      std::vector<int> GlobalIndx = {MyTask};
      int DecompID = IO::createDecomp(IO::IOTypeI4, 1, DimLengths, 1,
                                      GlobalIndx, IO::DefaultRearr);
      IO::destroyDecomp(DecompID);

      // Files must be written with the driver's netCDF format for both the
      // serial and parallel netCDF iotypes
      const I4 IOTypes[]      = {PIO_IOTYPE_NETCDF, PIO_IOTYPE_PNETCDF};
      const I4 NCFormats[][2] = {{0, 1},                // classic
                                 {PIO_64BIT_OFFSET, 2}, // CDF2
                                 {PIO_64BIT_DATA, 5}};  // CDF5
      for (I4 IOType : IOTypes)
         for (const auto &NCFormat : NCFormats)
            testWriteRead(DriverSysID, IOType, NCFormat[0], NCFormat[1],
                          DefComm, MyTask, NumTasks);

      PIOc_finalize(DriverSysID);

      // Exit environments
      MachEnv::removeAll();

      LOG_INFO("IOCoupledInitTest: Successful completion");
   }

   LOG_INFO("----- Driver-owned IO System Unit Tests Successful -----");
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0; // if we made it here, return successfully

} // end of main
//===-----------------------------------------------------------------------===/
