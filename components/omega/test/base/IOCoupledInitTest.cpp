//===-- Test driver for OMEGA coupled IO initialization ----------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA coupled IO initialization
///
/// This driver tests the IO::init overload that takes a SCORPIO iosystem id
/// and iotype. In a coupled run the iosystem is created by the driver
/// (CIME/shr_pio), so Omega must reuse it instead of creating its own. This
/// test creates an iosystem the way a driver would, using IO settings that
/// differ from the Omega config, and checks that IO::init adopts that
/// iosystem and its iotype.
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

#include <vector>

using namespace OMEGA;

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
      const I4 DriverNumIOTasks = 2;
      const I4 DriverStride     = NumTasks / DriverNumIOTasks;
      const I4 DriverBaseTask   = 0;
      const I4 DriverRearr      = PIO_REARR_SUBSET;
      const I4 DriverIOType     = PIO_IOTYPE_NETCDF; // not IO::FmtDefault

      int DriverSysID = -1;
      int PIOErr =
          PIOc_Init_Intracomm(DefComm, DriverNumIOTasks, DriverStride,
                              DriverBaseTask, DriverRearr, &DriverSysID);
      if (PIOErr != PIO_NOERR)
         ABORT_ERROR("IOCoupledInitTest: FAIL could not create the driver "
                     "iosystem");

      // Initialize IO with the driver-owned iosystem
      IO::init(DriverSysID, DriverIOType);

      if (IO::SysID != DriverSysID)
         ABORT_ERROR("IOCoupledInitTest: FAIL SysID is {} but the driver "
                     "supplied {}",
                     IO::SysID, DriverSysID);

      if (IO::DefaultFileFmt != IO::FmtNetCDF3)
         ABORT_ERROR("IOCoupledInitTest: FAIL default file format is {} but "
                     "the driver supplied iotype {}",
                     static_cast<int>(IO::DefaultFileFmt), DriverIOType);

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
