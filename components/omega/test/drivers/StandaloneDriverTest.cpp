//===-- Test driver for Omega standalone driver  -----------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for Omega standalone driver
///
/// This driver runs the standalone model to confirm that the initialization
/// phase builds all the necessary modules, the run phase successfully
/// integrates the model a few time steps, and the finalization phase
/// successfully cleans up all objects and exits without error.
///
//
//===-----------------------------------------------------------------------===/

#include "OceanDriver.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"

#include "mpi.h"

//------------------------------------------------------------------------------
// The test driver for the standalone driver
//
int main(int argc, char *argv[]) {

   OMEGA::I4 ErrAll;
   OMEGA::I4 Err1;
   OMEGA::I4 Err2;

   MPI_Init(&argc, &argv); // initialize MPI
   Kokkos::initialize();   // initialize Kokkos

   // Time management objects
   OMEGA::Calendar OmegaCal;
   OMEGA::TimeInstant CurrTime;
   OMEGA::Alarm EndAlarm;

   Err1 = OMEGA::ocnInit(MPI_COMM_WORLD, OmegaCal, CurrTime, EndAlarm);
   if (Err1 == 0) {
      LOG_INFO("DriverTest: Omega initialize PASS");
   } else {
      LOG_INFO("DriverTest: Omega initialize FAIL");
   }

   if (Err1 == 0) {
      Err1 = OMEGA::ocnRun(CurrTime, EndAlarm);
   }
   if (Err1 == 0) {
      LOG_INFO("DriverTest: Omega model run PASS");
   } else {
      LOG_INFO("DriverTest: Omega model run FAIL");
   }

   Err2 = OMEGA::ocnFinalize(CurrTime);
   if (Err2 == 0) {
      LOG_INFO("DriverTest: Omega finalize PASS");
   } else {
      LOG_INFO("DriverTest: Omega finalize FAIL");
   }

   ErrAll = abs(Err1) + abs(Err2);
   if (ErrAll == 0) {
      LOG_INFO("DriverTest: Successful completion");
   }

   Kokkos::finalize();
   MPI_Finalize();

   if (ErrAll >= 256)
      ErrAll = 255;
   return ErrAll;

} // end of main
//===-----------------------------------------------------------------------===/
