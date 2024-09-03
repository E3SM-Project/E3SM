//===-- drivers/standalone/OceanDriver.cpp - Standalone driver --*- C++ -*-===//
//
//
//===----------------------------------------------------------------------===//

#include "OceanDriver.h"
#include "DataTypes.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"

#include "mpi.h"

#include <iostream>

int main(int argc, char **argv) {

   int ErrAll;
   int Err1;
   int Err2;

   MPI_Init(&argc, &argv); // initialize MPI
   Kokkos::initialize();   // initialize Kokkos

   // Time management objects
   OMEGA::Calendar OmegaCal;
   OMEGA::TimeInstant CurrTime;
   OMEGA::Alarm EndAlarm;

   Err1 = OMEGA::ocnInit(MPI_COMM_WORLD, OmegaCal, CurrTime, EndAlarm);
   if (Err1 != 0)
      LOG_ERROR("Error initializing OMEGA");

   while (Err1 == 0 && !(EndAlarm.isRinging())) {

      Err1 = OMEGA::ocnRun(CurrTime, EndAlarm);

      if (Err1 != 0)
         LOG_ERROR("Error advancing Omega run interval");
   }

   Err2 = OMEGA::ocnFinalize(CurrTime);
   if (Err2 != 0)
      LOG_ERROR("Error finalizing OMEGA");

   ErrAll = abs(Err1) + abs(Err2);
   if (ErrAll == 0) {
      LOG_INFO("OMEGA successfully completed");
   } else {
      LOG_ERROR("OMEGA terminating due to error");
   }

   Kokkos::finalize();
   MPI_Finalize();

   if (ErrAll >= 256)
      ErrAll = 255;
   return ErrAll;
}

//===----------------------------------------------------------------------===//
