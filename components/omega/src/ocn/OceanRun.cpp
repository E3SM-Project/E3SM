//===-- ocn/OceanRun.cpp - Run Ocean Model ----------------------*- C++ -*-===//
//
// The ocnRun method advances the model forward from CurrTime until the
// EndAlarm rings.
//
//===----------------------------------------------------------------------===//

#include "OceanDriver.h"
#include "OceanState.h"
#include "TimeStepper.h"

namespace OMEGA {

int ocnRun(TimeInstant &CurrTime ///< [inout] current sim time
) {

   // error code
   I4 Err = 0;

   // fetch default OceanState and TimeStepper
   OceanState *DefOceanState   = OceanState::getDefault();
   TimeStepper *DefTimeStepper = TimeStepper::getDefault();

   // get simulation time and other time info
   Clock *OmegaClock     = DefTimeStepper->getClock();
   Alarm *EndAlarm       = DefTimeStepper->getEndAlarm();
   TimeInterval TimeStep = DefTimeStepper->getTimeStep();
   TimeInstant SimTime   = OmegaClock->getCurrentTime();

   I8 IStep = 0;

   // time loop, integrate until EndAlarm or error encountered
   while (Err == 0 && !(EndAlarm->isRinging())) {

      // track step count
      ++IStep;

      // call forcing routines, anything needed pre-timestep

      // do forward time step
      DefTimeStepper->doStep(DefOceanState, SimTime);

      // write restart file/output, anything needed post-timestep

      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               SimTime.getString(4, 4, "-"));
   }

   return Err;

} // end ocnRun

} // end namespace OMEGA
