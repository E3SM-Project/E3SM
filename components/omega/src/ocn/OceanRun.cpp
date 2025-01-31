//===-- ocn/OceanRun.cpp - Run Ocean Model ----------------------*- C++ -*-===//
//
// The ocnRun method advances the model forward from CurrTime until the
// EndAlarm rings.
//
//===----------------------------------------------------------------------===//

#include "IOStream.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "TimeMgr.h"
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

   // Get Simulation metadata field for later updates
   std::shared_ptr<Field> SimInfo = Field::get(SimMeta);

   // time loop, integrate until EndAlarm or error encountered
   I8 IStep = 0;
   while (Err == 0 && !(EndAlarm->isRinging())) {

      // track step count
      ++IStep;

      // call forcing routines, anything needed pre-timestep

      // do forward time step
      DefTimeStepper->doStep(DefOceanState, SimTime);

      // write restart file/output, anything needed post-timestep

      Err = IOStream::writeAll(OmegaClock);
      if (Err != 0) {
         LOG_CRITICAL("Error writing streams at end of step");
         break;
      }

      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               SimTime.getString(4, 4, "-"));
   }

   return Err;

} // end ocnRun

} // end namespace OMEGA
