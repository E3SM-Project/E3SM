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

   // EndAlarm must be set before calling ocnRun
   OMEGA_REQUIRE(DefTimeStepper->hasEndAlarm(), "ocnRun: no EndAlarm");

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
      // first call to doStep can sometimes take very long
      // we want to time it separately and disable child timers
      // for that timer
      if (IStep == 1) {
         Pacer::start("Stepper:firstDoStep", 1);
         Pacer::disableTiming();
         DefTimeStepper->doStep(DefOceanState, SimTime);
         Pacer::enableTiming();
         Pacer::stop("Stepper:firstDoStep", 1);
      } else {
         Pacer::start("Stepper:doStep", 1);
         DefTimeStepper->doStep(DefOceanState, SimTime);
         Pacer::stop("Stepper:doStep", 1);
      }

      // write restart file/output, anything needed post-timestep

      IOStream::writeAll(OmegaClock);

      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               SimTime.getString(4, 4, "-"));
   }

   return Err;

} // end ocnRun

} // end namespace OMEGA
