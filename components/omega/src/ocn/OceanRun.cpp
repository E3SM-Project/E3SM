//===-- ocn/OceanRun.cpp - Run Ocean Model ----------------------*- C++ -*-===//
//
// The ocnRun method advances the model forward from CurrTime until the
// EndAlarm rings.
//
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "Forcing.h"
#include "IOStream.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "SfcCoupling.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"

namespace OMEGA {

int ocnRun(TimeInstant &CurrTime ///< [inout] current sim time
) {

   // error code
   I4 Err = 0;

   // fetch default OceanState and TimeStepper
   OceanState *DefOceanState   = OceanState::getDefault();
   TimeStepper *DefTimeStepper = TimeStepper::getDefault();
   Forcing *DefForcing         = Forcing::getDefault();

   // get simulation time and other time info
   Clock *OmegaClock     = DefTimeStepper->getClock();
   Alarm *EndAlarm       = DefTimeStepper->getEndAlarm();
   TimeInterval TimeStep = DefTimeStepper->getTimeStep();
   TimeInstant SimTime   = OmegaClock->getCurrentTime();

   // Get Simulation metadata field for later updates
   std::shared_ptr<Field> SimInfo = Field::get(SimMeta);

   // time loop, integrate until EndAlarm or error encountered
   while (Err == 0 && !(EndAlarm->isRinging())) {

      // get step count
      const I8 IStep = DefTimeStepper->getStepCount();

      // placeholder: call needed pre-timestep compute here
      // (e.g. forcing routine)

      // do forward time step
      // first call to doStep can sometimes take very long
      // we want to time it separately and disable child timers
      // for that timer
      if (IStep == 0) {
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

      // Compute analysis fields whose alarms are ringing
      Analysis *DefAnalysis = Analysis::getDefault();
      DefAnalysis->computeAll();

      // write restart file/output, anything needed post-timestep

      IOStream::writeAll(OmegaClock);

      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               SimTime.getString(4, 4, "-"));
   }

   return Err;

} // end ocnRun

int ocnRun(TimeInstant &CurrTime, ///< [inout] current sim time
           bool WriteRestart ///< [in] write restart at end of coupling interval
) {

   // error code
   I4 Err = 0;

   // fetch default OceanState and TimeStepper
   OceanState *DefOceanState   = OceanState::getDefault();
   TimeStepper *DefTimeStepper = TimeStepper::getDefault();
   Forcing *DefForcing         = Forcing::getDefault();
   VertCoord *DefVertCoord     = VertCoord::getDefault();
   SfcCoupling *DefSfcCoupling = SfcCoupling::getDefault();
   Analysis *DefAnalysis       = Analysis::getDefault();

   // get simulation time and other time info
   Clock *OmegaClock     = DefTimeStepper->getClock();
   Alarm *CouplingAlarm  = DefSfcCoupling->getCouplingAlarm();
   TimeInterval TimeStep = DefTimeStepper->getTimeStep();
   TimeInstant SimTime   = OmegaClock->getCurrentTime();

   // Reset coupling alarm at the start of the coupling interval
   CouplingAlarm->reset(SimTime);

   DefSfcCoupling->importFromCoupler();
   DefSfcCoupling->applyImportFields(DefForcing, DefVertCoord);

   DefForcing->computeAll();

   // Recompute auxiliary variables and prepare for time stepping after
   // importing coupler data
   DefTimeStepper->prepareNextStep(DefOceanState, Tracers::getAll(0));

   // time loop, integrate until CouplingAlarm or error encountered
   while (Err == 0 && !(CouplingAlarm->isRinging())) {

      // get step count, over the whole simulation
      const I8 IStep = DefTimeStepper->getStepCount();

      // do forward time step
      // first call to doStep can sometimes take very long
      // we want to time it separately and disable child timers
      // for that timer
      if (IStep == 0) {
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

      // Update fields exported to the coupler
      DefSfcCoupling->updateExportFields(DefOceanState, Tracers::getAll(0));
      // Compute analysis fields whose alarms are ringing
      DefAnalysis->computeAll();
      // Write any IOStreams with their alarms ringing
      IOStream::writeAll(OmegaClock);

      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               SimTime.getString(4, 4, "-"));
   }

   // Minus 1 because we want to print completed coupling interval
   const I8 CoupledIntervalCount =
       DefTimeStepper->getStepCount() / DefSfcCoupling->getNAccumSteps() - 1;

   // export o2x fields to coupler at endof coupling interval
   DefSfcCoupling->exportToCoupler();

   // force write a restart at end of coupling interval if cpl tell us to
   if (WriteRestart) {
      IOStream::write("RestartWrite", OmegaClock, true);
   }

   LOG_INFO("ocnRun: ------ Completed coupling interval {} ------",
            CoupledIntervalCount);

   return Err;

} // end ocnRun
} // end namespace OMEGA
