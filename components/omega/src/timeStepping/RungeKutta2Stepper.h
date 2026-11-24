#ifndef OMEGA_TSRK2_H
#define OMEGA_TSRK2_H
//===-- RungeKutta2Stepper.h - 2nd-order Runge Kutta time step --*- C++ -*-===//
//
/// \file
/// \brief Contains the class for the midpoint Runge Kutta scheme
//
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"

namespace OMEGA {

class RungeKutta2Stepper : public TimeStepper {
 public:
   // name, tendencies, auxiliary state, mesh, and halo
   /// Constructor creates an instance of a midpoint Runge Kutta stepper and
   /// fills with some time information. Data pointers are added later.
   RungeKutta2Stepper(
       const std::string &InName,              ///< [in] name of time stepper
       const TimeInterval &InTimeStep,         ///< [in] time step
       const TimeStepperStartType InStartType, ///< [in] option for starting
       const TimeInstant &InStartTime,         ///< [in] start time for full sim
       const TimeStepperStopType InStopType,   ///< [in] option for stopping
       std::optional<TimeInstant> InStopTime,  ///< [in] stop time if opt AtTime
       std::optional<TimeInterval> InDuration  ///< [in] duration if AfterDur
   );

   /// Advance the state by one step of the midpoint Runge Kutta scheme
   void doStep(OceanState *State,   ///< [inout] model state
               TimeInstant &SimTime ///< [inout] current simulation time
   ) const override;
};

} // namespace OMEGA
#endif
