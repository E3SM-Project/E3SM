#ifndef OMEGA_TSFB_H
#define OMEGA_TSFB_H
//===-- ForwardBackwardStepper.h - forward-backward time step --*- C++ -*--===//
//
/// \file
/// \brief Contains the class for forward-backward time stepping scheme
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"

namespace OMEGA {

class ForwardBackwardStepper : public TimeStepper {
 public:
   /// Constructor creates an instance of a forward-backward stepper and
   /// fills with some time information. Data pointers are added later.
   ForwardBackwardStepper(
       const std::string &InName,      ///< [in] name of time stepper
       const TimeInstant &InStartTime, ///< [in] start time for time stepping
       const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
       const TimeInterval &InTimeStep  ///< [in] time step
   );

   /// Advance the state by one step of the forward-backward scheme
   void doStep(OceanState *State,   ///< [inout] model state
               TimeInstant &SimTime ///< [inout] current simulation time
   ) const override;
};

} // namespace OMEGA
#endif
