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
       const std::string &InName,      ///< [in] name of time stepper
       const TimeInstant &InStartTime, ///< [in] start time for time stepping
       const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
       const TimeInterval &InTimeStep  ///< [in] time step
   );

   /// Advance the state by one step of the midpoint Runge Kutta scheme
   void doStep(OceanState *State,   ///< [inout] model state
               TimeInstant &SimTime ///< [inout] current simulation time
   ) const override;
};

} // namespace OMEGA
#endif
