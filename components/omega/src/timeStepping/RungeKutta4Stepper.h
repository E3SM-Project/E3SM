#ifndef OMEGA_TSRK4_H
#define OMEGA_TSRK4_H
//===-- RungeKutta4Stepper.h - 4th-order Runge Kutta time step --*- C++ -*-===//
//
/// \file
/// \brief Contains the class for the classic fourth order Runge Kutta scheme
//
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"

namespace OMEGA {

class RungeKutta4Stepper : public TimeStepper {
 public:
   /// Constructor creates an instance of a fourth order Runge Kutta stepper and
   /// fills with some time information. Data pointers are added later.
   RungeKutta4Stepper(
       const std::string &InName,      ///< [in] name of time stepper
       const TimeInstant &InStartTime, ///< [in] start time for time stepping
       const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
       const TimeInterval &InTimeStep  ///< [in] time step
   );

   /// Advance the state by one step of the fourth-order Runge Kutta scheme
   void doStep(OceanState *State,   ///< [inout] model state
               TimeInstant &SimTime ///< [inout] current simulation time
   ) const override;

 protected:
   /// Performs some additional initialization for provisional fields
   void finalizeInit() override;

 private:
   /// Number of stages
   static constexpr int NStages = 4;

   /// Provisional state
   OceanState *ProvisState;

   /// Provisional tracers
   Array3DReal ProvisTracers;

   /// Runge-Kutta coefficients
   Real RKA[NStages];
   Real RKB[NStages];
   Real RKC[NStages];
};

} // namespace OMEGA
#endif
