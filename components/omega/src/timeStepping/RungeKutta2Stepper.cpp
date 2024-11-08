//===-- RungeKutta2Stepper.cpp - 2nd-order Runge Kutta methods --*- C++ -*-===//
//
// Contains methods for the midpoint Runge-Kutta time stepping scheme
//
//===----------------------------------------------------------------------===//

#include "RungeKutta2Stepper.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructor creates an instance of a midpoint Runge Kutta stepper and
// fills with some time information. Data pointers are added later.
// Mostly just passes info to the base constructor.
RungeKutta2Stepper::RungeKutta2Stepper(
    const std::string &InName,      ///< [in] name of time stepper
    const TimeInstant &InStartTime, ///< [in] start time for time stepping
    const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
    const TimeInterval &InTimeStep  ///< [in] time step
    )
    : TimeStepper(InName, TimeStepperType::RungeKutta2, 2, InStartTime,
                  InStopTime, InTimeStep) {}

//------------------------------------------------------------------------------
// Advance the state by one step of the midpoint Runge Kutta scheme
void RungeKutta2Stepper::doStep(OceanState *State,   // model state
                                TimeInstant &SimTime // current simulation time
) const {

   int Err = 0;

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Array3DReal CurTracerArray, NextTracerArray;
   Err = Tracers::getAll(CurTracerArray, CurLevel);
   Err = Tracers::getAll(NextTracerArray, NextLevel);

   // q = (h,u,phi)
   // R_q^{n} = RHS_q(u^{n}, h^{n}, phi^{n}, t^{n})
   Tend->computeAllTendencies(State, AuxState, CurTracerArray, CurLevel,
                              CurLevel, SimTime);

   // q^{n+0.5} = q^{n} + 0.5*dt*R_q^{n}
   updateStateByTend(State, NextLevel, State, CurLevel, 0.5 * TimeStep);
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel, State,
                       CurLevel, 0.5 * TimeStep);

   // R_q^{n+0.5} = RHS_q(u^{n+0.5}, h^{n+0.5}, phi^{n+0.5}, t^{n+0.5})
   Tend->computeAllTendencies(State, AuxState, NextTracerArray, NextLevel,
                              NextLevel, SimTime + 0.5 * TimeStep);

   // q^{n+1} = q^{n} + dt*R_q^{n+0.5}
   updateStateByTend(State, NextLevel, State, CurLevel, TimeStep);
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel, State,
                       CurLevel, TimeStep);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   State->updateTimeLevels();
   Tracers::updateTimeLevels();

   // Advance the clock and update the simulation time
   Err     = StepClock->advance();
   SimTime = StepClock->getCurrentTime();
}

} // namespace OMEGA
