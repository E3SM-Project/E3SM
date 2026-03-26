//===-- ForwardBackwardStepper.cpp - forward-backward methods --*- C++ -*--===//
//
// Contains methods for the Forward-Backward time stepper
//
//===----------------------------------------------------------------------===//

#include "ForwardBackwardStepper.h"
#include "Pacer.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructor creates an instance of a forward-backward stepper and
// fills with some time information. Data pointers are added later.
// Mostly passes relevant info to the base constructor.
ForwardBackwardStepper::ForwardBackwardStepper(
    const std::string &InName,      ///< [in] name of time stepper
    const TimeInstant &InStartTime, ///< [in] start time for time stepping
    const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
    const TimeInterval &InTimeStep  ///< [in] time step
    )
    : TimeStepper(InName, TimeStepperType::ForwardBackward, 2, InStartTime,
                  InStopTime, InTimeStep) {}

//------------------------------------------------------------------------------
// Advance the state by one step of the forward-backward scheme
void ForwardBackwardStepper::doStep(
    OceanState *State,   // input model state
    TimeInstant &SimTime // current simulation time
) const {

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Array3DReal CurTracerArray  = Tracers::getAll(CurLevel);
   Array3DReal NextTracerArray = Tracers::getAll(NextLevel);

   if (State == nullptr)
      LOG_CRITICAL("Invalid State");
   if (AuxState == nullptr)
      LOG_CRITICAL("Invalid AuxState");

   // R_h^{n} = RHS_h(u^{n}, h^{n}, t^{n})
   Tend->computeThicknessTendencies(State, AuxState, CurLevel, CurLevel,
                                    SimTime);

   // h^{n+1} = h^{n} + R_h^{n}
   updateThicknessByTend(State, NextLevel, State, CurLevel, TimeStep);

   // R_phi^{n} = RHS_phi(u^{n}, h^{n}, phi^{n}, t^{n})
   Tend->computeTracerTendencies(State, AuxState, CurTracerArray, CurLevel,
                                 CurLevel, SimTime);

   // phi^{n+1} = (phi^{n} * h^{n} + R_phi^{n}) / h^{n+1}
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel, State,
                       CurLevel, TimeStep);

   // R_u^{n+1} = RHS_u(u^{n}, h^{n+1}, t^{n+1})
   Tend->computeVelocityTendencies(State, AuxState, NextLevel, CurLevel,
                                   NextLevel, SimTime + TimeStep);

   // u^{n+1} = u^{n} + R_u^{n+1}
   updateVelocityByTend(State, NextLevel, State, CurLevel, TimeStep);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   const MPI_Comm Comm = MeshHalo->getComm();
   Pacer::timingBarrier("ForwardBackward:haloExchBarrier", 3, Comm);
   Pacer::start("ForwardBackward:haloExch", 3);
   State->updateTimeLevels();
   Tracers::updateTimeLevels();
   Pacer::stop("ForwardBackward:haloExch", 3);

   // Advance the clock and update the simulation time
   StepClock->advance();
   SimTime = StepClock->getCurrentTime();
}

} // namespace OMEGA
