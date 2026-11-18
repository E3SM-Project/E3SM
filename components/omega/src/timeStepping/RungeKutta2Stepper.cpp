//===-- RungeKutta2Stepper.cpp - 2nd-order Runge Kutta methods --*- C++ -*-===//
//
// Contains methods for the midpoint Runge-Kutta time stepping scheme
//
//===----------------------------------------------------------------------===//

#include "RungeKutta2Stepper.h"
#include "Pacer.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructor creates an instance of a midpoint Runge Kutta stepper and
// fills with some time information. Data pointers are added later.
// Mostly just passes info to the base constructor.
RungeKutta2Stepper::RungeKutta2Stepper(
    const std::string &InName,              // [in] name of time stepper
    const TimeInterval &InTimeStep,         // [in] time step
    const TimeStepperStartType InStartType, // [in] option for starting
    const TimeInstant &InStartTime,         // [in] start time for full sim
    const TimeStepperStopType InStopType,   // [in] option for stopping
    std::optional<TimeInstant> InStopTime,  // [in] stop time if option AtTime
    std::optional<TimeInterval> InDuration  // [in] duration if AfterDuration
    )
    : TimeStepper(InName, TimeStepperType::RungeKutta2, 2, InTimeStep,
                  InStartType, InStartTime, InStopType, InStopTime,
                  InDuration) {
   FirstProjFactor = 0.5_Real;
}

//------------------------------------------------------------------------------
// Advance the state by one step of the midpoint Runge Kutta scheme
void RungeKutta2Stepper::doStep(OceanState *State,   // model state
                                TimeInstant &SimTime // current simulation time
) const {

   OMEGA_REQUIRE(PrepareFirstStepCalled,
                 "prepareFirstStep was not called before doStep");

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Array3DReal CurTracerArray  = Tracers::getAll(CurLevel);
   Array3DReal NextTracerArray = Tracers::getAll(NextLevel);

   VertMix *VMix = VertMix::getInstance();

   resetFrazilFluxes();

   // q = (h,u,phi)
   // R_q^{n} = RHS_q(u^{n}, h^{n}, phi^{n}, t^{n})
   Tend->computeAllTendenciesOnly(State, AuxState, CurTracerArray, CurLevel,
                                  CurLevel, CurLevel, SimTime, 0.0_Real);

   // q^{n+0.5} = q^{n} + 0.5*dt*R_q^{n}
   updateStateByTend(State, NextLevel, State, CurLevel, 0.5 * TimeStep);
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel, State,
                       CurLevel, 0.5 * TimeStep);

   State->exchangeHalo(NextLevel);
   MeshHalo->exchangeFullArrayHalo(NextTracerArray, OnCell);

   prescribeState(State, NextLevel, State, CurLevel, SimTime + 0.5 * TimeStep);

   // R_q^{n+0.5} = RHS_q(u^{n+0.5}, h^{n+0.5}, phi^{n+0.5}, t^{n+0.5})
   Tend->computeAllTendencies(State, AuxState, NextTracerArray, NextLevel,
                              NextLevel, NextLevel, SimTime + 0.5 * TimeStep,
                              0.5 * TimeStep, 1.0_Real);

   // q^{n+1} = q^{n} + dt*R_q^{n+0.5}
   updateStateByTend(State, NextLevel, State, CurLevel, TimeStep);
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel, State,
                       CurLevel, TimeStep);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   const MPI_Comm Comm = MeshHalo->getComm();
   Pacer::timingBarrier("RK2:haloExch", 3, Comm);
   Pacer::start("RK2:haloExch", 3);
   State->updateTimeLevels();
   Tracers::updateTimeLevels();
   Pacer::stop("RK2:haloExch", 3);

   applyImplicitVerticalMixing(State, CurLevel, CurLevel, CurLevel, "RK2");

   // Advance the clock and update the simulation time
   StepClock->advance();
   SimTime = StepClock->getCurrentTime();
   ++StepCount;

   prepareNextStep(State, Tracers::getAll(CurLevel));
}

} // namespace OMEGA
