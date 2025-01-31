//===-- RungeKutta4Stepper.cpp - 4th-order Runge Kutta methods --*- C++ -*-===//
//
// Methods for the fourth-order Runge-Kutta time stepping scheme
//
//===----------------------------------------------------------------------===//

#include "RungeKutta4Stepper.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructor creates an instance of a fourth order Runge Kutta stepper and
// fills with some time information. Data pointers are added later.
// Uses the base constructor and adds some coefficients.
RungeKutta4Stepper::RungeKutta4Stepper(
    const std::string &InName,      ///< [in] name of time stepper
    const TimeInstant &InStartTime, ///< [in] start time for time stepping
    const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
    const TimeInterval &InTimeStep  ///< [in] time step
    )
    : TimeStepper(InName, TimeStepperType::RungeKutta4, 2, InStartTime,
                  InStopTime, InTimeStep) {

   RKA[0] = 0;
   RKA[1] = 1. / 2;
   RKA[2] = 1. / 2;
   RKA[3] = 1;

   RKB[0] = 1. / 6;
   RKB[1] = 1. / 3;
   RKB[2] = 1. / 3;
   RKB[3] = 1. / 6;

   RKC[0] = 0;
   RKC[1] = 1. / 2;
   RKC[2] = 1. / 2;
   RKC[3] = 1;
}

//------------------------------------------------------------------------------
// Performs some additional initialization for provisional fields
void RungeKutta4Stepper::finalizeInit() {

   if (!Tend)
      LOG_CRITICAL("Tendency not initialized");
   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!MeshHalo)
      LOG_CRITICAL("Invalid MeshHalo");

   int NVertLevels = Tend->LayerThicknessTend.extent_int(1);
   int NTracers    = Tracers::getNumTracers();
   int NCellsSize  = Mesh->NCellsSize;
   int NTimeLevels = 1; // for provisional tracer

   ProvisState = OceanState::create("Provis" + Name, Mesh, MeshHalo,
                                    NVertLevels, NTimeLevels);
   if (!ProvisState)
      LOG_CRITICAL("Error creating Provis state");

   ProvisTracers =
       Array3DReal("ProvisTracers", NTracers, NCellsSize, NVertLevels);
}

//------------------------------------------------------------------------------
// Advance the state by one step of the fourth-order Runge Kutta scheme
void RungeKutta4Stepper::doStep(OceanState *State,   // model state
                                TimeInstant &SimTime // current simulation time
) const {

   int Err = 0;

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Array3DReal NextTracerArray, CurTracerArray;
   Err = Tracers::getAll(CurTracerArray, CurLevel);
   Err = Tracers::getAll(NextTracerArray, NextLevel);

   for (int Stage = 0; Stage < NStages; ++Stage) {
      const TimeInstant StageTime = SimTime + RKC[Stage] * TimeStep;
      // first stage does:
      // R^{(0)} = RHS(q^{n}, t^{n})
      // q^{n+1} = q^{n} + dt * RKB[0] * R^{(0)}
      if (Stage == 0) {
         weightTracers(NextTracerArray, CurTracerArray, State, CurLevel);
         Tend->computeAllTendencies(State, AuxState, CurTracerArray, CurLevel,
                                    CurLevel, StageTime);
         updateStateByTend(State, NextLevel, State, CurLevel,
                           RKB[Stage] * TimeStep);
         accumulateTracersUpdate(NextTracerArray, RKB[Stage] * TimeStep);
      } else {
         // every other stage does:
         // q^{provis} = q^{n} + RKA[stage] * dt * R^{(s-1)}
         // R^{(s)} = RHS(q^{provis}, t^{n} + RKC[stage] * dt)
         // q^{n+1} += RKB[stage] * dt * R^{(s)}
         updateStateByTend(ProvisState, CurLevel, State, CurLevel,
                           RKA[Stage] * TimeStep);
         updateTracersByTend(ProvisTracers, CurTracerArray, ProvisState,
                             CurLevel, State, CurLevel, RKA[Stage] * TimeStep);

         // TODO(mwarusz) this depends on halo width actually
         if (Stage == 2) {
            ProvisState->exchangeHalo(CurLevel);
            MeshHalo->exchangeFullArrayHalo(ProvisTracers, OnCell);
         }

         Tend->computeAllTendencies(ProvisState, AuxState, ProvisTracers,
                                    CurLevel, CurLevel, StageTime);
         updateStateByTend(State, NextLevel, State, NextLevel,
                           RKB[Stage] * TimeStep);
         accumulateTracersUpdate(NextTracerArray, RKB[Stage] * TimeStep);
      }
   }

   finalizeTracersUpdate(NextTracerArray, State, NextLevel);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   State->updateTimeLevels();
   Tracers::updateTimeLevels();

   // Advance the clock and update the simulation time
   Err     = StepClock->advance();
   SimTime = StepClock->getCurrentTime();
}

} // namespace OMEGA
