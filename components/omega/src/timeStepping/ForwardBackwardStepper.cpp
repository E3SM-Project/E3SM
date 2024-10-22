#include "ForwardBackwardStepper.h"

namespace OMEGA {

// Constructor. Construct a forward-backward stepper from
// name, tendencies, auxiliary state, mesh, and halo
ForwardBackwardStepper::ForwardBackwardStepper(const std::string &Name,
                                               Tendencies *Tend,
                                               AuxiliaryState *AuxState,
                                               HorzMesh *Mesh, Halo *MeshHalo)
    : TimeStepper(Name, TimeStepperType::ForwardBackward, 2, Tend, AuxState,
                  Mesh, MeshHalo) {}

// Advance the state by one step of the forward-backward scheme
void ForwardBackwardStepper::doStep(OceanState *State, TimeInstant Time) const {

   int Err = 0;

   const int CurLevel  = 0;
   const int NextLevel = 1;
   // TODO: resolve time level indexing
   const int TrCurLevel  = -1;
   const int TrNextLevel = 0;

   Array3DReal CurTracerArray, NextTracerArray;
   Err = Tracers::getAll(CurTracerArray, TrCurLevel);
   Err = Tracers::getAll(NextTracerArray, TrNextLevel);

   // R_h^{n} = RHS_h(u^{n}, h^{n}, t^{n})
   Tend->computeThicknessTendencies(State, AuxState, CurLevel, CurLevel, Time);

   // h^{n+1} = h^{n} + R_h^{n}
   updateThicknessByTend(State, NextLevel, State, CurLevel, TimeStep);

   // R_phi^{n} = RHS_phi(u^{n}, h^{n}, phi^{n}, t^{n})
   Tend->computeTracerTendencies(State, AuxState, CurTracerArray, CurLevel,
                                 CurLevel, Time);

   // phi^{n+1} = (phi^{n} * h^{n} + R_phi^{n}) / h^{n+1}
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel, State,
                       CurLevel, TimeStep);

   // R_u^{n+1} = RHS_u(u^{n}, h^{n+1}, t^{n+1})
   Tend->computeVelocityTendencies(State, AuxState, NextLevel, CurLevel,
                                   Time + TimeStep);

   // u^{n+1} = u^{n} + R_u^{n+1}
   updateVelocityByTend(State, NextLevel, State, CurLevel, TimeStep);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   State->updateTimeLevels();
   Tracers::updateTimeLevels();
}

} // namespace OMEGA
