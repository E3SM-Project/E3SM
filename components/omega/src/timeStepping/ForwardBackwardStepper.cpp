#include "ForwardBackwardStepper.h"

namespace OMEGA {

// Constructor. Construct a forward-backward stepper from
// name, tendencies, auxiliary state, mesh, and halo
ForwardBackwardStepper::ForwardBackwardStepper(const std::string &Name,
                                               Tendencies *Tend,
                                               AuxiliaryState *AuxState,
                                               HorzMesh *Mesh, Halo *MeshHalo)
    : TimeStepper(Name, TimeStepperType::ForwardBackward, Tend, AuxState, Mesh,
                  MeshHalo) {}

// Advance the state by one step of the forward-backward scheme
void ForwardBackwardStepper::doStep(OceanState *State, Real Time,
                                    Real TimeStep) const {

   const int CurLevel  = 0;
   const int NextLevel = 1;

   // R_h^{n} = RHS_h(u^{n}, h^{n}, t^{n})
   Tend->computeThicknessTendencies(State, AuxState, CurLevel, Time);

   // h^{n+1} = h^{n} + R_h^{n}
   updateThicknessByTend(State, NextLevel, State, CurLevel, TimeStep);

   // R_u^{n+1} = RHS_u(u^{n}, h^{n+1}, t^{n+1})
   Tend->computeVelocityTendencies(State, AuxState, NextLevel, CurLevel,
                                   Time + TimeStep);
   // u^{n+1} = u^{n} + R_u^{n+1}
   updateVelocityByTend(State, NextLevel, State, CurLevel, TimeStep);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   State->updateTimeLevels();
}

} // namespace OMEGA
