#include "ForwardBackwardStepper.h"

namespace OMEGA {

ForwardBackwardStepper::ForwardBackwardStepper(const std::string &Name,
                                               Tendencies *Tend,
                                               AuxiliaryState *AuxState,
                                               HorzMesh *Mesh, Halo *MeshHalo)
    : TimeStepper(Name, TimeStepperType::ForwardBackward, Tend, AuxState, Mesh,
                  MeshHalo) {}

void ForwardBackwardStepper::doStep(OceanState *State, Real Time,
                                    Real TimeStep) const {

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Tend->computeThicknessTendencies(State, AuxState, CurLevel, Time);
   updateThicknessByTend(State, NextLevel, State, CurLevel, TimeStep);

   Tend->computeVelocityTendencies(State, AuxState, NextLevel, CurLevel,
                                   Time + TimeStep);
   updateVelocityByTend(State, NextLevel, State, CurLevel, TimeStep);

   State->updateTimeLevels();
}

} // namespace OMEGA
