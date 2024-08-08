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
   // TODO(mwarusz): this copy could be avoided with a more flexible interface
   deepCopy(State->NormalVelocity[NextLevel], State->NormalVelocity[CurLevel]);

   Tend->computeVelocityTendencies(State, AuxState, NextLevel, Time + TimeStep);
   updateVelocityByTend(State, NextLevel, TimeStep);

   State->updateTimeLevels();
}

} // namespace OMEGA
