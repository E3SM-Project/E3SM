#include "RungeKutta2Stepper.h"

namespace OMEGA {

RungeKutta2Stepper::RungeKutta2Stepper(const std::string &Name,
                                       Tendencies *Tend,
                                       AuxiliaryState *AuxState, HorzMesh *Mesh,
                                       Halo *MeshHalo)
    : TimeStepper(Name, TimeStepperType::RungeKutta2, Tend, AuxState, Mesh,
                  MeshHalo) {}

void RungeKutta2Stepper::doStep(OceanState *State, Real Time,
                                Real TimeStep) const {

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Tend->computeAllTendencies(State, AuxState, CurLevel, Time);
   updateStateByTend(State, NextLevel, State, CurLevel,  0.5*TimeStep);

   Tend->computeAllTendencies(State, AuxState, NextLevel, Time + 0.5*TimeStep);
   updateStateByTend(State, NextLevel, State, CurLevel, TimeStep);

   State->updateTimeLevels();
}

} // namespace OMEGA
