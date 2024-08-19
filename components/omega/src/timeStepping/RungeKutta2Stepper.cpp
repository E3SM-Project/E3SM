#include "RungeKutta2Stepper.h"

namespace OMEGA {

// Constructor. Construct a midpoint Runge Kutta stepper from
// name, tendencies, auxiliary state, mesh, and halo
RungeKutta2Stepper::RungeKutta2Stepper(const std::string &Name,
                                       Tendencies *Tend,
                                       AuxiliaryState *AuxState, HorzMesh *Mesh,
                                       Halo *MeshHalo)
    : TimeStepper(Name, TimeStepperType::RungeKutta2, Tend, AuxState, Mesh,
                  MeshHalo) {}

// Advance the state by one step of the midpoint Runge Kutta scheme
void RungeKutta2Stepper::doStep(OceanState *State, Real Time,
                                Real TimeStep) const {

   const int CurLevel  = 0;
   const int NextLevel = 1;

   // $k_1=\Delta f(\phi_n, t_n)$
   Tend->computeAllTendencies(State, AuxState, CurLevel, Time);

   // $\phi_{n+{1\over2}}=\phi_n+{k_1 \over2}$
   updateStateByTend(State, NextLevel, State, CurLevel, 0.5 * TimeStep);

   // $k_2=\Delta f(\phi_n+{1\over2}, t_n+{\Delta t\over2})$
   Tend->computeAllTendencies(State, AuxState, NextLevel,
                              Time + 0.5 * TimeStep);

   // $\phi_{n+1}=\phi_n+k_2$
   updateStateByTend(State, NextLevel, State, CurLevel, TimeStep);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   State->updateTimeLevels();
}

} // namespace OMEGA
