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

   // q = (h,u)
     
   // R_q^{n} = RHS_q(u^{n}, h^{n}, t^{n})
   Tend->computeAllTendencies(State, AuxState, CurLevel, Time);

   // q^{n+0.5} = q^{n} + 0.5*dt*R_q^{n}
   updateStateByTend(State, NextLevel, State, CurLevel, 0.5 * TimeStep);

   // R_q^{n+0.5} = RHS_q(u^{n+0.5}, h^{n+0.5}, t^{n+0.5})
   Tend->computeAllTendencies(State, AuxState, NextLevel,
                              Time + 0.5 * TimeStep);

   // q^{n+1} = q^{n} + dt*R_q^{n+0.5}
   updateStateByTend(State, NextLevel, State, CurLevel, TimeStep);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   State->updateTimeLevels();
}

} // namespace OMEGA
