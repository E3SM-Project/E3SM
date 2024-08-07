#include "RungeKutta4Stepper.h"

namespace OMEGA {

RungeKutta4Stepper::RungeKutta4Stepper(const std::string &Name,
                                       Tendencies *Tend,
                                       AuxiliaryState *AuxState, HorzMesh *Mesh,
                                       Halo *MeshHalo)
    : TimeStepper(Name, TimeStepperType::RungeKutta4, Tend, AuxState, Mesh,
                  MeshHalo) {

   auto NVertLevels = Tend->LayerThicknessTend.extent_int(1);

   ProvisState = OceanState::create("Provis", Mesh, MeshHalo, NVertLevels, 1);

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

void RungeKutta4Stepper::doStep(OceanState *State, Real Time,
                                Real TimeStep) const {

   const int CurLevel  = 0;
   const int NextLevel = 1;

   for (int Stage = 0; Stage < NStages; ++Stage) {
      const Real StageTime = Time + RKC[Stage] * TimeStep;

      if (Stage == 0) {
         Tend->computeAllTendencies(State, AuxState, CurLevel, StageTime);
         updateStateByTend(State, NextLevel, State, CurLevel,
                           RKB[Stage] * TimeStep);
      } else {
         updateStateByTend(ProvisState, CurLevel, State, CurLevel,
                           RKA[Stage] * TimeStep);

         // TODO(mwarusz) this depends on halo width actually
         if (Stage == 2) {
            ProvisState->exchangeHalo(CurLevel);
         }

         Tend->computeAllTendencies(ProvisState, AuxState, CurLevel, StageTime);
         updateStateByTend(State, NextLevel, RKB[Stage] * TimeStep);
      }
   }

   State->updateTimeLevels();
}

} // namespace OMEGA
