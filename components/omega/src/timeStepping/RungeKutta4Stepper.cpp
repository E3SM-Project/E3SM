#include "RungeKutta4Stepper.h"

namespace OMEGA {

// Constructor. Construct a fourth order Runge Kutta stepper from
// name, tendencies, auxiliary state, mesh, and halo
RungeKutta4Stepper::RungeKutta4Stepper(const std::string &Name,
                                       Tendencies *Tend,
                                       AuxiliaryState *AuxState, HorzMesh *Mesh,
                                       Halo *MeshHalo)
    : TimeStepper(Name, TimeStepperType::RungeKutta4, 2, Tend, AuxState, Mesh,
                  MeshHalo) {

   if (Tend)
      finalizeInit();

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

void RungeKutta4Stepper::finalizeInit() {

   auto NVertLevels = Tend->LayerThicknessTend.extent_int(1);
   auto NTracers    = Tracers::getNumTracers();

   ProvisState = OceanState::create("Provis", Mesh, MeshHalo, NVertLevels, 1);
   ProvisTracers =
       Array3DReal("ProvisTracers", NTracers, Mesh->NCellsSize, NVertLevels);
}

// Advance the state by one step of the fourth-order Runge Kutta scheme
void RungeKutta4Stepper::doStep(OceanState *State, TimeInstant Time) const {

   int Err = 0;

   const int CurLevel  = 0;
   const int NextLevel = 1;
   // TODO: resolve time level indexing
   const int TrCurLevel  = -1;
   const int TrNextLevel = 0;

   Array3DReal NextTracerArray, CurTracerArray;
   Err = Tracers::getAll(CurTracerArray, TrCurLevel);
   Err = Tracers::getAll(NextTracerArray, TrNextLevel);

   for (int Stage = 0; Stage < NStages; ++Stage) {
      const TimeInstant StageTime = Time + RKC[Stage] * TimeStep;
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
            auto ProvisTracersH = createHostMirrorCopy(ProvisTracers);
            MeshHalo->exchangeFullArrayHalo(ProvisTracersH, OnCell);
            deepCopy(ProvisTracers, ProvisTracersH);
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
}

} // namespace OMEGA
