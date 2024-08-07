#ifndef OMEGA_TSRK4_H
#define OMEGA_TSRK4_H

#include "TimeStepper.h"

namespace OMEGA {

class RungeKutta4Stepper : public TimeStepper {
 public:
   RungeKutta4Stepper(const std::string &Name, Tendencies *Tend,
                      AuxiliaryState *AuxState, HorzMesh *Mesh, Halo *MeshHalo);

   void doStep(OceanState *State, Real Time, Real TimeStep) const override;

 private:
   static constexpr int NStages = 4;
   OceanState *ProvisState;
   Real RKA[NStages];
   Real RKB[NStages];
   Real RKC[NStages];
};

} // namespace OMEGA
#endif
