#ifndef OMEGA_TSRK2_H
#define OMEGA_TSRK2_H

#include "TimeStepper.h"

namespace OMEGA {

class RungeKutta2Stepper : public TimeStepper {
 public:
   RungeKutta2Stepper(const std::string &Name, Tendencies *Tend,
                      AuxiliaryState *AuxState, HorzMesh *Mesh, Halo *MeshHalo);

   void doStep(OceanState *State, Real Time, Real TimeStep) const override;

};

} // namespace OMEGA
#endif
