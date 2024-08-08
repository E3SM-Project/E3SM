#ifndef OMEGA_TSFB_H
#define OMEGA_TSFB_H

#include "TimeStepper.h"

namespace OMEGA {

class ForwardBackwardStepper : public TimeStepper {
 public:
   ForwardBackwardStepper(const std::string &Name, Tendencies *Tend,
                          AuxiliaryState *AuxState, HorzMesh *Mesh,
                          Halo *MeshHalo);

   void doStep(OceanState *State, Real Time, Real TimeStep) const override;
};

} // namespace OMEGA
#endif
