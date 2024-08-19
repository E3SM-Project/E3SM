#ifndef OMEGA_TSFB_H
#define OMEGA_TSFB_H
//===-- timeStepping/ForwardBackwardStepper.h - forward-backward time stepper
//--------------------*- C++
//-*-===//
//
/// \file
/// \brief Contains the class for forward-backward time stepping scheme
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"

namespace OMEGA {

class ForwardBackwardStepper : public TimeStepper {
 public:
   // Constructor. Construct a forward-backward stepper from
   // name, tendencies, auxiliary state, mesh, and halo
   ForwardBackwardStepper(const std::string &Name, Tendencies *Tend,
                          AuxiliaryState *AuxState, HorzMesh *Mesh,
                          Halo *MeshHalo);

   // Advance the state by one step of the forward-backward scheme
   void doStep(OceanState *State, Real Time, Real TimeStep) const override;
};

} // namespace OMEGA
#endif
