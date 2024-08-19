#ifndef OMEGA_TSRK2_H
#define OMEGA_TSRK2_H
//===-- timeStepping/RungeKutta2Stepper.h - second-order Runge Kutta time
// stepper --------------------*- C++
//-*-===//
//
/// \file
/// \brief Contains the class for the midpoint Runge Kutta scheme
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"

namespace OMEGA {

class RungeKutta2Stepper : public TimeStepper {
 public:
   // Constructor. Construct a midpoint Runge Kutta stepper from
   // name, tendencies, auxiliary state, mesh, and halo
   RungeKutta2Stepper(const std::string &Name, Tendencies *Tend,
                      AuxiliaryState *AuxState, HorzMesh *Mesh, Halo *MeshHalo);

   // Advance the state by one step of the midpoint Runge Kutta scheme
   void doStep(OceanState *State, Real Time, Real TimeStep) const override;
};

} // namespace OMEGA
#endif
