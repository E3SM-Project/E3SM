#ifndef OMEGA_TSRK4_H
#define OMEGA_TSRK4_H
//===-- timeStepping/RungeKutta4Stepper.h - fourth-order Runge Kutta time
// stepper --------------------*- C++
//-*-===//
//
/// \file
/// \brief Contains the class for the classic fourth order Runge Kutta scheme
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"

namespace OMEGA {

class RungeKutta4Stepper : public TimeStepper {
 public:
   // Constructor. Construct a fourth order Runge Kutta stepper from
   // name, tendencies, auxiliary state, mesh, and halo
   RungeKutta4Stepper(const std::string &Name, Tendencies *Tend,
                      AuxiliaryState *AuxState, HorzMesh *Mesh, Halo *MeshHalo);

   // Advance the state by one step of the fourth-order Runge Kutta scheme
   void doStep(OceanState *State, Real Time, Real TimeStep) const override;

 private:
   // Number of stages
   static constexpr int NStages = 4;

   // Provisional state
   OceanState *ProvisState;

   // Runge-Kutta coefficients
   Real RKA[NStages];
   Real RKB[NStages];
   Real RKC[NStages];
};

} // namespace OMEGA
#endif
