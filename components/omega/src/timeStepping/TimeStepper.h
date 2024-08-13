#ifndef OMEGA_TIMESTEPPER_H
#define OMEGA_TIMESTEPPER_H

#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "TendencyTerms.h"

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

enum class TimeStepperType { ForwardBackward, RungeKutta4 };

class TimeStepper {
 public:
   virtual ~TimeStepper() = default;
   virtual void doStep(OceanState *State, Real Time, Real TimeStep) const = 0;

   static int init();

   static TimeStepper *create(const std::string &Name, TimeStepperType Type,
                              Tendencies *Tend, AuxiliaryState *AuxState,
                              HorzMesh *Mesh, Halo *MeshHalo);

   /// Get the default time stepper
   static TimeStepper *getDefault();

   /// Get time stepper by name
   static TimeStepper *get(const std::string &Name);

   /// Remove time stepper by name
   static void erase(const std::string &Name);

   /// Remove all time steppers
   static void clear();

   std::string getName() const;
   TimeStepperType getType() const;

   // these should be protected, they are public only because of CUDA
   // limitations
   void updateStateByTend(OceanState *State1, int TimeLevel1,
                          OceanState *State2, int TimeLevel2, Real Coeff) const;
   void updateThicknessByTend(OceanState *State1, int TimeLevel1,
                              OceanState *State2, int TimeLevel2,
                              Real Coeff) const;
   void updateVelocityByTend(OceanState *State1, int TimeLevel1,
                             OceanState *State2, int TimeLevel2,
                             Real Coeff) const;

   void updateStateByTend(OceanState *State, int TimeLevel, Real Coeff) const;
   void updateThicknessByTend(OceanState *State, int TimeLevel,
                              Real Coeff) const;
   void updateVelocityByTend(OceanState *State, int TimeLevel,
                             Real Coeff) const;

 protected:
   std::string Name;
   TimeStepperType Type;

   Tendencies *Tend;
   AuxiliaryState *AuxState;
   HorzMesh *Mesh;
   Halo *MeshHalo;

   TimeStepper(const std::string &Name, TimeStepperType Type, Tendencies *Tend,
               AuxiliaryState *AuxState, HorzMesh *Mesh, Halo *MeshHalo);

   TimeStepper(const TimeStepper &) = delete;
   TimeStepper(TimeStepper &&)      = delete;

 private:
   static TimeStepper *DefaultTimeStepper;
   static std::map<std::string, std::unique_ptr<TimeStepper>> AllTimeSteppers;
};

} // namespace OMEGA
#endif
