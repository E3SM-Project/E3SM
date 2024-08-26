#ifndef OMEGA_TIMESTEPPER_H
#define OMEGA_TIMESTEPPER_H
//===-- timeStepping/TimeStepper.h - time stepper --------------------*- C++
//-*-===//
//
/// \file
/// \brief Contains the base class for all Omega time steppers
///
/// The TimeStepper class defines the interface of a time stepper and contains
/// data and methods common to all time steppers
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

/// An enum for every time stepper type
/// needs to extended every time a new time stepper is added
enum class TimeStepperType { ForwardBackward, RungeKutta4, RungeKutta2 };

//------------------------------------------------------------------------------
/// A base class for Omega time steppers
///
/// The TimeStepper class defines the interface of a time stepper, manages
/// time stepper objects and contains common routines for state updates
class TimeStepper {
 public:
   virtual ~TimeStepper() = default;

   // The main method that every time stepper needs to define. Advances state by
   // by one time step, from Time to Time + TimeStep
   virtual void doStep(OceanState *State, TimeInstant Time) const = 0;

   /// Initialize the default time stepper
   static int init();

   // Create a time stepper from name, type, tendencies, auxiliary state, mesh
   // and halo
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

   /// Get name of time stepper from instance
   std::string getName() const;

   /// Get type (enum) of time stepper from instance
   TimeStepperType getType() const;

   /// Get time step
   TimeInterval getTimeStep() const;

   /// Set time step
   void setTimeStep(const TimeInterval &TimeStepIn);

   // these should be protected, they are public only because of CUDA
   // limitations

   // State1(TimeLevel1) = State2(TimeLevel2) + Coeff * Tend
   void updateStateByTend(OceanState *State1, int TimeLevel1,
                          OceanState *State2, int TimeLevel2,
                          TimeInterval Coeff) const;

   // LayerThickness1(TimeLevel1) = LayerThickness2(TimeLevel2) + Coeff *
   // LayerThicknessTend
   void updateThicknessByTend(OceanState *State1, int TimeLevel1,
                              OceanState *State2, int TimeLevel2,
                              TimeInterval Coeff) const;

   // NormalVelocity1(TimeLevel1) = NormalVelocity2(TimeLevel2) + Coeff *
   // NormalVelocityTend
   void updateVelocityByTend(OceanState *State1, int TimeLevel1,
                             OceanState *State2, int TimeLevel2,
                             TimeInterval Coeff) const;

 protected:
   // Name of time stepper
   std::string Name;

   // Type of time stepper
   TimeStepperType Type;

   // Time step
   TimeInterval TimeStep;

   // Pointers to objects needed by every time stepper
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
