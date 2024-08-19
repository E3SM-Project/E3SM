//===-- timeStepping/TimeStepper.cpp - time stepper methods -------------*- C++
//-*-===//

#include "TimeStepper.h"
#include "ForwardBackwardStepper.h"
#include "RungeKutta2Stepper.h"
#include "RungeKutta4Stepper.h"

namespace OMEGA {

// create the static class members
TimeStepper *TimeStepper::DefaultTimeStepper = nullptr;
std::map<std::string, std::unique_ptr<TimeStepper>>
    TimeStepper::AllTimeSteppers;

// Constructor. Construct a time stepper from name, type, tendencies, auxiliary
// state, mesh and halo
TimeStepper::TimeStepper(const std::string &Name, TimeStepperType Type,
                         Tendencies *Tend, AuxiliaryState *AuxState,
                         HorzMesh *Mesh, Halo *MeshHalo)
    : Name(Name), Type(Type), Tend(Tend), AuxState(AuxState), Mesh(Mesh),
      MeshHalo(MeshHalo) {}

// Create a time stepper from name, type, tendencies, auxiliary state, mesh and
// halo
TimeStepper *TimeStepper::create(const std::string &Name, TimeStepperType Type,
                                 Tendencies *Tend, AuxiliaryState *AuxState,
                                 HorzMesh *Mesh, Halo *MeshHalo) {

   if (AllTimeSteppers.find(Name) != AllTimeSteppers.end()) {
      LOG_ERROR("Attempted to create a new TimeStepper with name {} but it "
                "already exists",
                Name);
      return nullptr;
   }

   TimeStepper *NewTimeStepper;

   switch (Type) {
   case TimeStepperType::ForwardBackward:
      NewTimeStepper =
          new ForwardBackwardStepper(Name, Tend, AuxState, Mesh, MeshHalo);
      break;
   case TimeStepperType::RungeKutta4:
      NewTimeStepper =
          new RungeKutta4Stepper(Name, Tend, AuxState, Mesh, MeshHalo);
      break;
   case TimeStepperType::RungeKutta2:
      NewTimeStepper =
          new RungeKutta2Stepper(Name, Tend, AuxState, Mesh, MeshHalo);
      break;
   }

   AllTimeSteppers.emplace(Name, NewTimeStepper);

   return NewTimeStepper;
}

// Initialize the default time stepper
int TimeStepper::init() {
   int Err           = 0;
   auto *DefMesh     = HorzMesh::getDefault();
   auto *DefAuxState = AuxiliaryState::getDefault();
   auto *DefHalo     = Halo::getDefault();
   auto *DefTend     = Tendencies::getDefault();

   TimeStepper::DefaultTimeStepper =
       create("Default", TimeStepperType::ForwardBackward, DefTend, DefAuxState,
              DefMesh, DefHalo);
   return Err;
}

// Get the default time stepper
TimeStepper *TimeStepper::getDefault() {
   return TimeStepper::DefaultTimeStepper;
}

// Get time stepper by name
TimeStepper *TimeStepper::get(const std::string &Name) {
   // look for an instance of this name
   auto it = AllTimeSteppers.find(Name);

   // if found, return the pointer
   if (it != AllTimeSteppers.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("TimeStepper::get: Attempt to retrieve non-existent "
                "time stepper:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
}

// Remove time stepper by name
void TimeStepper::erase(const std::string &Name) {
   AllTimeSteppers.erase(Name);
}

// Remove all time steppers
void TimeStepper::clear() { AllTimeSteppers.clear(); }

// Get time stepper name
std::string TimeStepper::getName() const { return Name; }

// Get time stepper type
TimeStepperType TimeStepper::getType() const { return Type; }

// LayerThickness1(TimeLevel1) = LayerThickness(TimeLevel2) + Coeff *
// LayerThicknessTend
void TimeStepper::updateThicknessByTend(OceanState *State1, int TimeLevel1,
                                        OceanState *State2, int TimeLevel2,
                                        Real Coeff) const {

   const auto &LayerThick1    = State1->LayerThickness[TimeLevel1];
   const auto &LayerThick2    = State2->LayerThickness[TimeLevel2];
   const auto &LayerThickTend = Tend->LayerThicknessTend;
   const int NVertLevels      = LayerThickTend.extent_int(1);

   parallelFor(
       "updateThickByTend", {Mesh->NCellsAll, NVertLevels},
       KOKKOS_LAMBDA(int ICell, int K) {
          LayerThick1(ICell, K) =
              LayerThick2(ICell, K) + Coeff * LayerThickTend(ICell, K);
       });
}

// NormalVelocity(TimeLevel1) = NormalVelocity(TimeLevel2) + Coeff *
// NormalVelocityTend
void TimeStepper::updateVelocityByTend(OceanState *State1, int TimeLevel1,
                                       OceanState *State2, int TimeLevel2,
                                       Real Coeff) const {

   const auto &NormalVel1    = State1->NormalVelocity[TimeLevel1];
   const auto &NormalVel2    = State2->NormalVelocity[TimeLevel2];
   const auto &NormalVelTend = Tend->NormalVelocityTend;
   const int NVertLevels     = NormalVelTend.extent_int(1);

   parallelFor(
       "updateVelByTend", {Mesh->NEdgesAll, NVertLevels},
       KOKKOS_LAMBDA(int IEdge, int K) {
          NormalVel1(IEdge, K) =
              NormalVel2(IEdge, K) + Coeff * NormalVelTend(IEdge, K);
       });
}

// State1(TimeLevel1) = State2(TimeLevel2) + Coeff * Tend
void TimeStepper::updateStateByTend(OceanState *State1, int TimeLevel1,
                                    OceanState *State2, int TimeLevel2,
                                    Real Coeff) const {
   updateThicknessByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
   updateVelocityByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
}

// LayerThickness(TimeLevel) += Coeff * LayerThicknessTend
void TimeStepper::updateThicknessByTend(OceanState *State, int TimeLevel,
                                        Real Coeff) const {
   updateThicknessByTend(State, TimeLevel, State, TimeLevel, Coeff);
}

// NormalVelocity(TimeLevel) += Coeff * NormalVelocityTend
void TimeStepper::updateVelocityByTend(OceanState *State, int TimeLevel,
                                       Real Coeff) const {
   updateVelocityByTend(State, TimeLevel, State, TimeLevel, Coeff);
}

// State(TimeLevel) += Coeff * Tend
void TimeStepper::updateStateByTend(OceanState *State, int TimeLevel,
                                    Real Coeff) const {
   updateStateByTend(State, TimeLevel, State, TimeLevel, Coeff);
}

} // namespace OMEGA
