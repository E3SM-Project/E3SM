//===--- TimeStepper.cpp - time stepper methods -------------*- C++ -*-----===//
//
// The TimeStepper defines how a solution advances forward in time.
//
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"
#include "Config.h"
#include "ForwardBackwardStepper.h"
#include "RungeKutta2Stepper.h"
#include "RungeKutta4Stepper.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// create the static class members
// Default model time stepper
TimeStepper *TimeStepper::DefaultTimeStepper = nullptr;
// All defined time steppers
std::map<std::string, std::unique_ptr<TimeStepper>>
    TimeStepper::AllTimeSteppers;

//------------------------------------------------------------------------------
// utility functions
// convert string into TimeStepperType enum
TimeStepperType getTimeStepperFromStr(const std::string &InString) {

   // Initialize TimeStepperChoice with Invalid
   TimeStepperType TimeStepperChoice = TimeStepperType::Invalid;

   if (InString == "Forward-Backward") {
      TimeStepperChoice = TimeStepperType::ForwardBackward;
   } else if (InString == "RungeKutta4") {
      TimeStepperChoice = TimeStepperType::RungeKutta4;
   } else if (InString == "RungeKutta2") {
      TimeStepperChoice = TimeStepperType::RungeKutta2;
   } else {
      LOG_CRITICAL("TimeStepper should be one of 'Forward-Backward', "
                   "'RungeKutta4' or 'RungeKutta2' but got {}:",
                   InString);
   }

   return TimeStepperChoice;
}

//------------------------------------------------------------------------------
// Constructors and creation methods.

// Constructor creates a new instance and fills in the time
// related data. attachData function is used to add the data pointers
TimeStepper::TimeStepper(
    const std::string &InName,      // [in] name of time stepper
    TimeStepperType InType,         // [in] type (time stepping method)
    I4 InNTimeLevels,               // [in] num time levels needed by method
    const TimeInstant &InStartTime, // [in] start time for time stepping
    const TimeInstant &InStopTime,  // [in] stop  time for time stepping
    const TimeInterval &InTimeStep  // [in] time step
    )
    : Name(InName), Type(InType), NTimeLevels(InNTimeLevels),
      StartTime(InStartTime), StopTime(InStopTime), TimeStep(InTimeStep) {
   // Most variables initialized via initializer list

   // Set up clock associated with this time stepper
   StepClock = std::make_unique<Clock>(Clock(InStartTime, InTimeStep));

   // Create an EndAlarm associated with the StopTime
   std::string AlarmName = "EndAlarm";
   if (InName != "Default")
      AlarmName += InName;
   EndAlarm = std::make_unique<Alarm>(Alarm(AlarmName, InStopTime));
   int Err  = StepClock->attachAlarm(EndAlarm.get());
   if (Err != 0)
      LOG_CRITICAL("Error attaching EndAlarm to TimeStep Clock");
}

//------------------------------------------------------------------------------
// Create a time stepper when all components are known
TimeStepper *TimeStepper::create(
    const std::string &InName,      ///< [in] name of time stepper
    TimeStepperType InType,         ///< [in] type (time stepping method)
    const TimeInstant &InStartTime, ///< [in] start time for time stepping
    const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
    const TimeInterval &InTimeStep, ///< [in] time step
    Tendencies *InTend,             ///< [in] ptr to tendencies
    AuxiliaryState *InAuxState,     ///< [in] ptr to aux state variables
    HorzMesh *InMesh,               ///< [in] ptr to mesh information
    Halo *InMeshHalo                ///< [in] ptr to halos
) {

   // Check for duplicates
   if (AllTimeSteppers.find(InName) != AllTimeSteppers.end()) {
      LOG_ERROR("Attempted to create a new TimeStepper with name {} but it "
                "already exists",
                InName);
      return nullptr;
   }

   TimeStepper *NewTimeStepper;

   // Call specific constructor with time info
   // These constructors mostly just call the constructor above with some
   // additional info (NTimeLevels)
   switch (InType) {
   case TimeStepperType::ForwardBackward:
      NewTimeStepper = new ForwardBackwardStepper(InName, InStartTime,
                                                  InStopTime, InTimeStep);
      break;
   case TimeStepperType::RungeKutta4:
      NewTimeStepper =
          new RungeKutta4Stepper(InName, InStartTime, InStopTime, InTimeStep);
      break;
   case TimeStepperType::RungeKutta2:
      NewTimeStepper =
          new RungeKutta2Stepper(InName, InStartTime, InStopTime, InTimeStep);
      break;
   }

   // Attach data pointers
   NewTimeStepper->attachData(InTend, InAuxState, InMesh, InMeshHalo);

   // Store instance
   AllTimeSteppers.emplace(InName, NewTimeStepper);

   return NewTimeStepper;
}

//------------------------------------------------------------------------------
// Create a time stepper when time information is needed before state
// and tendencies are defined. It creates an instance and only fills
// the time information. Data pointers are attached later.
TimeStepper *TimeStepper::create(
    const std::string &InName, // [in] name of time stepper
    TimeStepperType InType,    // [in] type (time stepping method)
    TimeInstant &InStartTime,  // [in] start time for time stepping
    TimeInstant &InStopTime,   // [in] stop  time for time stepping
    TimeInterval &InTimeStep   // [in] time step
) {

   // Check for duplicates
   if (AllTimeSteppers.find(InName) != AllTimeSteppers.end()) {
      LOG_ERROR("Attempted to create a new TimeStepper with name {} but it "
                "already exists",
                InName);
      return nullptr;
   }

   TimeStepper *NewTimeStepper;

   // Call specific constructor with time info
   switch (InType) {
   case TimeStepperType::ForwardBackward:
      NewTimeStepper = new ForwardBackwardStepper(InName, InStartTime,
                                                  InStopTime, InTimeStep);
      break;
   case TimeStepperType::RungeKutta4:
      NewTimeStepper =
          new RungeKutta4Stepper(InName, InStartTime, InStopTime, InTimeStep);
      break;
   case TimeStepperType::RungeKutta2:
      NewTimeStepper =
          new RungeKutta2Stepper(InName, InStartTime, InStopTime, InTimeStep);
      break;
   }

   // Store instance
   AllTimeSteppers.emplace(InName, NewTimeStepper);

   return NewTimeStepper;
}

//------------------------------------------------------------------------------
// For 2-step creation, this attaches all the data pointers to an instance
// once the data and tendencies have been created.
void TimeStepper::attachData(
    Tendencies *InTend,         // [in] ptr to tendencies (right hand side)
    AuxiliaryState *InAuxState, // [in] ptr to needed aux state variables
    HorzMesh *InMesh,           // [in] ptr to mesh information
    Halo *InMeshHalo            // [in] ptr to halos
) {

   if (!InTend)
      LOG_CRITICAL("Tend pointer not defined");
   if (!InAuxState)
      LOG_CRITICAL("AuxState pointer not defined");
   if (!InMesh)
      LOG_CRITICAL("HorzMesh pointer not defined");
   if (!InMeshHalo)
      LOG_CRITICAL("MeshHalo pointer not defined");

   Tend     = InTend;
   AuxState = InAuxState;
   Mesh     = InMesh;
   MeshHalo = InMeshHalo;

   // Some time steppers have additional tasks to finalize
   finalizeInit();
}

//------------------------------------------------------------------------------
// Destructors or delete functions

// Remove time stepper by name
void TimeStepper::erase(const std::string &Name) {
   AllTimeSteppers.erase(Name);
}

// Remove all time steppers
void TimeStepper::clear() { AllTimeSteppers.clear(); }

//------------------------------------------------------------------------------
// Initialize the default time stepper in two phases

// Begin initialization of the default time stepper (phase 1)
// This is primarily the time information.
int TimeStepper::init1() {
   int Err = 0;

   // Retrieve TimeStepper options from Config if available
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimeIntConfig("TimeIntegration");
   Err = OmegaConfig->get(TimeIntConfig);
   if (Err != 0) {
      LOG_CRITICAL("TimeIntegration group not found in Config");
      return Err;
   }

   // Must initialize the calendar first
   std::string CalendarStr;
   Err = TimeIntConfig.get("CalendarType", CalendarStr);
   if (Err != 0) {
      LOG_CRITICAL("CalendarType not found in TimeIntegration Config");
      return Err;
   }
   Calendar::init(CalendarStr);

   // Initialize choice of time stepper
   std::string TimeStepperStr;
   Err = TimeIntConfig.get("TimeStepper", TimeStepperStr);
   if (Err != 0) {
      LOG_CRITICAL("TimeStepper not found in TimeIntegration Config");
      return Err;
   }
   TimeStepperType TimeStepperChoice = getTimeStepperFromStr(TimeStepperStr);

   // Initialize time step
   std::string TimeStepStr;
   Err = TimeIntConfig.get("TimeStep", TimeStepStr);
   if (Err != 0) {
      LOG_CRITICAL("TimeStep not found in TimeIntegration Config");
      return Err;
   }
   TimeInterval TimeStep(TimeStepStr);

   // Initialize start time
   std::string StartTimeStr;
   Err = TimeIntConfig.get("StartTime", StartTimeStr);
   if (Err != 0) {
      LOG_CRITICAL("StartTime not found in TimeIntConfig");
      return Err;
   }
   TimeInstant StartTime(StartTimeStr);

   // Either the StopTime or RunDuration will be used to set the StopTime
   int Err1 = 0;
   std::string StopTimeStr;
   std::string DurationStr;
   Err  = TimeIntConfig.get("StopTime", StopTimeStr);
   Err1 = TimeIntConfig.get("RunDuration", DurationStr);

   // Check for empty or none strings for either choice
   if (Err == 0) {
      if (StopTimeStr == "" or StopTimeStr == " " or StopTimeStr == "none" or
          StopTimeStr == "None")
         Err = 1;
   }
   if (Err1 == 0) {
      if (DurationStr == "" or DurationStr == " " or DurationStr == "none" or
          DurationStr == "None")
         Err1 = 1;
   }
   if (Err != 0 and Err1 != 0) {
      LOG_CRITICAL("Either StopTime or RunDuration must be supplied in"
                   "TimeIntegration Config");
      return Err + Err1;
   }

   // Set stop time if a valid value is present. If both are present and
   // valid, we use the RunDuration, so compute stop time based on that first
   TimeInstant StopTime;
   if (Err1 == 0) { // valid RunDuration supplied
      TimeInterval Duration(DurationStr);
      StopTime = StartTime + Duration;
   }
   if (Err == 0 and Err1 != 0) { // only valid StopTime supplied
      TimeInstant StopTime2(StopTimeStr);
      StopTime = StopTime2;
   }
   Err = 0; // reset return code

   // Now that all the inputs are defined, create the default time stepper
   // Use the partial creation function for only the time info. Data
   // pointers will be attached in phase 2 initialization
   TimeStepper::DefaultTimeStepper =
       create("Default", TimeStepperChoice, StartTime, StopTime, TimeStep);

   return Err;
}

//------------------------------------------------------------------------------
// Finish initialization of the default time stepper (phase 2)
int TimeStepper::init2() {
   int Err = 0;

   // Get default pointers
   HorzMesh *DefMesh        = HorzMesh::getDefault();
   Halo *DefHalo            = Halo::getDefault();
   Tendencies *DefTend      = Tendencies::getDefault();
   AuxiliaryState *AuxState = AuxiliaryState::getDefault();

   // Attach data pointers
   DefaultTimeStepper->attachData(DefTend, AuxState, DefMesh, DefHalo);

   return Err;
}

//------------------------------------------------------------------------------
// Change time step
void TimeStepper::changeTimeStep(const TimeInterval &TimeStepIn) {
   TimeStep = TimeStepIn;
   int Err  = StepClock->changeTimeStep(TimeStepIn);
   if (Err != 0)
      LOG_CRITICAL("Error changing clock time step");
}

//------------------------------------------------------------------------------
// Retrieval functions

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

// Get time stepper name
std::string TimeStepper::getName() const { return Name; }

// Get time stepper type
TimeStepperType TimeStepper::getType() const { return Type; }

// Get number of time level
int TimeStepper::getNTimeLevels() const { return NTimeLevels; }

// Get time step
TimeInterval TimeStepper::getTimeStep() const { return TimeStep; }

// Get start time
TimeInstant TimeStepper::getStartTime() const { return StartTime; }

// Get stop time from instance
TimeInstant TimeStepper::getStopTime() const { return StopTime; }

// Get clock (ptr) from instance
Clock *TimeStepper::getClock() { return StepClock.get(); }

// Get end alarm (ptr) from instance
Alarm *TimeStepper::getEndAlarm() { return EndAlarm.get(); }

//------------------------------------------------------------------------------
// Update functions

//------------------------------------------------------------------------------
// Updates layer thickness using tendency terms
// LayerThickness1(TimeLevel1) = LayerThickness2(TimeLevel2) +
//                               Coeff * LayerThicknessTend
void TimeStepper::updateThicknessByTend(OceanState *State1, int TimeLevel1,
                                        OceanState *State2, int TimeLevel2,
                                        TimeInterval Coeff) const {

   Array2DReal LayerThick1;
   Array2DReal LayerThick2;
   I4 Err;
   Err = State1->getLayerThickness(LayerThick1, TimeLevel1);
   Err = State2->getLayerThickness(LayerThick2, TimeLevel2);
   const auto &LayerThickTend = Tend->LayerThicknessTend;
   const int NVertLevels      = LayerThickTend.extent_int(1);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   parallelFor(
       "updateThickByTend", {Mesh->NCellsAll, NVertLevels},
       KOKKOS_LAMBDA(int ICell, int K) {
          LayerThick1(ICell, K) =
              LayerThick2(ICell, K) + CoeffSeconds * LayerThickTend(ICell, K);
       });
}

//------------------------------------------------------------------------------
// Updates velocity using tendency terms
// NormalVelocity1(TimeLevel1) = NormalVelocity2(TimeLevel2) +
//                               Coeff * NormalVelocityTend
void TimeStepper::updateVelocityByTend(OceanState *State1, int TimeLevel1,
                                       OceanState *State2, int TimeLevel2,
                                       TimeInterval Coeff) const {

   Array2DReal NormalVel1;
   Array2DReal NormalVel2;
   I4 Err;
   Err = State1->getNormalVelocity(NormalVel1, TimeLevel1);
   Err = State2->getNormalVelocity(NormalVel2, TimeLevel2);
   const auto &NormalVelTend = Tend->NormalVelocityTend;
   const int NVertLevels     = NormalVelTend.extent_int(1);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   parallelFor(
       "updateVelByTend", {Mesh->NEdgesAll, NVertLevels},
       KOKKOS_LAMBDA(int IEdge, int K) {
          NormalVel1(IEdge, K) =
              NormalVel2(IEdge, K) + CoeffSeconds * NormalVelTend(IEdge, K);
       });
}

//------------------------------------------------------------------------------
// Updates full non-tracer state
// State1(TimeLevel1) = State2(TimeLevel2) + Coeff * Tend
void TimeStepper::updateStateByTend(OceanState *State1, int TimeLevel1,
                                    OceanState *State2, int TimeLevel2,
                                    TimeInterval Coeff) const {
   updateThicknessByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
   updateVelocityByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
}

//------------------------------------------------------------------------------
// Updates tracers
// NextTracers = (CurTracers * LayerThickness2(TimeLevel2)) +
//               Coeff * TracersTend) / LayerThickness1(TimeLevel1)
void TimeStepper::updateTracersByTend(const Array3DReal &NextTracers,
                                      const Array3DReal &CurTracers,
                                      OceanState *State1, int TimeLevel1,
                                      OceanState *State2, int TimeLevel2,
                                      TimeInterval Coeff) const {
   int Err = 0;

   const auto &LayerThick1 = State1->LayerThickness[TimeLevel1];
   const auto &LayerThick2 = State2->LayerThickness[TimeLevel2];
   const auto &TracerTend  = Tend->TracerTend;
   const int NTracers      = TracerTend.extent(0);
   const int NVertLevels   = TracerTend.extent(2);

   R8 CoeffSeconds;
   Err = Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   parallelFor(
       "updateTracersByTend", {NTracers, Mesh->NCellsAll, NVertLevels},
       KOKKOS_LAMBDA(int L, int ICell, int K) {
          NextTracers(L, ICell, K) =
              (CurTracers(L, ICell, K) * LayerThick2(ICell, K) +
               CoeffSeconds * TracerTend(L, ICell, K)) /
              LayerThick1(ICell, K);
       });
}

//------------------------------------------------------------------------------
// couple tracer array to layer thickness
void TimeStepper::weightTracers(const Array3DReal &NextTracers,
                                const Array3DReal &CurTracers,
                                OceanState *CurState, int TimeLevel1) const {

   const Array2DReal &CurThickness = CurState->LayerThickness[TimeLevel1];
   const int NTracers              = NextTracers.extent(0);
   const int NVertLevels           = NextTracers.extent(2);

   parallelFor(
       "weightTracers", {NTracers, Mesh->NCellsAll, NVertLevels},
       KOKKOS_LAMBDA(int L, int ICell, int K) {
          NextTracers(L, ICell, K) =
              CurTracers(L, ICell, K) * CurThickness(ICell, K);
       });
}

//------------------------------------------------------------------------------
// accumulate contributions to the tracer array at the next time level from
// each Runge-Kutta stage
void TimeStepper::accumulateTracersUpdate(const Array3DReal &AccumTracer,
                                          TimeInterval Coeff) const {

   const auto &TracerTend = Tend->TracerTend;
   const int NTracers     = TracerTend.extent(0);
   const int NVertLevels  = TracerTend.extent(2);

   R8 CoeffSeconds;
   int Err = Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   parallelFor(
       "accumulateTracersUpdate", {NTracers, Mesh->NCellsAll, NVertLevels},
       KOKKOS_LAMBDA(int L, int ICell, int K) {
          AccumTracer(L, ICell, K) += CoeffSeconds * TracerTend(L, ICell, K);
       });
}

//------------------------------------------------------------------------------
// normalize tracer array so final array stores concentrations
void TimeStepper::finalizeTracersUpdate(const Array3DReal &NextTracers,
                                        OceanState *State,
                                        int TimeLevel) const {

   const Array2DReal &NextThick = State->LayerThickness[TimeLevel];
   const int NTracers           = NextTracers.extent(0);
   const int NVertLevels        = NextTracers.extent(2);

   parallelFor(
       "finalizeTracersUpdate", {NTracers, Mesh->NCellsAll, NVertLevels},
       KOKKOS_LAMBDA(int L, int ICell, int K) {
          NextTracers(L, ICell, K) /= NextThick(ICell, K);
       });
}

} // namespace OMEGA
