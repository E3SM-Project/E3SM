//===--- TimeStepper.cpp - time stepper methods -------------*- C++ -*-----===//
//
// The TimeStepper defines how a solution advances forward in time.
//
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"
#include "Config.h"
#include "Error.h"
#include "ForwardBackwardStepper.h"
#include "Logging.h"
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
PrescribeStateType TimeStepper::DefaultPrescribeThicknessMode =
    PrescribeStateType::None;
PrescribeStateType TimeStepper::DefaultPrescribeVelocityMode =
    PrescribeStateType::None;

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
      ABORT_ERROR("TimeStepper should be one of 'Forward-Backward', "
                  "'RungeKutta4' or 'RungeKutta2' but got {}:",
                  InString);
   }

   return TimeStepperChoice;
}

PrescribeStateType
getPrescribeThicknessTypeFromStr(const std::string &InString) {

   if (InString == "None") {
      return PrescribeStateType::None;
   }
   if (InString == "Init") {
      return PrescribeStateType::Init;
   }

   ABORT_ERROR(
       "PrescribeStateType should be 'None' or 'Init' for thickness but got {}",
       InString);
   return PrescribeStateType::Invalid;
}
PrescribeStateType
getPrescribeVelocityTypeFromStr(const std::string &InString) {

   if (InString == "None") {
      return PrescribeStateType::None;
   } else if (InString == "Init") {
      return PrescribeStateType::Init;
   } else if (InString == "NonDivergent") {
      return PrescribeStateType::NonDivergent;
   } else if (InString == "Divergent") {
      return PrescribeStateType::Divergent;
   }

   ABORT_ERROR("PrescribeStateType should be 'None', 'Init', 'NonDivergent' or "
               "'Divergent' for velocity but got {}",
               InString);
   return PrescribeStateType::Invalid;
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
   StepClock->attachAlarm(EndAlarm.get());
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
    VertCoord *InVCoord,            ///< [in] ptr to vertical coordinate
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
   case TimeStepperType::Invalid:
      ABORT_ERROR("Invalid time stepping method");
   default:
      ABORT_ERROR("Unknown time stepping method");
   }

   NewTimeStepper->PrescribeThicknessMode = DefaultPrescribeThicknessMode;
   NewTimeStepper->PrescribeVelocityMode  = DefaultPrescribeVelocityMode;

   // Attach data pointers
   NewTimeStepper->attachData(InTend, InAuxState, InMesh, InVCoord, InMeshHalo);

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
   case TimeStepperType::Invalid:
      ABORT_ERROR("Invalid time stepping method");
   default:
      ABORT_ERROR("Unknown time stepping method");
   }

   NewTimeStepper->PrescribeThicknessMode = DefaultPrescribeThicknessMode;
   NewTimeStepper->PrescribeVelocityMode  = DefaultPrescribeVelocityMode;

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
    VertCoord *InVCoord,        // [in] ptr to vertical coordinate
    Halo *InMeshHalo            // [in] ptr to halos
) {

   if (!InTend)
      ABORT_ERROR("Tend pointer not defined");
   if (!InAuxState)
      ABORT_ERROR("AuxState pointer not defined");
   if (!InMesh)
      ABORT_ERROR("HorzMesh pointer not defined");
   if (!InVCoord)
      ABORT_ERROR("VertCoord pointer not defined");
   if (!InMeshHalo)
      ABORT_ERROR("MeshHalo pointer not defined");

   Tend     = InTend;
   AuxState = InAuxState;
   Mesh     = InMesh;
   VCoord   = InVCoord;
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
void TimeStepper::clear() {
   AllTimeSteppers.clear();
   DefaultTimeStepper = nullptr; // prevent dangling pointer
}

//------------------------------------------------------------------------------
// Initialize the default time stepper in two phases

// Begin initialization of the default time stepper (phase 1)
// This is primarily the time information.
void TimeStepper::init1() {

   Error Err; // error code - default to success

   // Retrieve TimeStepper options from Config if available
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimeIntConfig("TimeIntegration");
   Err = OmegaConfig->get(TimeIntConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration group not found in Config");

   // Must initialize the calendar first
   std::string CalendarStr;
   Err += TimeIntConfig.get("CalendarType", CalendarStr);
   CHECK_ERROR_ABORT(Err, "CalendarType not found in TimeIntegration Config");
   Calendar::init(CalendarStr);

   // Initialize choice of time stepper
   std::string TimeStepperStr;
   Err += TimeIntConfig.get("TimeStepper", TimeStepperStr);
   CHECK_ERROR_ABORT(Err, "TimeStepper not found in TimeIntegration Config");
   TimeStepperType TimeStepperChoice = getTimeStepperFromStr(TimeStepperStr);

   // Initialize time step
   std::string TimeStepStr;
   Err += TimeIntConfig.get("TimeStep", TimeStepStr);
   CHECK_ERROR_ABORT(Err, "TimeStep not found in TimeIntegration Config");
   TimeInterval TimeStep(TimeStepStr);

   // Initialize start time
   std::string StartTimeStr;
   Err = TimeIntConfig.get("StartTime", StartTimeStr);
   CHECK_ERROR_ABORT(Err, "StartTime not found in TimeIntConfig");
   TimeInstant StartTime(StartTimeStr);

   // Either the StopTime or RunDuration will be used to set the StopTime
   std::string StopTimeStr;
   std::string DurationStr;
   Err += TimeIntConfig.get("StopTime", StopTimeStr);
   Error Err1 = TimeIntConfig.get("RunDuration", DurationStr);

   // Check for empty or none strings for either choice
   bool ValidStopTime = false;
   bool ValidDuration = false;
   if (Err.isSuccess()) { // stop time was read from config
      ValidStopTime = true;
      if (StopTimeStr == "" or StopTimeStr == " " or StopTimeStr == "none" or
          StopTimeStr == "None")
         ValidStopTime = false;
   }
   if (Err1.isSuccess()) { // duration was read from config
      ValidDuration = true;
      if (DurationStr == "" or DurationStr == " " or DurationStr == "none" or
          DurationStr == "None")
         ValidDuration = false;
   }
   if (!ValidStopTime and !ValidDuration) {
      ABORT_ERROR("Either StopTime or RunDuration must be supplied in"
                  "TimeIntegration Config");
   }

   // Set stop time if a valid value is present. If both are present and
   // valid, we use the RunDuration, so compute stop time based on that first
   TimeInstant StopTime;
   if (ValidDuration) { // valid RunDuration supplied
      TimeInterval Duration(DurationStr);
      StopTime = StartTime + Duration;
   } else { // only valid StopTime supplied
      TimeInstant StopTime2(StopTimeStr);
      StopTime = StopTime2;
   }

   Config StateConfig("State");
   Error StateErr = OmegaConfig->get(StateConfig);
   if (StateErr.isSuccess()) {
      std::string ThicknessMode;
      if (StateConfig.get("PrescribeThicknessType", ThicknessMode)
              .isSuccess()) {
         TimeStepper::DefaultPrescribeThicknessMode =
             getPrescribeThicknessTypeFromStr(ThicknessMode);
      }

      std::string VelocityMode;
      if (StateConfig.get("PrescribeVelocityType", VelocityMode).isSuccess()) {
         TimeStepper::DefaultPrescribeVelocityMode =
             getPrescribeVelocityTypeFromStr(VelocityMode);
      }
   }

   // Now that all the inputs are defined, create the default time stepper
   // Use the partial creation function for only the time info. Data
   // pointers will be attached in phase 2 initialization
   TimeStepper::DefaultTimeStepper =
       create("Default", TimeStepperChoice, StartTime, StopTime, TimeStep);
}

//------------------------------------------------------------------------------
// Finish initialization of the default time stepper (phase 2)
void TimeStepper::init2() {

   // Get default pointers
   HorzMesh *DefMesh        = HorzMesh::getDefault();
   VertCoord *DefVCoord     = VertCoord::getDefault();
   Halo *DefHalo            = Halo::getDefault();
   Tendencies *DefTend      = Tendencies::getDefault();
   AuxiliaryState *AuxState = AuxiliaryState::getDefault();

   // Attach data pointers
   DefaultTimeStepper->attachData(DefTend, AuxState, DefMesh, DefVCoord,
                                  DefHalo);
}

//------------------------------------------------------------------------------
// Change time step
void TimeStepper::changeTimeStep(const TimeInterval &TimeStepIn) {
   TimeStep = TimeStepIn;
   StepClock->changeTimeStep(TimeStepIn);
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

   Array2DReal LayerThick1 = State1->getLayerThickness(TimeLevel1);
   Array2DReal LayerThick2 = State2->getLayerThickness(TimeLevel2);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   OMEGA_SCOPE(LayerThickTend, Tend->LayerThicknessTend);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   parallelForOuter(
       "updateThickByTend", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 LayerThick1(ICell, K) =
                     LayerThick2(ICell, K) +
                     CoeffSeconds * LayerThickTend(ICell, K);
              });
       });
}

//------------------------------------------------------------------------------
// Updates velocity using tendency terms
// NormalVelocity1(TimeLevel1) = NormalVelocity2(TimeLevel2) +
//                               Coeff * NormalVelocityTend
void TimeStepper::updateVelocityByTend(OceanState *State1, int TimeLevel1,
                                       OceanState *State2, int TimeLevel2,
                                       TimeInterval Coeff) const {

   Array2DReal NormalVel1 = State1->getNormalVelocity(TimeLevel1);
   Array2DReal NormalVel2 = State2->getNormalVelocity(TimeLevel2);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   OMEGA_SCOPE(NormalVelTend, Tend->NormalVelocityTend);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelForOuter(
       "updateVelByTend", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin = MinLayerEdgeBot(IEdge);
          const int KMax = MaxLayerEdgeTop(IEdge);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 NormalVel1(IEdge, K) = NormalVel2(IEdge, K) +
                                        CoeffSeconds * NormalVelTend(IEdge, K);
              });
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
// Reset state variables to their initial values
void TimeStepper::prescribeThickness(OceanState *State1, int TimeLevel1,
                                     OceanState *State2, int TimeLevel2) const {

   if (PrescribeThicknessMode == PrescribeStateType::None) {
      return;
   }

   if (PrescribeThicknessMode == PrescribeStateType::Init) {
      Array2DReal LayerThick1 = State1->getLayerThickness(TimeLevel1);
      Array2DReal LayerThick2 = State2->getLayerThickness(TimeLevel2);

      OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
      OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

      parallelForOuter(
          "prescribeThickness", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRange(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const int K           = KMin + KChunk;
                    LayerThick1(ICell, K) = LayerThick2(ICell, K);
                 });
          });
      return;
   }
}

//------------------------------------------------------------------------------
void TimeStepper::prescribeVelocity(OceanState *State1, int TimeLevel1,
                                    OceanState *State2, int TimeLevel2,
                                    const TimeInstant &SimTime) const {

   if (PrescribeVelocityMode == PrescribeStateType::None) {
      return;
   }

   if (PrescribeVelocityMode == PrescribeStateType::Init) {
      Array2DReal NormalVel1 = State1->getNormalVelocity(TimeLevel1);
      Array2DReal NormalVel2 = State2->getNormalVelocity(TimeLevel2);

      OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
      OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

      parallelForOuter(
          "prescribeVelocity", {Mesh->NEdgesAll},
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRange(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const int K          = KMin + KChunk;
                    NormalVel2(IEdge, K) = NormalVel1(IEdge, K);
                 });
          });
      return;
   } else if (PrescribeVelocityMode == PrescribeStateType::NonDivergent) {
      Array2DReal NormalVel2 = State2->getNormalVelocity(TimeLevel2);

      OMEGA_SCOPE(LatEdge, Mesh->LatEdgeH);
      OMEGA_SCOPE(LonEdge, Mesh->LonEdgeH);
      OMEGA_SCOPE(AngleEdge, Mesh->AngleEdgeH);
      OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBotH);
      OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTopH);

      const Clock *ModelClock = StepClock.get();
      R8 ElapsedTimeSec;
      TimeInterval ElapsedTimeInterval = SimTime - ModelClock->getStartTime();
      ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

      const R8 Tau  = 12. * Day2Sec; // 12 days in seconds
      const R8 TSim = ElapsedTimeSec;

      parallelForOuter(
          "prescribeVelocityNonDivergent", {Mesh->NEdgesAll},
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRange(KMin, KMax);

             const R8 lon_p = LonEdge(IEdge) - 2.0 * Pi * TSim / Tau;
             const R8 u     = (1 / Tau) * (10.0 * Kokkos::pow(sin(lon_p), 2) *
                                           sin(2.0 * LatEdge(IEdge)) *
                                           cos(Pi * TSim / Tau) +
                                       2.0 * Pi * cos(LatEdge(IEdge)));
             const R8 v     = (10.0 / Tau) * sin(2.0 * lon_p) *
                          cos(LatEdge(IEdge)) * cos(Pi * TSim / Tau);
             const R8 normalVel = REarth * (u * cos(AngleEdge(IEdge)) +
                                            v * sin(AngleEdge(IEdge)));

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const int K          = KMin + KChunk;
                    NormalVel2(IEdge, K) = normalVel;
                 });
          });
      return;
   } else if (PrescribeVelocityMode == PrescribeStateType::Divergent) {
      Array2DReal NormalVel2 = State2->getNormalVelocity(TimeLevel2);

      OMEGA_SCOPE(LatEdge, Mesh->LatEdgeH);
      OMEGA_SCOPE(LonEdge, Mesh->LonEdgeH);
      OMEGA_SCOPE(AngleEdge, Mesh->AngleEdgeH);
      OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBotH);
      OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTopH);

      const Clock *ModelClock = StepClock.get();
      R8 ElapsedTimeSec;
      TimeInterval ElapsedTimeInterval = SimTime - ModelClock->getStartTime();
      ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

      const R8 Tau  = 12. * Day2Sec; // 14 days in seconds
      const R8 TSim = ElapsedTimeSec;

      parallelForOuter(
          "prescribeVelocityDivergent", {Mesh->NEdgesAll},
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRange(KMin, KMax);

             const R8 lon_p = LonEdge(IEdge) - 2.0 * Pi * TSim / Tau;
             const R8 u =
                 (1.0 / Tau) * (-5.0 * Kokkos::pow(sin(lon_p / 2), 2) *
                                    sin(2.0 * LatEdge(IEdge)) *
                                    Kokkos::pow(cos(LatEdge(IEdge)), 2) *
                                    cos(Pi * TSim / Tau) +
                                2.0 * Pi * cos(LatEdge(IEdge)));
             const R8 v =
                 ((2.5 / Tau) * sin(lon_p) *
                  Kokkos::pow(cos(LatEdge(IEdge)), 3) * cos(Pi * TSim / Tau));
             const R8 normalVel = REarth * (u * cos(AngleEdge(IEdge)) +
                                            v * sin(AngleEdge(IEdge)));

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const int K          = KMin + KChunk;
                    NormalVel2(IEdge, K) = normalVel;
                 });
          });
      return;
   }
}

//------------------------------------------------------------------------------
void TimeStepper::prescribeState(OceanState *State1, int TimeLevel1,
                                 OceanState *State2, int TimeLevel2,
                                 const TimeInstant &SimTime) const {
   prescribeThickness(State1, TimeLevel1, State2, TimeLevel2);
   prescribeVelocity(State1, TimeLevel1, State2, TimeLevel2, SimTime);
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

   Array2DReal LayerThick1 = State1->getLayerThickness(TimeLevel1);
   Array2DReal LayerThick2 = State2->getLayerThickness(TimeLevel2);

   OMEGA_SCOPE(TracerTend, Tend->TracerTend);
   const int NTracers = TracerTend.extent(0);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   parallelForOuter(
       "updateTracersByTend", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 NextTracers(L, ICell, K) =
                     (CurTracers(L, ICell, K) * LayerThick2(ICell, K) +
                      CoeffSeconds * TracerTend(L, ICell, K)) /
                     LayerThick1(ICell, K);
              });
       });
}

//------------------------------------------------------------------------------
// couple tracer array to layer thickness
void TimeStepper::weightTracers(const Array3DReal &NextTracers,
                                const Array3DReal &CurTracers,
                                OceanState *CurState, int TimeLevel1) const {

   Array2DReal CurThickness = CurState->getLayerThickness(TimeLevel1);
   const int NTracers       = NextTracers.extent(0);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   parallelForOuter(
       "weightTracers", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 NextTracers(L, ICell, K) =
                     CurTracers(L, ICell, K) * CurThickness(ICell, K);
              });
       });
}

//------------------------------------------------------------------------------
// accumulate contributions to the tracer array at the next time level from
// each Runge-Kutta stage
void TimeStepper::accumulateTracersUpdate(const Array3DReal &AccumTracer,
                                          TimeInterval Coeff) const {

   const auto &TracerTend = Tend->TracerTend;
   const int NTracers     = TracerTend.extent(0);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   parallelForOuter(
       "accumulateTracersUpdate", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 AccumTracer(L, ICell, K) +=
                     CoeffSeconds * TracerTend(L, ICell, K);
              });
       });
}

//------------------------------------------------------------------------------
// normalize tracer array so final array stores concentrations
void TimeStepper::finalizeTracersUpdate(const Array3DReal &NextTracers,
                                        OceanState *State,
                                        int TimeLevel) const {

   Array2DReal NextThick = State->getLayerThickness(TimeLevel);
   const int NTracers    = NextTracers.extent(0);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   parallelForOuter(
       "finalizeTracersUpdate", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 NextTracers(L, ICell, K) /= NextThick(ICell, K);
              });
       });
}

} // namespace OMEGA
