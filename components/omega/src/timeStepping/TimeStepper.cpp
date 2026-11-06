//===--- TimeStepper.cpp - time stepper methods -------------*- C++ -*-----===//
//
// The TimeStepper defines how a solution advances forward in time.
//
//===----------------------------------------------------------------------===//

#include "TimeStepper.h"
#include "Config.h"
#include "Error.h"
#include "ForwardBackwardStepper.h"
#include "Frazil.h"
#include "KPPMix.h"
#include "Logging.h"
#include "Pacer.h"
#include "RungeKutta2Stepper.h"
#include "RungeKutta4Stepper.h"
#include "SplitExplicitRK2Stepper.h"
#include "VertMix.h"

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
//------------------------------------------------------------------------------
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
   } else if (InString == "SplitExplicitRK2") {
      TimeStepperChoice = TimeStepperType::SplitExplicitRK2;
   } else if (InString == "UnsplitRK2") {
      TimeStepperChoice = TimeStepperType::UnsplitRK2;
   } else {
      ABORT_ERROR("TimeStepper should be one of 'Forward-Backward', "
                  "'RungeKutta4', 'RungeKutta2', 'SplitExplicitRK2' or "
                  "'UnsplitRK2' but got {}:",
                  InString);
   }

   return TimeStepperChoice;
}

//------------------------------------------------------------------------------
// convert string into StartType enum
TimeStepperStartType
getTimeStepperStartTypeFromStr(const std::string &InString) {

   // Convert input string to lowercase for easier and more robust comparison
   std::string StartStr = InString;
   std::transform(StartStr.begin(), StartStr.end(), StartStr.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   TimeStepperStartType StartChoice;
   if (StartStr == "startup" or StartStr == "init" or StartStr == "initial") {
      StartChoice = TimeStepperStartType::StartUp;
   } else if (StartStr == "continue" or StartStr == "restart") {
      StartChoice = TimeStepperStartType::Continue;
   } else if (StartStr == "branch") {
      StartChoice = TimeStepperStartType::Branch;
   } else if (StartStr == "hybrid") {
      StartChoice = TimeStepperStartType::Hybrid;
   } else if (StartStr == "coupled") {
      StartChoice = TimeStepperStartType::Coupled;
   } else {
      StartChoice = TimeStepperStartType::Invalid;
      ABORT_ERROR("Invalid StartType {}", InString);
   }

   return StartChoice;
}

//------------------------------------------------------------------------------
// convert E3SM start type into StartType enum
TimeStepperStartType getTimeStepperStartTypeFromE3SM(const int E3SMOption) {

   // Translate the integer start option from the E3SM coupler to the
   // internal enum. For coupled simulations, the run type and continue
   // options have been translated by CIME into these three startup choices.

   TimeStepperStartType StartChoice;
   switch (E3SMOption) {
   case 0:
      StartChoice = TimeStepperStartType::StartUp;
      break;
   case 1:
      StartChoice = TimeStepperStartType::Continue;
      break;
   case 2:
      StartChoice = TimeStepperStartType::Branch;
      break;
   default:
      StartChoice = TimeStepperStartType::Invalid;
      ABORT_ERROR("Invalid E3SM start type value: {}", E3SMOption);
   }
   return StartChoice;
}

//------------------------------------------------------------------------------
// convert string into StopType enum
TimeStepperStopType getTimeStepperStopTypeFromStr(const std::string &InString) {

   // Convert input string to lowercase for easier and more robust comparison
   std::string StopStr = InString;
   std::transform(StopStr.begin(), StopStr.end(), StopStr.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   TimeStepperStopType StopChoice;
   if (StopStr == "attime") {
      StopChoice = TimeStepperStopType::AtTime;
   } else if (StopStr == "afterduration") {
      StopChoice = TimeStepperStopType::AfterDuration;
   } else if (StopStr == "coupled") {
      StopChoice = TimeStepperStopType::Coupled;
   } else {
      StopChoice = TimeStepperStopType::Invalid;
      ABORT_ERROR("Invalid StopType {}", InString);
   }

   return StopChoice;
}

//------------------------------------------------------------------------------
// Convert string into PrescribeStateType enum
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

// Constructor creates a new instance and fills in most of the time related
// data. The attachData function is used to add the data pointers once they
// are known. While an initial StopTime and EndAlarm are created here, they
// must be reset later if the initial time for the segment is updated (eg during
// the restart read).
TimeStepper::TimeStepper(
    const std::string &InName,              // [in] name of time stepper
    TimeStepperType InType,                 // [in] type (time stepping method)
    I4 InNTimeLevels,                       // [in] num time levels for method
    const TimeInterval &InTimeStep,         // [in] time step
    const TimeStepperStartType InStartType, // [in] option for starting
    const TimeInstant &InStartTime,         // [in] start time for full sim
    const TimeStepperStopType InStopType,   // [in] option for stopping
    std::optional<TimeInstant> InStopTime,  // [in] stop time if opt AtTime
    std::optional<TimeInterval> InDuration
    // [in] duration of simulation segment if StopType is AfterDuration
    )
    : Name(InName), Type(InType), NTimeLevels(InNTimeLevels),
      TimeStep(InTimeStep), StartType(InStartType), StartTime(InStartTime),
      StopType(InStopType) {
   // Many variables initialized via initializer list

   // Set up clock associated with this time stepper
   StepClock = std::make_unique<Clock>(Clock(InStartTime, InTimeStep));

   // Create an initial stop time and end alarm and attach the alarm to the
   // StepClock. For the duration case, these are only valid for the first run
   // segment and must be reset later after the current time is updated on
   // restart.

   // Create alarm name based on TimeStepper instance name
   std::string AlarmName = "EndAlarm";
   if (InName != "Default")
      AlarmName += InName;

   switch (StopType) {
   case TimeStepperStopType::AtTime:
      if (InStopTime.has_value()) { // StopTime must be provided
         StopTime = InStopTime.value();
         EndAlarm = std::make_unique<Alarm>(Alarm(AlarmName, StopTime));
         StepClock->attachAlarm(EndAlarm.get());
         Duration = StopTime - StartTime;
      } else {
         ABORT_ERROR("While creating TimeStepper {}, the StopType AtTime"
                     " was requested but a StopTime was not provided",
                     Name);
      }
      break;
   case TimeStepperStopType::AfterDuration:
      // This initialization of the end alarm is only valid for a simulation
      // starting from scratch or branching with a time reset.
      // If the simulation is being continued from a restart file, then
      // these must be reset later using the resetEndAlarm function.
      if (InDuration.has_value()) { // Duration must be provided
         TimeInstant CurrentTime = StepClock->getCurrentTime();
         Duration                = InDuration.value();
         StopTime                = CurrentTime + Duration;
         EndAlarm = std::make_unique<Alarm>(Alarm(AlarmName, StopTime));
         StepClock->attachAlarm(EndAlarm.get());
      } else {
         ABORT_ERROR("While creating TimeStepper {}, the StopType"
                     " AfterDuration was requested but a Duration was not"
                     " provided",
                     Name);
      }
      break;
   case TimeStepperStopType::Coupled: {
      // Simulation will stop on an external signal so no StopTime
      // or Duration are needed.
      TimeInstant CurrentTime = StepClock->getCurrentTime();
      // Set Duration and StopTime with long future
      // values for a dummy EndAlarm
      Duration = TimeInterval(1.e16, TimeUnits::Seconds);
      StopTime = CurrentTime + Duration;
      EndAlarm = std::make_unique<Alarm>(Alarm(AlarmName, StopTime));
      StepClock->attachAlarm(EndAlarm.get());
   } break;
   default:
      ABORT_ERROR("Invalid StopType encountered creating TimeStepper {}", Name);
   }
}

//------------------------------------------------------------------------------
// Create a time stepper when all components are known
// Note that if StopType is AfterDuration, the StopTime and EndAlarm
// are computed based on the clock's current time and must be reset if
// the current time is reset (eg by reading a restart) using resetEndAlarm
TimeStepper *TimeStepper::create(
    const std::string &InName,              // [in] name of time stepper
    TimeStepperType InType,                 // [in] type (time stepping method)
    Tendencies *InTend,                     // [in] ptr to tendencies
    AuxiliaryState *InAuxState,             // [in] ptr to aux state variables
    HorzMesh *InMesh,                       // [in] ptr to mesh information
    VertCoord *InVCoord,                    // [in] ptr to vertical coordinate
    Halo *InMeshHalo,                       // [in] ptr to halos
    const TimeInterval &InTimeStep,         // [in] time step
    const TimeStepperStartType InStartType, // [in] option to start sim
    const TimeInstant &InStartTime,         // [in] full simulation start time
    const TimeStepperStopType InStopType,   // [in] option to stop
    std::optional<TimeInstant> InStopTime,  // [in] stop time if option AtTime
    std::optional<TimeInterval> InDuration  // [in] duration if opt AfterDur
) {

   OMEGA_REQUIRE(
       InTend, "Null Tendencies pointer in TimeStepper::create with Name = {}",
       InName);
   OMEGA_REQUIRE(
       InAuxState,
       "Null AuxiliaryState pointer in TimeStepper::create with Name = {}",
       InName);
   OMEGA_REQUIRE(InMesh,
                 "Null HorzMesh pointer in TimeStepper::create with Name = {}",
                 InName);
   OMEGA_REQUIRE(InVCoord,
                 "Null VertCoord pointer in TimeStepper::create with Name = {}",
                 InName);
   OMEGA_REQUIRE(InMeshHalo,
                 "Null Halo pointer in TimeStepper::create with Name = {}",
                 InName);

   // Start by calling the two-phase create function
   TimeStepper *NewTimeStepper =
       create(InName, InType, InTimeStep, InStartType, InStartTime, InStopType,
              InStopTime, InDuration);

   NewTimeStepper->PrescribeThicknessMode = DefaultPrescribeThicknessMode;
   NewTimeStepper->PrescribeVelocityMode  = DefaultPrescribeVelocityMode;

   // Attach data pointers
   NewTimeStepper->attachData(InTend, InAuxState, InMesh, InVCoord, InMeshHalo);

   return NewTimeStepper;
}

//------------------------------------------------------------------------------
// Create a time stepper when time information is needed before state
// and tendencies are defined. It creates an instance and only fills
// the time information. Data pointers are attached later.
// Note that if StopType is AfterDuration, the StopTime and EndAlarm
// are computed based on the clock's current time and must be reset if
// the current time is reset (eg by reading a restart) using resetEndAlarm
TimeStepper *TimeStepper::create(
    const std::string &InName,              // [in] name of time stepper
    TimeStepperType InType,                 // [in] type (time stepping method)
    const TimeInterval &InTimeStep,         // [in] time step
    const TimeStepperStartType InStartType, // [in] option for starting
    const TimeInstant &InStartTime,         // [in] start time for full sim
    const TimeStepperStopType InStopType,   // [in] option for stopping
    std::optional<TimeInstant> InStopTime,  // [in] stop time option is AtTime
    std::optional<TimeInterval> InDuration  // [in] duration opt AfterDuration
) {

   // Check for duplicates
   if (AllTimeSteppers.find(InName) != AllTimeSteppers.end()) {
      ABORT_ERROR("Attempted to create a new TimeStepper with name {} but it "
                  "already exists",
                  InName);
   }

   TimeStepper *NewTimeStepper;

   // Call specific constructor with time info
   switch (InType) {
   case TimeStepperType::ForwardBackward:
      NewTimeStepper = new ForwardBackwardStepper(
          InName, InTimeStep, InStartType, InStartTime, InStopType, InStopTime,
          InDuration);
      break;
   case TimeStepperType::RungeKutta4:
      NewTimeStepper =
          new RungeKutta4Stepper(InName, InTimeStep, InStartType, InStartTime,
                                 InStopType, InStopTime, InDuration);
      break;
   case TimeStepperType::RungeKutta2:
      NewTimeStepper =
          new RungeKutta2Stepper(InName, InTimeStep, InStartType, InStartTime,
                                 InStopType, InStopTime, InDuration);
      break;
   // Both share an implementation; the type selects whether the barotropic
   // mode is split off (SplitFactor 1) or carried in the baroclinic solve
   // (SplitFactor 0).
   case TimeStepperType::SplitExplicitRK2:
   case TimeStepperType::UnsplitRK2:
      NewTimeStepper = new SplitExplicitRK2Stepper(
          InName, InType, InTimeStep, InStartType, InStartTime, InStopType,
          InStopTime, InDuration);
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
// 1st stage of two-stage initialization for the default time stepper
// This version with no input parameters is for standalone mode
// with parameters from configuration input
void TimeStepper::init1() {

   Error Err; // error code - default to success

   // Retrieve TimeStepper options from Config if available
   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in TimeStepper::init1");
   Config TimeIntConfig("TimeIntegration");
   Err = OmegaConfig->get(TimeIntConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration group not found in Config");

   // Must initialize the calendar first, before any TimeInstant objects
   std::string CalendarStr;
   Err += TimeIntConfig.get("CalendarType", CalendarStr);
   CHECK_ERROR_ABORT(Err, "CalendarType not found in TimeIntegration Config");
   Calendar::init(CalendarStr);

   // Initialize start option
   std::string StartTypeStr;
   Err = TimeIntConfig.get("StartType", StartTypeStr);
   CHECK_ERROR_ABORT(Err, "StartType not found in TimeIntegration Config");
   TimeStepperStartType InStartType =
       getTimeStepperStartTypeFromStr(StartTypeStr);

   // Initialize start time
   std::string StartTimeStr;
   Err = TimeIntConfig.get("StartTime", StartTimeStr);
   CHECK_ERROR_ABORT(Err, "StartTime not found in TimeIntConfig");
   TimeInstant StartTime(StartTimeStr);

   // Extract Prescribe options from config
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

   init1(InStartType, StartTime);
}

//------------------------------------------------------------------------------
// 1st stage of two-stage initialization for the default time stepper
// This version with start parameters is for coupled mode but is also
// called by the init1 routine above to complete standalone init
void TimeStepper::init1(const TimeStepperStartType InStartType,
                        const TimeInstant &StartTime) {

   // Calendar must be initialized before this is called — the no-arg init1()
   // does this internally. In coupled mode, the caller (ocnInit) is responsible
   // for calling Calendar::init() before constructing TimeInitParams.
   OMEGA_REQUIRE(Calendar::isDefined(),
                 "Calendar must be initialized before TimeStepper::init1");

   Error Err; // error code - default to success

   // TimeStepper and TimeStep are always read from the Config
   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in TimeStepper::init1");
   Config TimeIntConfig("TimeIntegration");
   Err = OmegaConfig->get(TimeIntConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration group not found in Config");

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

   // Initialize the option for stopping the simulation
   std::string StopTypeStr;
   Err += TimeIntConfig.get("StopType", StopTypeStr);
   CHECK_ERROR_ABORT(Err, "StopType not found in TimeIntegrationConfig");
   TimeStepperStopType InStopType = getTimeStepperStopTypeFromStr(StopTypeStr);

   // Depending on the stop option, extract the StopTime, Duration variables
   // from the StopCriterion config string
   std::string StopCriterionStr;
   std::optional<TimeInstant> InStopTime  = std::nullopt;
   std::optional<TimeInterval> InDuration = std::nullopt;
   switch (InStopType) {
   case TimeStepperStopType::AtTime: {
      // The StopCriterion string contains the StopTime string
      Err += TimeIntConfig.get("StopCriterion", StopCriterionStr);
      CHECK_ERROR_ABORT(Err, "StopCriterion not found in TimeIntConfig");
      TimeInstant StopTimeTmp(StopCriterionStr); // extract StopTime
      InStopTime = StopTimeTmp;
      break;
   }
   case TimeStepperStopType::AfterDuration: {
      // The StopCriterion string contains the Duration string
      Err += TimeIntConfig.get("StopCriterion", StopCriterionStr);
      CHECK_ERROR_ABORT(Err, "StopCriterion not found in TimeIntConfig");
      TimeInterval DurationTmp(StopCriterionStr); // extract Duration
      InDuration = DurationTmp;
      break;
   }
   case TimeStepperStopType::Coupled:
      // Neither StopTime or Duration needed
      break;
   default:
      ABORT_ERROR("Unknown StopType {} while initializing TimeStepper",
                  StopTypeStr);
   }

   // Now that all the inputs are defined, create the default time stepper
   // Use the partial creation function for only the time info. Data
   // pointers will be attached in phase 2 initialization
   TimeStepper::DefaultTimeStepper =
       create("Default", TimeStepperChoice, TimeStep, InStartType, StartTime,
              InStopType, InStopTime, InDuration);
}

//------------------------------------------------------------------------------
// Finish initialization of the default time stepper (phase 2)
void TimeStepper::init2() {

   OMEGA_REQUIRE(DefaultTimeStepper,
                 "Null default TimeStepper pointer in TimeStepper::init2");

   // Get default pointers
   HorzMesh *DefMesh = HorzMesh::getDefault();
   OMEGA_REQUIRE(DefMesh,
                 "Null default HorzMesh pointer in TimeStepper::init2");
   VertCoord *DefVCoord = VertCoord::getDefault();
   OMEGA_REQUIRE(DefVCoord,
                 "Null default VertCoord pointer in TimeStepper::init2");
   Halo *DefHalo = Halo::getDefault();
   OMEGA_REQUIRE(DefHalo, "Null default Halo pointer in TimeStepper::init2");
   Tendencies *DefTend = Tendencies::getDefault();
   OMEGA_REQUIRE(DefTend,
                 "Null default Tendencies pointer in TimeStepper::init2");
   AuxiliaryState *AuxState = AuxiliaryState::getDefault();
   OMEGA_REQUIRE(AuxState,
                 "Null default AuxiliaryState pointer in TimeStepper::init2");

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
// Get number of doStep calls made on this instance
I8 TimeStepper::getStepCount() const { return StepCount; }

void TimeStepper::resetFrazilFluxes() const {
   auto *DefaultFrazil = Frazil::getDefault();
   if (DefaultFrazil) {
      DefaultFrazil->resetOcnStepFluxes();
   }
}

//------------------------------------------------------------------------------
// Is this a split time stepper. False by default
bool TimeStepper::isSplit() const { return false; }

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

// Get start option from instance
TimeStepperStartType TimeStepper::getStartType() const { return StartType; }

// Get start time
TimeInstant TimeStepper::getStartTime() const { return StartTime; }

// Get stop option from instance
TimeStepperStopType TimeStepper::getStopType() const { return StopType; }

// Get stop time from instance
TimeInstant TimeStepper::getStopTime() const { return StopTime; }

// Get duration of run segment
TimeInterval TimeStepper::getDuration() const { return Duration; }

// Get clock (ptr) from instance
Clock *TimeStepper::getClock() { return StepClock.get(); }

// Get end alarm (ptr) from instance
Alarm *TimeStepper::getEndAlarm() { return EndAlarm.get(); }

//------------------------------------------------------------------------------
// Update functions
//------------------------------------------------------------------------------
// If the current segment start time is reset (eg by reading a restart), the
// relevant stop times and alarms must be reset.  It is assumed that the
// ModelClock CurrentTime has been reset already and uses this updated time
// as the initial time for the duration. This only modifies StopTime and
// EndAlarm if StopType is AfterDuration.
void TimeStepper::resetEndAlarm() {

   // The end alarm is only reset from the StopType AfterDuration based on
   // the new current time in the step clock
   switch (StopType) {
   case TimeStepperStopType::AtTime:
      // Stop time is fixed, so no reset
      break;
   case TimeStepperStopType::AfterDuration: {
      // We assume the CurrentTime has been updated and compute a new
      // stop time and alarm based on that CurrentTime and Duration
      TimeInstant CurrentTime = StepClock->getCurrentTime();
      StopTime                = CurrentTime + Duration;
      EndAlarm->reset(StopTime);
      break;
   }
   case TimeStepperStopType::Coupled:
      // Simulation will stop on an external signal and EndTime has been
      // set far in the future. No reset is needed.
      break;
   default:
      ABORT_ERROR("Invalid StopType encountered creating TimeStepper {}", Name);
   }
}

//------------------------------------------------------------------------------
// Updates pseudo-thickness using tendency terms
// PseudoThickness1(TimeLevel1) = PseudoThickness2(TimeLevel2) +
//                               Coeff * PseudoThicknessTend
void TimeStepper::updateThicknessByTend(OceanState *State1, int TimeLevel1,
                                        OceanState *State2, int TimeLevel2,
                                        TimeInterval Coeff) const {

   Array2DReal PseudoThick1 = State1->getPseudoThickness(TimeLevel1);
   Array2DReal PseudoThick2 = State2->getPseudoThickness(TimeLevel2);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   OMEGA_SCOPE(PseudoThickTend, Tend->PseudoThicknessTend);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   parallelForOuter(
       "updateThickByTend", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 PseudoThick1(ICell, K) =
                     PseudoThick2(ICell, K) +
                     CoeffSeconds * PseudoThickTend(ICell, K);
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
   // Zero boundary layers; updateVelocityByTend only writes [KMin,KMax], so
   // layers outside that range retain fill values that would corrupt
   // tendencies.
   VCoord->applyEdgeLayerMask(NormalVel1, Mesh->NEdgesAll);
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
      Array2DReal PseudoThick1 = State1->getPseudoThickness(TimeLevel1);
      Array2DReal PseudoThick2 = State2->getPseudoThickness(TimeLevel2);

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
                    const int K            = KMin + KChunk;
                    PseudoThick1(ICell, K) = PseudoThick2(ICell, K);
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

      OMEGA_SCOPE(LatEdge, Mesh->LatEdge);
      OMEGA_SCOPE(LonEdge, Mesh->LonEdge);
      OMEGA_SCOPE(AngleEdge, Mesh->AngleEdge);
      OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
      OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

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

      OMEGA_SCOPE(LatEdge, Mesh->LatEdge);
      OMEGA_SCOPE(LonEdge, Mesh->LonEdge);
      OMEGA_SCOPE(AngleEdge, Mesh->AngleEdge);
      OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
      OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

      const Clock *ModelClock = StepClock.get();
      R8 ElapsedTimeSec;
      TimeInterval ElapsedTimeInterval = SimTime - ModelClock->getStartTime();
      ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

      const R8 Tau  = 12. * Day2Sec; // 12 days in seconds
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
// NextTracers = (CurTracers * PseudoThickness2(TimeLevel2)) +
//               Coeff * TracersTend) / PseudoThickness1(TimeLevel1)
void TimeStepper::updateTracersByTend(const Array3DReal &NextTracers,
                                      const Array3DReal &CurTracers,
                                      OceanState *State1, int TimeLevel1,
                                      OceanState *State2, int TimeLevel2,
                                      TimeInterval Coeff) const {

   Array2DReal PseudoThick1 = State1->getPseudoThickness(TimeLevel1);
   Array2DReal PseudoThick2 = State2->getPseudoThickness(TimeLevel2);

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
                     (CurTracers(L, ICell, K) * PseudoThick2(ICell, K) +
                      CoeffSeconds * TracerTend(L, ICell, K)) /
                     PseudoThick1(ICell, K);
              });
       });
}

//------------------------------------------------------------------------------
// couple tracer array to pseudo-thickness
void TimeStepper::weightTracers(const Array3DReal &NextTracers,
                                const Array3DReal &CurTracers,
                                OceanState *CurState, int TimeLevel1) const {

   Array2DReal CurThickness = CurState->getPseudoThickness(TimeLevel1);
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

   Array2DReal NextThick = State->getPseudoThickness(TimeLevel);
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

//------------------------------------------------------------------------------
// Compute KPP fields once per time step, before any tendency evaluation.
void TimeStepper::updateKPPFields(OceanState *State, int TracerTimeLevel,
                                  int ThickTimeLevel, int VelTimeLevel) const {

   if (!Tend->VelVertMixTendencyEnable && !Tend->TracerVertMixTendencyEnable)
      return;

   KPPMix *KPPInstance = KPPMix::getInstance();
   if (!KPPInstance || !KPPInstance->Enabled)
      return;

   Array3DReal CurTracerArray = Tracers::getAll(TracerTimeLevel);
   const bool UseTracerForcing =
       Tend->SfcTracerForcing.Enabled || Tend->KPPNonLocalTracerFlux.Enabled;
   KPPInstance->computeKPPFields(State, CurTracerArray, ThickTimeLevel,
                                 VelTimeLevel, UseTracerForcing);
}

//------------------------------------------------------------------------------
// Apply implicit vertical mixing after state/tracer time levels are updated.
void TimeStepper::applyImplicitVerticalMixing(
    OceanState *State, int TracerTimeLevel, int ThickTimeLevel,
    int VelTimeLevel, const std::string &TimerPrefix) const {

   Array3DReal CurTracerArray = Tracers::getAll(TracerTimeLevel);

   VertMix *VMix = VertMix::getInstance();
   if (!VMix)
      return;

   if (VMix->VelVertMixSetup.Enabled or VMix->TracerVertMixSetup.Enabled) {
      const int NTracers = Tracers::getNumTracers();
      VMix->VertMixImplicit(State, AuxState, CurTracerArray, NTracers,
                            VelTimeLevel);

      // Re-exchange halos after vertical mixing
      const MPI_Comm Comm = MeshHalo->getComm();
      Pacer::timingBarrier(TimerPrefix + ":vMixHaloExchBarrier", 3, Comm);
      Pacer::start(TimerPrefix + ":vMixHaloExch", 3);
      State->exchangeHalo(VelTimeLevel);
      Tracers::exchangeHalo(VelTimeLevel);
      Pacer::stop(TimerPrefix + ":vMixHaloExch", 3);
   }
}

void TimeStepper::prepareFirstStep(OceanState *State,
                                   const Array3DReal &CurTracerArray,
                                   bool ReadRestart) {
   prepareNextStep(State, CurTracerArray);

   PrepareFirstStepCalled = true;
}

//------------------------------------------------------------------------------
// Computes all auxiliary variables and performs other setup (e.g.
// prescribes state) for the first stage of the next time step. It also
// validates the current state.
void TimeStepper::prepareNextStep(OceanState *State,
                                  const Array3DReal &CurTracerArray) const {
   const int CurLevel        = 0;
   const TimeInstant SimTime = StepClock->getCurrentTime();

   prescribeState(State, CurLevel, State, CurLevel, SimTime);

   AuxState->computeAll(State, CurTracerArray, CurLevel, CurLevel,
                        FirstProjFactor * TimeStep);

   validateOceanState(State, AuxState, VertCoord::getDefault(), CurLevel);

   updateKPPFields(State, CurLevel, CurLevel, CurLevel);
}

} // namespace OMEGA
