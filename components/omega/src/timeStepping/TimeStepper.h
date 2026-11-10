#ifndef OMEGA_TIMESTEPPER_H
#define OMEGA_TIMESTEPPER_H
//===--- TimeStepper.h - time stepper --------------------*- C++ --*-------===//
//
/// \file
/// \brief Contains the base class for all Omega time steppers
///
/// The TimeStepper class defines the interface of a time stepper and contains
/// data and methods common to all time steppers
/// TimeStepper options are set in the TimeIntegration configuration:
/// \ConfigInput
/// # Example time integration configuration (for default config)
/// TimeIntegration:
///    # Calendar to use for integration. Supported options are:
///    # Gregorian, No Leap, Julian, Julian Day, Modified Julian Day,
///    # 360 Day, Custom, No Calendar
///    # These options are described in more detail in the UserGuide
///    # Default is the No Leap (Gregorian calendar with no leap years)
///    CalendarType: No Leap
///    # Algorithm to use for the dynamics time integration
///    # Supported options are Forward-Backward (default), RungeKutta2,
///    # RungeKutta4, SplitExplicitRK2 and UnsplitRK2
///    TimeStepper: Forward-Backward
///    # Time step to use, in form of DDDD_hh:mm:ss (days, hours, minutes, secs)
///    TimeStep: 0000_00:10:00
///    # Start type. Options are StartUp (for starting a simulation from scratch
///    # with an initial state), Continue (for continuing a simulation from a
///    # restart file, Branch (for starting a simulation from a restart file
///    # but possibly changing case name) or Hybrid (for starting a simulation
///    # from restart but changing the simulation start time to a new time)
///    StartType: StartUp
///    # Start time of full simulation (YYYY-MM-DD_hh:mm:ss)
///    StartTime: 0001-01-01_00:00:00
///    # StopType is a string that determines, together with the criterion
///    # below, when to stop the simulation. Supported options are AtTime (to
///    # stop at a specific time instant), AfterDuration (to stop after a
///    # specified time interval - typical of a production simulation when the
///    # simulation is advanced for a time that fits within a queue limit), and
///    # Coupled (to stop based on a signal from the coupled model driver).
///    StopType: AtTime
///    # Stop criterion is a time string that determines, based on the
///    # StopType above, when to stop the simulation. If StopType is
///    # AtTime, this is a string formatted as YYYY-MM-DD_hh:mm:ss that
///    # specifies the time instant at which to stop. If StopType is
///    # AfterDuration, the string is any supported TimeInterval string
///    # (typically DDDD_HH:MM:SS) that represents the duration of the segment
///    # of a simulation. For the Coupled StopType, the criterion string is
///    # ignored and not required.
///    StopCriterion: 0001-01-01_02:00:00
///    # Options shared by the split time steppers SplitExplicitRK2 and
///    # UnsplitRK2. The subgroup is required whenever one of those steppers
///    # is selected. BtrTimeStep is required within it for SplitExplicitRK2
///    # and must be smaller than TimeStep; the remaining keys are optional
///    # and fall back to the SplitExplicitConfig defaults.
///    ModeSplitShare:
///      # BtrTimeStepper and BtrTimeStep are ignored by UnsplitRK2, which
///      # has no barotropic subcycle
///      BtrTimeStepper: Predictor-Corrector
///      BtrTimeStep: 0000_00:00:20
///      # Iterations over stages 1-3. Two are needed for the RK2 predictor-
///      # corrector structure; one degrades the scheme to forward Euler.
///      NTimeStepIteration: 2
///      NBclCoriolisIteration: 2
///      # Recompute the velocity split from NormalVelocity at every slow step
///      ReinitSplitVelocity: false
/// \EndConfigInput
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "StateValidation.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "Tracers.h"
#include "VertCoord.h"

#include <map>
#include <memory>
#include <optional>
#include <string>

namespace OMEGA {

/// An enum for every time stepper type
/// needs to extended every time a new time stepper is added
enum class TimeStepperType {
   ForwardBackward,
   RungeKutta4,
   RungeKutta2,
   SplitExplicitRK2,
   UnsplitRK2,
   Invalid
};

/// An enum for the time stepper start option
enum class TimeStepperStartType {
   StartUp,  // Simulation will start from an initial state at sim start time
   Continue, // Simulation will start from a restart at the restart time
   Branch,   // Identical to continue, though case or sim name can change
   Hybrid,   // Simulation will start from a restart but reset to start time
   Coupled,  // Omega will receive start type from coupler based on E3SM opts
   Invalid   // Invalid or undefined stop option
};

/// An enum for the time stepper stop option
enum class TimeStepperStopType {
   AtTime,        // Simulation will stop at specified time
   AfterDuration, // Simulation will stop after a specified time interval
   Coupled,       // Simulation will stop based on an a signal from coupler
   Invalid        // Invalid or undefined stop option
};

/// An enum describing how a state variable should be prescribed from the
/// reference time level
enum class PrescribeStateType { None, Init, NonDivergent, Divergent, Invalid };

//------------------------------------------------------------------------------
// Utility routines
/// Translate string for time stepper type into enum
TimeStepperType getTimeStepperFromStr(
    const std::string &InString ///< [in] choice of time stepping method
);

/// Translate the input start option string to enum option
TimeStepperStartType getTimeStepperStartTypeFromStr(
    const std::string &InString ///< [in] choice for starting simulation leg
);

/// Translate the input E3SM start option integer to enum start option
TimeStepperStartType getTimeStepperStartTypeFromE3SM(
    const int E3SMOpt ///< [in] E3SM integer start type
);

/// Translate the input stop option string to enum option
TimeStepperStopType getTimeStepperStopTypeFromStr(
    const std::string &InString ///< [in] choice for stopping simulation
);

/// Translate string for prescribe state type into enum
PrescribeStateType getPrescribeStateTypeFromStr(
    const std::string &InString ///< [in] choice of prescribe method
);

//------------------------------------------------------------------------------
/// A base class for Omega time steppers
///
/// The TimeStepper class defines the interface of a time stepper, manages
/// time stepper objects and contains common routines for state updates
class TimeStepper {
 public:
   virtual ~TimeStepper() = default;

   /// The main method that every time stepper needs to define. Advances state
   /// by one time step, from Time to Time + TimeStep
   virtual void doStep(OceanState *State,   ///< [inout] model state
                       TimeInstant &SimTime ///< [inout] current simulation time
   ) const = 0;

   /// Method-specific initialization to set up the first time step after
   /// initial/restart read. Calling doStep before calling this function is an
   /// error.
   virtual void prepareFirstStep(
       OceanState *State, const Array3DReal &CurTracerArray,
       bool ReadRestart ///< [in] true if restart input initialized the state
   );

   /// 1st stage of two-stage initialization for the default time stepper
   /// This version with no input parameters is for standalone mode
   /// with parameters from configuration input
   static void init1();

   /// 1st stage of two-stage initialization for the default time stepper
   /// This version with start parameters is for coupled mode but is also
   /// called by the init1 routine above to complete standalone init
   static void
   init1(TimeStepperStartType InStartType, ///< [in] start choice for this leg
         const TimeInstant &StartTime);    ///< [in] start time for full sim

   /// 2nd phase of Initialization for the default time stepper
   static void init2();

   /// Create a time stepper when all components are known
   /// Note that if StopType is AfterDuration, the StopTime and EndAlarm
   /// are computed based on the clock's current time and must be reset if
   /// the current time is reset (eg by reading a restart) using resetEndAlarm
   static TimeStepper *
   create(const std::string &InName,      ///< [in] name of time stepper
          TimeStepperType InType,         ///< [in] type (time stepping method)
          Tendencies *InTend,             ///< [in] ptr to tendencies
          AuxiliaryState *InAuxState,     ///< [in] ptr to aux state variables
          HorzMesh *InMesh,               ///< [in] ptr to mesh information
          VertCoord *InVCoord,            ///< [in] ptr to vertical coordinate
          Halo *InMeshHalo,               ///< [in] ptr to halos
          const TimeInterval &InTimeStep, ///< [in] time step
          const TimeStepperStartType InStartType, ///< [in] option for start
          const TimeInstant &InStartTime, ///< [in] full simulation start time
          const TimeStepperStopType InStopType, ///< [in] option to stop
          std::optional<TimeInstant> InStopTime = std::nullopt,
          ///< [in] stop time if StopType is AtTime
          std::optional<TimeInterval> InDuration = std::nullopt
          ///< [in] duration of current run segment if StopType is AfterDuration
   );

   /// Create a time stepper when time information is needed before state
   /// and tendencies are defined. It creates an instance and only fills
   /// the time information. Data pointers are attached later.
   /// Note that if StopType is AfterDuration, the StopTime and EndAlarm
   /// are computed based on the clock's current time and must be reset if
   /// the current time is reset (eg by reading a restart) using resetEndAlarm
   static TimeStepper *
   create(const std::string &InName,              ///< [in] name of time stepper
          TimeStepperType InType,                 ///< [in] type (method)
          const TimeInterval &InTimeStep,         ///< [in] time step
          const TimeStepperStartType InStartType, ///< [in] option for starting
          const TimeInstant &InStartTime,         ///< [in] full sim start time
          const TimeStepperStopType InStopType,   ///< [in] option for stopping
          std::optional<TimeInstant> InStopTime = std::nullopt,
          ///< [in] stop time if StopType is AtTime
          std::optional<TimeInterval> InDuration = std::nullopt
          ///< [in] duration of current run segment if StopType is AfterDuration
   );

   /// If the current segment start time is reset (eg by reading a restart), the
   /// relevant stop times and alarms must be reset.  This only modifies
   /// StopTime and EndAlarm if StopType is AfterDuration
   void resetEndAlarm();

   /// For 2-step creation, this attaches all the data pointers to an instance
   /// once the data and tendencies have been created.
   void attachData(
       Tendencies *InTend,         ///< [in] ptr to tendencies
       AuxiliaryState *InAuxState, ///< [in] ptr to needed aux state variables
       HorzMesh *InMesh,           ///< [in] ptr to mesh information
       VertCoord *InVCoord,        ///< [in] ptr to vertical coordinate
       Halo *InMeshHalo            ///< [in] ptr to halos
   );

   // Delete/destroy functions
   /// Remove time stepper by name
   static void erase(const std::string &Name);

   /// Remove all time steppers
   static void clear();

   // Retrieval functions
   /// Get the default time stepper
   static TimeStepper *getDefault();

   /// Get time stepper by name
   static TimeStepper *
   get(const std::string &Name ///< [in] name of stepper to retrieve
   );

   /// Get name of time stepper from instance
   std::string getName() const;

   /// Get type (enum) of time stepper from instance
   TimeStepperType getType() const;

   /// Get number of time levels from instance
   int getNTimeLevels() const;

   /// Get time step
   TimeInterval getTimeStep() const;

   /// Get start option
   TimeStepperStartType getStartType() const;

   /// Get start time
   TimeInstant getStartTime() const;

   /// Get stop option
   TimeStepperStopType getStopType() const;

   /// Get stop time
   TimeInstant getStopTime() const;

   /// Get duration of run segment
   TimeInterval getDuration() const;

   /// Get a pointer to the clock
   Clock *getClock();

   /// Get a pointer to the end alarm
   Alarm *getEndAlarm();

   /// Change time step
   void changeTimeStep(const TimeInterval &TimeStepIn ///< [in] new time step
   );

   /// Get number of doStep calls made on this instance
   I8 getStepCount() const;

   /// Is this a split time stepper. False by default
   virtual bool isSplit() const;

   // these should be protected, they are public only because of CUDA
   // limitations

   /// Updates state using tendency terms
   /// State1(TimeLevel1) = State2(TimeLevel2) + Coeff * Tend
   void updateStateByTend(
       OceanState *State1, ///< [out] updated state
       int TimeLevel1,     ///< [in] time level index for new time
       OceanState *State2, ///< [in] state data for current time
       int TimeLevel2,     ///< [in] time level index for current time
       TimeInterval Coeff  ///< [in] time-related coeff for tendency
   ) const;

   /// Updates pseudo-thickness using tendency terms
   /// PseudoThickness1(TimeLevel1) = PseudoThickness2(TimeLevel2) +
   ///                               Coeff * PseudoThicknessTend
   void updateThicknessByTend(
       OceanState *State1, ///< [out] updated pseudo-thickness in state
       int TimeLevel1,     ///< [in] time level index for new time
       OceanState *State2, ///< [in] state (thickness) for current time
       int TimeLevel2,     ///< [in] time level index for current time
       TimeInterval Coeff  ///< [in] time-related coeff for tendency
   ) const;

   /// Updates velocity using tendency terms
   /// NormalVelocity1(TimeLevel1) = NormalVelocity2(TimeLevel2) +
   ///                               Coeff * NormalVelocityTend
   void updateVelocityByTend(
       OceanState *State1, ///< [out] updated state (velocity)
       int TimeLevel1,     ///< [in] time level index for new time
       OceanState *State2, ///< [in] state (velocity) for current time
       int TimeLevel2,     ///< [in] time level index for current time
       TimeInterval Coeff  ///< [in] time-related coeff for tendency
   ) const;

   /// Resets the pseudo-thickness at the working time level to the initial
   /// condition stored at the reference time level.
   void prescribeThickness(
       OceanState *State1, ///< [out] destination state
       int TimeLevel1,     ///< [in] time level index of the reference data
       OceanState *State2, ///< [in] source state (initial condition)
       int TimeLevel2      ///< [in] time level index of the destination data
   ) const;

   /// Resets the velocity at the working time level to the initial condition.
   void prescribeVelocity(
       OceanState *State1, ///< [out] destination state
       int TimeLevel1,     ///< [in] time level index of the reference data
       OceanState *State2, ///< [in] source state (initial condition)
       int TimeLevel2,     ///< [in] time level index of the destination data
       const TimeInstant &SimTime ///< [in] current simulation time
   ) const;

   /// Reset thickness and velocity to their initial values
   void prescribeState(
       OceanState *State1, ///< [out] destination state
       int TimeLevel1,     ///< [in] time level index of the reference data
       OceanState *State2, ///< [in] source state (initial condition)
       int TimeLevel2,     ///< [in] time level index of the destination data
       const TimeInstant &SimTime ///< [in] current simulation time
   ) const;

   /// Updates tracers
   /// NextTracers = (CurTracers * PseudoThickness2(TimeLevel2)) +
   ///               Coeff * TracersTend) / PseudoThickness1(TimeLevel1)
   void updateTracersByTend(
       const Array3DReal &NextTracers, ///< [out] updated tracers
       const Array3DReal &CurTracers,  ///< [in]  current tracers
       OceanState *State1, ///< [in] state (thickness) at updated time
       int TimeLevel1,     ///< [in] time level index for new time
       OceanState *State2, ///< [in] state (thickness) for current time
       int TimeLevel2,     ///< [in] time level index for current time
       TimeInterval Coeff  ///< [in] time-related coeff for tendency
   ) const;

   /// couple tracer array to pseudo-thickness
   void weightTracers(
       const Array3DReal &NextTracers, ///< [inout] tracers to modify
       const Array3DReal &CurTracers,  ///< [inout] tracers at current time
       OceanState *CurState,           ///< [in] state (thick) at current time
       int TimeLevel1                  ///< [in] time index
   ) const;

   /// accumulate contributions to the tracer array at the next time level from
   /// each Runge-Kutta stage
   void accumulateTracersUpdate(
       const Array3DReal &AccumTracer, ///< [inout] accumulated tracers
       TimeInterval Coeff ///< [in] time-related coeff for accumulation
   ) const;

   /// normalize tracer array so final array stores concentrations
   void finalizeTracersUpdate(
       const Array3DReal &NextTracers, ///< [inout] tracers to normalize
       OceanState *State,              ///< [in] state with thickness data
       int TimeLevel                   ///< [in] time level index
   ) const;

   /// Compute KPP boundary layer depth, mixing coefficients and non-local
   /// flux once for the current time step. Must be called before the first
   /// tendency evaluation of the step so that every stage and the end-of-step
   /// implicit mixing share the same KPP fields.
   void
   updateKPPFields(OceanState *State,   ///< [in] model state
                   int TracerTimeLevel, ///< [in] tracer time level
                   int ThickTimeLevel,  ///< [in] pseudo-thickness time level
                   int VelTimeLevel     ///< [in] velocity time level
   ) const;

   /// Apply implicit vertical mixing after state/tracer time levels are
   /// updated, using the KPP fields computed at the start of the step.
   void applyImplicitVerticalMixing(
       OceanState *State,             ///< [inout] model state
       int TracerTimeLevel,           ///< [in] tracer time level
       int ThickTimeLevel,            ///< [in] pseudo-thickness time level
       int VelTimeLevel,              ///< [in] velocity time level
       const std::string &TimerPrefix ///< [in] timer name prefix
   ) const;

   /// Computes all auxiliary variables and performs other setup (e.g.
   /// prescribes state ) for the first stage of the next time step. It also
   /// validates the current state.
   void prepareNextStep(OceanState *State,
                        const Array3DReal &CurTracerArray) const;

 protected:
   /// Name of time stepper
   std::string Name;

   /// Type of time stepper
   TimeStepperType Type;

   /// Number of time levels required
   int NTimeLevels;

   /// Time step
   TimeInterval TimeStep;

   /// Start option
   TimeStepperStartType StartType;

   /// Start time
   TimeInstant StartTime;

   /// Stop option
   TimeStepperStopType StopType;

   /// Stop time for this run segment
   TimeInstant StopTime;

   /// Duration of this run segment
   TimeInterval Duration;

   /// Alarm that rings at StopTime if needed
   std::unique_ptr<Alarm> EndAlarm;

   /// Clock for this time stepper
   /// For the default time stepper, this is the model clock
   std::unique_ptr<Clock> StepClock;

   /// Number of doStep calls made on this instance since creation
   mutable I8 StepCount = 0;

   /// Prescribe state configuration
   PrescribeStateType PrescribeThicknessMode;
   PrescribeStateType PrescribeVelocityMode;

   // Pointers to objects needed by every time stepper
   Tendencies *Tend;         /// Ptr to tendency terms
   AuxiliaryState *AuxState; /// Ptr to auxiliary state data
   HorzMesh *Mesh;           /// Ptr to horizontal mesh info
   VertCoord *VCoord;        /// Ptr to vertical coordinate
   Halo *MeshHalo;           /// Ptr to defined halos

   // Flag that indicates prepareFirstStep has been called
   bool PrepareFirstStepCalled = false;

   // Time step projection factor of the first time stepper stage
   Real FirstProjFactor;

   /// Reset frazil fluxes accumulated across the current outer ocean step.
   void resetFrazilFluxes() const;

   /// Function for any method-specific modifications for the default stepper
   virtual void finalizeInit() {}

   /// Constructor creates a new instance and fills in most of the time
   /// related data. The stop time is added later if needed and the attachData
   /// function is used to add the data pointers
   TimeStepper(
       const std::string &InName,      ///< [in] name of time stepper
       TimeStepperType InType,         ///< [in] type (time stepping method)
       I4 InNTimeLevels,               ///< [in] num time levels for method
       const TimeInterval &InTimeStep, ///< [in] time step
       const TimeStepperStartType InStartType, ///< [in] option for starting
       const TimeInstant &InStartTime,         ///< [in] start time for full sim
       const TimeStepperStopType InStopType,   ///< [in] option for stopping
       std::optional<TimeInstant> InStopTime = std::nullopt,
       ///< [in] stop time if StopType is AtTime
       std::optional<TimeInterval> InDuration = std::nullopt
       ///< [in] duration of current run segment if StopType is AfterDuration
   );

   // Disable copy constructor
   TimeStepper(const TimeStepper &) = delete;
   TimeStepper(TimeStepper &&)      = delete;

 private:
   /// Default model time stepper
   static TimeStepper *DefaultTimeStepper;
   /// All defined time steppers
   static std::map<std::string, std::unique_ptr<TimeStepper>> AllTimeSteppers;
   /// Prescribe modes parsed from config
   static PrescribeStateType DefaultPrescribeThicknessMode;
   static PrescribeStateType DefaultPrescribeVelocityMode;
};

} // namespace OMEGA
#endif
