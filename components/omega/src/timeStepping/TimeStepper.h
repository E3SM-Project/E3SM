#ifndef OMEGA_TIMESTEPPER_H
#define OMEGA_TIMESTEPPER_H
//===--- TimeStepper.h - time stepper --------------------*- C++ --*-------===//
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
#include "Tendencies.h"
#include "TimeMgr.h"
#include "Tracers.h"

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

/// An enum for every time stepper type
/// needs to extended every time a new time stepper is added
enum class TimeStepperType {
   ForwardBackward,
   RungeKutta4,
   RungeKutta2,
   Invalid
};

//------------------------------------------------------------------------------
// Utility routine
/// Translate string for time stepper type into enum
TimeStepperType getTimeStepperFromStr(
    const std::string &InString ///< [in] choice of time stepping method
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

   /// 1st phase of Initialization for the default time stepper
   static int init1();

   /// 2nd phase of Initialization for the default time stepper
   static int init2();

   /// Create a time stepper when all components are known
   static TimeStepper *
   create(const std::string &InName,      ///< [in] name of time stepper
          TimeStepperType InType,         ///< [in] type (time stepping method)
          const TimeInstant &InStartTime, ///< [in] start time for time stepping
          const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
          const TimeInterval &InTimeStep, ///< [in] time step
          Tendencies *InTend,             ///< [in] ptr to tendencies
          AuxiliaryState *InAuxState,     ///< [in] ptr to aux state variables
          HorzMesh *InMesh,               ///< [in] ptr to mesh information
          Halo *InMeshHalo                ///< [in] ptr to halos
   );

   /// Create a time stepper when time information is needed before state
   /// and tendencies are defined. It creates an instance and only fills
   /// the time information. Data pointers are attached later.
   static TimeStepper *
   create(const std::string &InName, ///< [in] name of time stepper
          TimeStepperType InType,    ///< [in] type (time stepping method)
          TimeInstant &InStartTime,  ///< [in] start time for time stepping
          TimeInstant &InStopTime,   ///< [in] stop  time for time stepping
          TimeInterval &InTimeStep   ///< [in] time step
   );

   /// For 2-step creation, this attaches all the data pointers to an instance
   /// once the data and tendencies have been created.
   void attachData(
       Tendencies *InTend,         ///< [in] ptr to tendencies
       AuxiliaryState *InAuxState, ///< [in] ptr to needed aux state variables
       HorzMesh *InMesh,           ///< [in] ptr to mesh information
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

   /// Get start time
   TimeInstant getStartTime() const;

   /// Get stop time
   TimeInstant getStopTime() const;

   /// Get a pointer to the clock
   Clock *getClock();

   /// Get a pointer to the end alarm
   Alarm *getEndAlarm();

   /// Change time step
   void changeTimeStep(const TimeInterval &TimeStepIn ///< [in] new time step
   );

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

   /// Updates layer thickness using tendency terms
   /// LayerThickness1(TimeLevel1) = LayerThickness2(TimeLevel2) +
   ///                               Coeff * LayerThicknessTend
   void updateThicknessByTend(
       OceanState *State1, ///< [out] updated layer thickness in state
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

   /// Updates tracers
   /// NextTracers = (CurTracers * LayerThickness2(TimeLevel2)) +
   ///               Coeff * TracersTend) / LayerThickness1(TimeLevel1)
   void updateTracersByTend(
       const Array3DReal &NextTracers, ///< [out] updated tracers
       const Array3DReal &CurTracers,  ///< [in]  current tracers
       OceanState *State1, ///< [in] state (thickness) at updated time
       int TimeLevel1,     ///< [in] time level index for new time
       OceanState *State2, ///< [in] state (thickness) for current time
       int TimeLevel2,     ///< [in] time level index for current time
       TimeInterval Coeff  ///< [in] time-related coeff for tendency
   ) const;

   /// couple tracer array to layer thickness
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

 protected:
   /// Name of time stepper
   std::string Name;

   /// Type of time stepper
   TimeStepperType Type;

   /// Number of time levels required
   int NTimeLevels;

   /// Time step
   TimeInterval TimeStep;

   /// Start time
   TimeInstant StartTime;

   /// Stop time
   TimeInstant StopTime;

   /// Alarm that rings at StopTime
   std::unique_ptr<Alarm> EndAlarm;

   /// Clock for this time stepper
   /// For the default time stepper, this is the model clock
   std::unique_ptr<Clock> StepClock;

   // Pointers to objects needed by every time stepper
   Tendencies *Tend;         /// Ptr to tendency terms
   AuxiliaryState *AuxState; /// Ptr to auxiliary state data
   HorzMesh *Mesh;           /// Ptr to horizontal mesh info
   Halo *MeshHalo;           /// Ptr to defined halos

   /// Function for any method-specific modifications for the default stepper
   virtual void finalizeInit() {}

   /// Constructor creates a new instance and fills in the time
   /// related data. attachData function is used to add the data pointers
   TimeStepper(
       const std::string &InName,      ///< [in] name of time stepper
       TimeStepperType InType,         ///< [in] type (time stepping method)
       I4 InNTimeLevels,               ///< [in] num time levels for method
       const TimeInstant &InStartTime, ///< [in] start time for time stepping
       const TimeInstant &InStopTime,  ///< [in] stop  time for time stepping
       const TimeInterval &InTimeStep  ///< [in] time step
   );

   // Disable copy constructor
   TimeStepper(const TimeStepper &) = delete;
   TimeStepper(TimeStepper &&)      = delete;

 private:
   /// Default model time stepper
   static TimeStepper *DefaultTimeStepper;
   /// All defined time steppers
   static std::map<std::string, std::unique_ptr<TimeStepper>> AllTimeSteppers;
};

} // namespace OMEGA
#endif
