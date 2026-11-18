#ifndef OMEGA_TS_SPLIT_EXPLICIT_RK2_H
#define OMEGA_TS_SPLIT_EXPLICIT_RK2_H
//===-- SplitExplicitRK2Stepper.h - split-explicit RK2 step ---*- C++ -*-===//
//
/// \file
/// \brief Contains the framework for split-explicit RK2 time stepping.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitBarotropicPCStepper.h"
#include "SplitExplicitTypes.h"
#include "TimeStepper.h"

namespace OMEGA {

class SplitExplicitRK2Stepper : public TimeStepper {
 public:
   SplitExplicitRK2Stepper(
       const std::string &InName,              ///< [in] name of time stepper
       TimeStepperType InType,                 ///< [in] SplitExplicit type
       const TimeInterval &InTimeStep,         ///< [in] time step
       const TimeStepperStartType InStartType, ///< [in] option for starting sim
       const TimeInstant &InStartTime,         ///< [in] start time for full sim
       const TimeStepperStopType InStopType,   ///< [in] option for stopping
       std::optional<TimeInstant> InStopTime,  ///< [in] stop time if AtTime
       std::optional<TimeInterval> InDuration  ///< [in] duration if AfterDur
   );

   /// Indicate that this is a split time stepper
   bool isSplit() const override;

   void prepareFirstStep(OceanState *State, const Array3DReal &CurTracerArray,
                         bool ReadRestart) override;

   /// Advance the state by one split-explicit RK2 step.
   void doStep(OceanState *State,   ///< [inout] model state
               TimeInstant &SimTime ///< [inout] current simulation time
   ) const override;

   /// Initialize split-explicit state after initial/restart input is read.
   void initializeStateFromInput(
       OceanState *State, ///< [inout] model state after input has been read
       bool ReadRestart   ///< [in] true if restart input initialized the state
   ) const;

   void
   computeTransportVelocity(OceanState *State, ///< [inout] model state
                            I4 TimeLevel ///< [in] time level for split velocity
   ) const;

   /// Removes the barotropic forcing from the iterated baroclinic tendency and
   /// advances the baroclinic velocity to the stage midpoint, in one pass.
   void updateBaroclinicVelocityWithBarotropicForcing(
       OceanState *State,                ///< [inout] model state
       I4 CurLevel,                      ///< [in] current time level
       I4 NextLevel,                     ///< [in] next time level
       const TimeInterval &StageTimeStep ///< [in] current stage time step
   ) const;

   void updateTracersToMidpoint(
       const Array3DReal &MidpointTracers, ///< [out] midpoint tracers
       const Array3DReal &CurTracers,      ///< [in] current tracers
       OceanState *State,                  ///< [in] model state
       I4 CurLevel,                        ///< [in] current state time level
       const TimeInterval &StageTimeStep   ///< [in] full stage time step
   ) const;

   void initializeNextState(OceanState *State, ///< [inout] model state
                            I4 CurLevel,       ///< [in] current time level
                            I4 NextLevel,      ///< [in] next time level
                            Real SplitFactor,  ///< [in] split factor
                            bool ReinitSplitVelocity ///< [in] whether to
                                                     ///< recompute the split
   ) const;

   /// Extrapolates NormalVelocity from the stage midpoint to the new time
   /// level. Only needed on the final time-step iteration.
   void reconstructNormalVelocity(OceanState *State, ///< [inout] model state
                                  I4 CurLevel, ///< [in] current time level
                                  I4 NextLevel ///< [in] next time level
   ) const;

   void finalizeTimeStepIterationState(
       OceanState *State,  ///< [inout] model state
       I4 CurLevel,        ///< [in] current time level
       I4 NextLevel,       ///< [in] next time level
       bool FinalIteration ///< [in] true on the final time-step iteration
   ) const;

 protected:
   /// Performs additional initialization for split-explicit scratch fields.
   void finalizeInit() override;

 private:
   void doBaroclinicVelocityUpdate(
       OceanState *State,                      ///< [inout] model state
       const Array3DReal &TendencyTracerArray, ///< [in] tracers for tendencies
       I4 CurLevel,                            ///< [in] current time level
       I4 NextLevel,                           ///< [in] next time level
       const TimeInstant &StageTime,           ///< [in] current stage time
       const TimeInterval &StageTimeStep,      ///< [in] current stage time step
       I4 Iter                                 ///< [in] current outer iteration
   ) const;

   void doThicknessTracerUpdate(
       OceanState *State,                  ///< [inout] model state
       const Array3DReal &CurTracerArray,  ///< [in] current tracers
       const Array3DReal &NextTracerArray, ///< [out] next tracers
       I4 CurLevel,                        ///< [in] current time level
       I4 NextLevel,                       ///< [in] next time level
       const TimeInstant &StageTime,       ///< [in] current stage time
       const TimeInterval &StageTimeStep,  ///< [in] current stage time step
       bool FinalIteration ///< [in] true on the final time-step iteration
   ) const;

   void doBaroclinicCoriolisIteration(
       OceanState *State,                ///< [inout] model state
       I4 CurLevel,                      ///< [in] current time level
       I4 NextLevel,                     ///< [in] next time level
       const TimeInterval &StageTimeStep ///< [in] current stage time step
   ) const;

   SplitExplicitConfig SEConfig;
   mutable SplitExplicitScratch SEScratch;
   SplitExplicitBarotropicPCStepper BarotropicPCStepper;
};

} // namespace OMEGA

#endif
