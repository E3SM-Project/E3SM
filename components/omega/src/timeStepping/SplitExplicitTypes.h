#ifndef OMEGA_SPLIT_EXPLICIT_TYPES_H
#define OMEGA_SPLIT_EXPLICIT_TYPES_H
//===-- SplitExplicitTypes.h - split-explicit shared types -----*- C++ -*-===//
//
/// \file
/// \brief Contains shared configuration and scratch arrays for SE time
/// stepping.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "TimeMgr.h"

namespace OMEGA {

enum class SplitExplicitBarotropicStepperType { PredictorCorrector, Invalid };

struct SplitExplicitConfig {
   TimeInterval BtrTimeStep;
   SplitExplicitBarotropicStepperType BtrTimeStepper =
       SplitExplicitBarotropicStepperType::PredictorCorrector;
   I4 NBtrSubcycles = 1;
   // Two iterations are required for the RK2 predictor-corrector structure;
   // a single iteration leaves only the predictor and degrades the scheme to
   // forward Euler.
   I4 NTimeStepIteration    = 2;
   I4 NBclCoriolisIteration = 2;
   bool ReinitSplitVelocity = false;
   Real SplitFactor         = 1._Real;
};

struct SplitExplicitScratch {
   // The barotropic subcycle keeps three buffers per field:
   //    - the state at the start of the subcycle (Cur)
   //    - the predictor output (Pre)
   //    - the corrector output (Cor).
   Array1DReal NormalBarotropicVelocitySubcycleCur;
   Array1DReal NormalBarotropicVelocitySubcyclePre;
   Array1DReal NormalBarotropicVelocitySubcycleCor;
   Array1DReal BarotropicPressureAnomalySubcycleCur;
   Array1DReal BarotropicPressureAnomalySubcyclePre;
   Array1DReal BarotropicPressureAnomalySubcycleCor;
   Array1DReal BarotropicForcing;
   Array1DReal BarotropicFlux;
   Array1DReal BaroclinicPseudoThicknessEdge;
   Array2DReal IterVelocityTend;
   Array2DReal TransportVelocityAdd;
};

} // namespace OMEGA

#endif
