#ifndef OMEGA_SPLIT_EXPLICIT_INIT_H
#define OMEGA_SPLIT_EXPLICIT_INIT_H
//===-- SplitExplicitInit.h - split-explicit initialization ----*- C++ -*-===//
//
/// \file
/// \brief Initialization helpers for split-explicit time stepping.
//
//===----------------------------------------------------------------------===//

#include "HorzMesh.h"
#include "OceanState.h"
#include "SplitExplicitTypes.h"
#include "VertCoord.h"

namespace OMEGA {

class SplitExplicitInit {
 public:
   static SplitExplicitConfig readConfigOptions(
       const TimeInterval &TimeStep, ///< [in] full baroclinic time step
       bool IsUnsplit ///< [in] true for UnsplitRK2 (no barotropic split)
   );

   static void allocateScratch(
       SplitExplicitScratch &Scratch, ///< [inout] split-explicit scratch data
       const HorzMesh *Mesh,          ///< [in] horizontal mesh
       const VertCoord *VCoord,       ///< [in] vertical coordinate
       const std::string &Name        ///< [in] owning time stepper name
   );

   static void
   computeVelocitySplit(OceanState *State,       ///< [inout] ocean state
                        const HorzMesh *Mesh,    ///< [in] horizontal mesh
                        const VertCoord *VCoord, ///< [in] vertical coordinate
                        I4 TimeLevel ///< [in] state time level to split
   );

   static void computeUnsplitVelocitySplit(
       OceanState *State,       ///< [inout] ocean state
       const HorzMesh *Mesh,    ///< [in] horizontal mesh
       const VertCoord *VCoord, ///< [in] vertical coordinate
       I4 CurLevel,             ///< [in] state current level to initialize
       I4 NextLevel             ///< [in] state next level to initialize
   );

   static void initializeBarotropicPressure(
       SplitExplicitScratch &Scratch, ///< [inout] split-explicit scratch data
       OceanState *State,             ///< [inout] ocean state
       const HorzMesh *Mesh,          ///< [in] horizontal mesh
       const VertCoord *VCoord,       ///< [in] vertical coordinate
       I4 TimeLevel                   ///< [in] state time level to initialize
   );

 private:
   static SplitExplicitBarotropicStepperType getBtrTimeStepperFromStr(
       const std::string &InString ///< [in] barotropic stepping method
   );

   static I4 computeSubcycleCount(
       const TimeInterval &TimeStep,   ///< [in] full baroclinic time step
       const TimeInterval &BtrTimeStep ///< [in] requested barotropic time step
   );
};

} // namespace OMEGA

#endif
