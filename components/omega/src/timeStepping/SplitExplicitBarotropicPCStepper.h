#ifndef OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
#define OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
//===-- SplitExplicitBarotropicPCStepper.h - SE stage 2 ------*- C++ -*-===//
//
/// \file
/// \brief Forward-backward predictor-corrector barotropic subcycle framework.
//
//===----------------------------------------------------------------------===//

#include "Halo.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "SplitExplicitTypes.h"
#include "TimeMgr.h"
#include "VertCoord.h"

namespace OMEGA {

class AuxiliaryState;

class SplitExplicitBarotropicPCStepper {
 public:
   void
   init(const AuxiliaryState *InAuxState, ///< [in] provisional auxiliary state
        SplitExplicitScratch *InScratch,  ///< [inout] split-explicit scratch
        const SplitExplicitConfig *InConfig, ///< [in] split-explicit options
        const HorzMesh *InMesh,              ///< [in] horizontal mesh
        Halo *InMeshHalo,                    ///< [in] mesh halo exchange
        const VertCoord *InVCoord            ///< [in] vertical coordinate
   );

   void doBarotropicVelocityUpdate(
       OceanState *State,                ///< [inout] model state
       I4 CurLevel,                      ///< [in] state time level to update
       I4 NextLevel,                     ///< [in] state time level to update
       const TimeInterval &StageTimeStep ///< [in] current stage time step
   ) const;

 private:
   const AuxiliaryState *AuxState      = nullptr;
   SplitExplicitScratch *Scratch       = nullptr;
   const SplitExplicitConfig *SEConfig = nullptr;
   const HorzMesh *Mesh                = nullptr;
   Halo *MeshHalo                      = nullptr;
   const VertCoord *VCoord             = nullptr;
};

} // namespace OMEGA

#endif
