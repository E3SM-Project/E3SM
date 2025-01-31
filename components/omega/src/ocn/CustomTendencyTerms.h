#ifndef OMEGA_CUSTOMTENDENCYTERMS_H
#define OMEGA_CUSTOMTENDENCYTERMS_H
//===-- ocn/CustomTendencyTerms.h - Custom tendency terms -------*- C++ -*-===//
//
/// \file
/// \brief Contains customized tendency terms for the thickness and momentum
///        equations
///
/// For details on the manufactured solution class, see Bishnu et al. (2024)
/// (https://doi.org/10.1029/2022MS003545) and the manufactured solution test
/// case in Polaris. The Polaris package leverages this feature to validate
/// an expected order of convergence of Omega.
//
//===----------------------------------------------------------------------===//

#include "HorzMesh.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"

namespace OMEGA {

//===-----------------------------------------------------------------------===/
// A class for the manufactured solution tendency terms
//===-----------------------------------------------------------------------===/
class ManufacturedSolution {

 public:
   //===--------------------------------------------------------------------===/
   // Manufactured tendency term for the thickness equation
   //===--------------------------------------------------------------------===/
   struct ManufacturedThicknessTendency {

      // Constants defined in 'init'
      TimeInstant ReferenceTime;
      R8 H0;
      R8 Eta0;
      R8 Kx;
      R8 Ky;
      R8 AngFreq;

      void operator()(Array2DReal ThicknessTend, const OceanState *State,
                      const AuxiliaryState *AuxState, int ThickTimeLevel,
                      int VelTimeLevel, TimeInstant Time) const;
   }; // end struct ManufacturedThicknessTendency

   //===--------------------------------------------------------------------===/
   // Manufactured tendency term for the momentum equation
   //===--------------------------------------------------------------------===/
   struct ManufacturedVelocityTendency {

      // Constants defined in 'init'
      TimeInstant ReferenceTime;
      R8 Grav;
      R8 Eta0;
      R8 Kx;
      R8 Ky;
      R8 AngFreq;
      R8 ViscDel2;
      R8 ViscDel4;
      bool VelDiffTendencyEnable;
      bool VelHyperDiffTendencyEnable;

      void operator()(Array2DReal NormalVelTend, const OceanState *State,
                      const AuxiliaryState *AuxState, int ThickTimeLevel,
                      int VelTimeLevel, TimeInstant Time) const;

   }; // end struct ManufacturedVelocityTendency

   // Instances of manufactured tendencies
   ManufacturedThicknessTendency ManufacturedThickTend;
   ManufacturedVelocityTendency ManufacturedVelTend;

   int init();

}; // end class ManufacturedSolution

} // end namespace OMEGA

//===-----------------------------------------------------------------------===/

#endif
