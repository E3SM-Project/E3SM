//===-- ocn/OceanFinalize.cpp -----------------------------------*- C++ -*-===//
//
// The ocnFinalize method writes a restart file if necessary, and then cleans
// up all Omega objects
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Decomp.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"

namespace OMEGA {

int ocnFinalize(const TimeInstant &CurrTime ///< [in] current sim time
) {

   // error code
   I4 RetVal = 0;

   // Write restart file if necessary

   // clean up all objects
   Tracers::clear();
   TimeStepper::clear();
   Tendencies::clear();
   AuxiliaryState::clear();
   OceanState::clear();
   Dimension::clear();
   Field::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();

   return RetVal;
} // end ocnFinalize

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
