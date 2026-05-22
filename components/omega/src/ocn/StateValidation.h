#ifndef OMEGA_STATEVALIDATION_H
#define OMEGA_STATEVALIDATION_H
//===-- ocn/StateValidation.h - ocean state validation ----------*- C++ -*-===//
//
/// \file
/// \brief Declares state validation functions for ocean state validation
///
/// Provides two functions:
///   - checkOceanState: checks for NaN and out-of-bounds conditions and
///     returns the total error count without aborting. Suitable for testing.
///   - validateOceanState: calls checkOceanState and aborts via MPI_Abort if
///     any errors are found.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "OceanState.h"
#include "Tracers.h"
#include "VertCoord.h"

namespace OMEGA {

/// Check ocean state fields for NaN values and out-of-bounds conditions
/// without aborting, returning the total count of errors found.
///
/// Only active ocean cells (where CellMask > 0) are checked. Critical log
/// messages are emitted for each type of error found.
///
/// The following fields are validated:
///   - PseudoThickness     : [1e-10, 1000]  (from OceanState)
///   - KineticEnergyCell   : [0, 10]
///                           (from AuxiliaryState::KineticAux)
///   - Temperature tracer  : [-10, 50]      (from Tracers)
///   - Salinity tracer     : [-2, 60]       (from Tracers)
///
/// \param[in] State       Ocean state to validate
/// \param[in] AuxState    Auxiliary state containing KineticEnergyCell
/// \param[in] VCoord      Vertical coordinate containing the CellMask
/// \param[in] TimeLevel   Time level index to validate (typically 0 = current)
/// \return I4 total count of errors found across all checked fields (0 = valid)
I4 checkOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                   const VertCoord *VCoord, I4 TimeLevel);

/// Check ocean state fields for NaN values and out-of-bounds conditions.
///
/// Calls checkOceanState and aborts via MPI_Abort if any errors are found.
/// Only active ocean cells (where CellMask > 0) are checked.
///
/// The following fields are validated:
///   - PseudoThickness     : [1e-10, 1000]  (from OceanState)
///   - KineticEnergyCell   : [0, 10]
///                           (from AuxiliaryState::KineticAux)
///   - Temperature tracer  : [-10, 50]      (from Tracers)
///   - Salinity tracer     : [-2, 60]       (from Tracers)
///
/// If any check fails a critical error is logged with an informative message
/// and a stack backtrace, and the run is aborted via MPI_Abort on the
/// communicator obtained from the default MachEnv.
///
/// \param[in] State       Ocean state to validate
/// \param[in] AuxState    Auxiliary state containing KineticEnergyCell
/// \param[in] VCoord      Vertical coordinate containing the CellMask
/// \param[in] TimeLevel   Time level index to validate (typically 0 = current)
void validateOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                        const VertCoord *VCoord, I4 TimeLevel);

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_STATEVALIDATION_H
