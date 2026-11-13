#ifndef OMEGA_STATEVALIDATION_H
#define OMEGA_STATEVALIDATION_H
//===-- ocn/StateValidation.h - ocean state validation ----------*- C++ -*-===//
//
/// \file
/// \brief Declares state validation functions for ocean state validation
///
/// Provides two functions:
///   - checkOceanState: checks for NaN and out-of-bounds conditions and
///     returns the NaN count and out-of-bounds count as a pair without
///     aborting. Suitable for testing.
///   - validateOceanState: calls checkOceanState and conditionally aborts
///     via MPI_Abort based on runtime YAML options AbortOnNan and AbortOnOOB
///     (under Omega::State::Validation). When both flags are false a
///     CRITICAL log message is emitted and the run continues.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "OceanState.h"
#include "Tracers.h"
#include "VertCoord.h"

namespace OMEGA {

/// Check ocean state fields for NaN values and out-of-bounds conditions
/// without aborting, returning counts of each error type found.
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
/// \return std::pair<I4,I4> {NaN count, out-of-bounds count}; both 0 if valid
std::pair<I4, I4> checkOceanState(const OceanState *State,
                                  const AuxiliaryState *AuxState,
                                  const VertCoord *VCoord, I4 TimeLevel);

/// Check ocean state fields for NaN values and out-of-bounds conditions.
///
/// Calls checkOceanState and conditionally aborts via MPI_Abort based on the
/// runtime YAML options AbortOnNan and AbortOnOOB (under
/// Omega::State::Validation). Only active ocean cells (where CellMask > 0)
/// are checked.
///
/// The following fields are validated:
///   - PseudoThickness     : [1e-10, 1000]  (from OceanState)
///   - KineticEnergyCell   : [0, 10]
///                           (from AuxiliaryState::KineticAux)
///   - Temperature tracer  : [-10, 50]      (from Tracers)
///   - Salinity tracer     : [-2, 60]       (from Tracers)
///
/// If AbortOnNan is true and NaN values are detected, or AbortOnOOB is true
/// and out-of-bounds values are detected, a critical error is logged with a
/// stack backtrace and the run is aborted via MPI_Abort. When both flags are
/// false but errors are found, a CRITICAL log message is emitted and execution
/// continues.
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
