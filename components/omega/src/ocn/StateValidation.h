#ifndef OMEGA_STATEVALIDATION_H
#define OMEGA_STATEVALIDATION_H
//===-- ocn/StateValidation.h - ocean state validation ----------*- C++ -*-===//
//
/// \file
/// \brief Declares the validateOceanState function for ocean state validation
///
/// Provides a function that validates the ocean prognostic state and selected
/// auxiliary/tracer fields by checking for NaN values and out-of-bounds
/// conditions. If any check fails the function logs a critical error with a
/// backtrace and aborts via MPI_Abort on the local MPI communicator.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "OceanState.h"
#include "Tracers.h"

namespace OMEGA {

/// Check ocean state fields for NaN values and out-of-bounds conditions.
///
/// The following fields are validated:
///   - LayerThickness      : [1e-10, 1000]  (from OceanState)
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
/// \param[in] TimeLevel   Time level index to validate (typically 0 = current)
void validateOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                        I4 TimeLevel);

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_STATEVALIDATION_H
