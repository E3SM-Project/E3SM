#ifndef SCREAM_MODEL_NAME_HPP
#define SCREAM_MODEL_NAME_HPP

#include <string>

namespace scream {

// Default model name for EAMxx output files
// NOTE: In the future, this could be provided by the MCT coupling layer
// via a component-specific identifier. For now, hardcoded to distinguish
// EAMxx output from EAM output in E3SM archives.
constexpr const char* EAMXX_MODEL_NAME_DEFAULT = "eamxx";

/**
 * Get the model name to use in output filenames.
 *
 * Currently returns the hardcoded EAMXX_MODEL_NAME_DEFAULT ("eamxx").
 *
 * FUTURE: When E3SM's MCT coupling layer provides component-specific names
 * (rather than generic "ATM"), this function should be updated to:
 * 1. Query the coupling layer for the model-specific name
 * 2. Fall back to EAMXX_MODEL_NAME_DEFAULT if not available
 *
 * HOW TO UPDATE WHEN MCT PROVIDES MODEL NAME:
 *
 * 1. The MCT coupling layer (driver-mct) would need to provide a function
 *    to query the model-specific name (e.g., "eamxx", "eam", "cam").
 *    Currently, seq_comm_name() only returns "ATM" (generic component type).
 *
 * 2. Once available, update get_eamxx_model_name() to call the MCT function.
 *    This might look like:
 *    - Add C-Fortran interface in eamxx_f2c_mod.F90
 *    - Call seq_comm_get_model_name() or similar from MCT
 *    - Return the result if non-empty, else use default
 *
 * 3. The change will be isolated to this single function - no other code
 *    needs modification.
 *
 * 4. Ensure backward compatibility: if MCT name is unavailable or empty,
 *    always fall back to EAMXX_MODEL_NAME_DEFAULT.
 *
 * Example future implementation:
 *   std::string model_name = get_mct_model_name();  // hypothetical MCT function
 *   return model_name.empty() ? EAMXX_MODEL_NAME_DEFAULT : model_name;
 */
inline std::string get_eamxx_model_name() {
  // TODO: Check for MCT-provided model name when available
  // For now, return hardcoded name
  return std::string(EAMXX_MODEL_NAME_DEFAULT);
}

} // namespace scream

#endif // SCREAM_MODEL_NAME_HPP
