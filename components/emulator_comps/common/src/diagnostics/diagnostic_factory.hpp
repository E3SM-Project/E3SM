/**
 * @file diagnostic_factory.hpp
 * @brief Factory functions for creating diagnostics.
 */

#ifndef EMULATOR_DIAGNOSTIC_FACTORY_HPP
#define EMULATOR_DIAGNOSTIC_FACTORY_HPP

#include "derived_diagnostic.hpp"
#include <memory>
#include <mpi.h>

namespace emulator {

/**
 * @brief Metadata for creating diagnostics.
 */
struct DiagnosticMetadata {
  std::vector<double> area_weights; ///< Area weights for horiz averaging
  MPI_Comm comm = MPI_COMM_WORLD;   ///< MPI communicator
  int nlevs = 1;                    ///< Default number of levels
};

/**
 * @brief Parse a field name and create appropriate diagnostic if needed.
 *
 * Recognized patterns:
 * - "{field}_horiz_avg" or "{field}_global_mean" → HorizAvgDiagnostic
 * - "{field}_at_lev{N}" → VertSliceDiagnostic at level N
 *
 * @param diag_name Requested diagnostic/field name
 * @param metadata Context for creating diagnostics
 * @return Unique pointer to diagnostic, or nullptr if name is not a diagnostic
 */
std::unique_ptr<DerivedDiagnostic>
create_diagnostic(const std::string &diag_name,
                  const DiagnosticMetadata &metadata);

/**
 * @brief Check if a field name is a derived diagnostic pattern.
 * @param name Field name to check
 * @return true if name matches a diagnostic pattern
 */
bool is_derived_diagnostic(const std::string &name);

/**
 * @brief Extract base field name from diagnostic name.
 *
 * Examples:
 * - "T_horiz_avg" → "T"
 * - "T_at_lev3" → "T"
 *
 * @param diag_name Diagnostic name
 * @return Base field name
 */
std::string get_base_field_name(const std::string &diag_name);

} // namespace emulator

#endif // EMULATOR_DIAGNOSTIC_FACTORY_HPP
