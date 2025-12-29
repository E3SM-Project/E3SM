/**
 * @file derived_diagnostic.hpp
 * @brief Base class for derived diagnostics.
 *
 * Part of the modular diagnostics infrastructure in common/src/diagnostics/.
 */

#ifndef EMULATOR_DERIVED_DIAGNOSTIC_HPP
#define EMULATOR_DERIVED_DIAGNOSTIC_HPP

#include "../emulator_output_stream.hpp"
#include <memory>
#include <string>
#include <vector>

namespace emulator {

/**
 * @brief Base class for derived diagnostics.
 *
 * Derived diagnostics compute output values from one or more input fields.
 * Examples: horizontal averages, vertical slices, pressure interpolation.
 *
 * ## Implementing Custom Diagnostics
 *
 * To add a new diagnostic type:
 *
 * 1. Create a subclass of DerivedDiagnostic
 * 2. Implement all virtual methods
 * 3. Add pattern matching in diagnostic_factory.cpp
 */
class DerivedDiagnostic {
public:
  virtual ~DerivedDiagnostic() = default;

  /**
   * @brief Get the diagnostic name.
   */
  virtual std::string name() const = 0;

  /**
   * @brief Get the source field name.
   */
  virtual std::string source_field() const = 0;

  /**
   * @brief Compute the diagnostic.
   * @param fields Input field provider
   * @param output Output data buffer (caller must allocate)
   */
  virtual void compute(const FieldDataProvider &fields,
                       std::vector<double> &output) = 0;

  /**
   * @brief Get output size.
   * @param ncols Number of local columns
   * @param nlevs Number of vertical levels (for source field)
   */
  virtual int output_size(int ncols, int nlevs) const = 0;
};

} // namespace emulator

#endif // EMULATOR_DERIVED_DIAGNOSTIC_HPP
