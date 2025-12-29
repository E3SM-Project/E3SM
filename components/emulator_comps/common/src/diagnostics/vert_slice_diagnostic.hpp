/**
 * @file vert_slice_diagnostic.hpp
 * @brief Vertical slicing diagnostic.
 */

#ifndef EMULATOR_VERT_SLICE_DIAGNOSTIC_HPP
#define EMULATOR_VERT_SLICE_DIAGNOSTIC_HPP

#include "derived_diagnostic.hpp"

namespace emulator {

/**
 * @brief Extracts a single vertical level from a 3D field.
 *
 * Reduces a field from [ncols, nlevs] to [ncols].
 *
 * ## Pattern
 * Triggered by field names with "_at_lev{N}" suffix.
 *
 * ## Example
 * - Input: "air_temperature" with shape [ncols, 8]
 * - Request: "air_temperature_at_lev3"
 * - Output: "air_temperature" at level index 3, shape [ncols]
 */
class VertSliceDiagnostic : public DerivedDiagnostic {
public:
  /**
   * @brief Construct vertical slicing diagnostic.
   * @param field_name Source field name
   * @param level_idx Level index to extract (0-based)
   * @param nlevs Total number of levels in source field
   */
  VertSliceDiagnostic(const std::string &field_name, int level_idx, int nlevs);

  std::string name() const override { return m_name; }
  std::string source_field() const override { return m_source_field; }

  void compute(const FieldDataProvider &fields,
               std::vector<double> &output) override;

  int output_size(int ncols, int nlevs) const override {
    (void)nlevs;
    return ncols; // Returns one value per column
  }

private:
  std::string m_name;
  std::string m_source_field;
  int m_level_idx;
  int m_nlevs;
};

} // namespace emulator

#endif // EMULATOR_VERT_SLICE_DIAGNOSTIC_HPP
