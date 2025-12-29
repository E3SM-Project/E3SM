/**
 * @file horiz_avg_diagnostic.hpp
 * @brief Horizontal averaging diagnostic.
 */

#ifndef EMULATOR_HORIZ_AVG_DIAGNOSTIC_HPP
#define EMULATOR_HORIZ_AVG_DIAGNOSTIC_HPP

#include "derived_diagnostic.hpp"
#include <mpi.h>

namespace emulator {

/**
 * @brief Computes area-weighted horizontal average of a field.
 *
 * Reduces a field from [ncols] or [ncols, nlevs] to a scalar (or [nlevs]).
 *
 * ## Pattern
 * Triggered by field names ending in "_horiz_avg" or "_global_mean".
 *
 * ## Example
 * - Input: "surface_temperature" with shape [ncols]
 * - Request: "surface_temperature_horiz_avg"
 * - Output: scalar (global mean)
 */
class HorizAvgDiagnostic : public DerivedDiagnostic {
public:
  /**
   * @brief Construct horizontal averaging diagnostic.
   * @param field_name Source field name
   * @param area_weights Area weights for each column (normalized to sum to 1)
   * @param comm MPI communicator for global reduction
   */
  HorizAvgDiagnostic(const std::string &field_name,
                     const std::vector<double> &area_weights, MPI_Comm comm);

  std::string name() const override { return m_name; }
  std::string source_field() const override { return m_source_field; }

  void compute(const FieldDataProvider &fields,
               std::vector<double> &output) override;

  int output_size(int ncols, int nlevs) const override {
    (void)ncols;
    return nlevs; // Returns one value per level (or 1 for 2D fields)
  }

private:
  std::string m_name;
  std::string m_source_field;
  std::vector<double> m_area_weights;
  MPI_Comm m_comm;
};

} // namespace emulator

#endif // EMULATOR_HORIZ_AVG_DIAGNOSTIC_HPP
