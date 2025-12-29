/**
 * @file horiz_avg_diagnostic.cpp
 * @brief Implementation of horizontal averaging diagnostic.
 */

#include "horiz_avg_diagnostic.hpp"
#include <algorithm>

namespace emulator {

HorizAvgDiagnostic::HorizAvgDiagnostic(const std::string &field_name,
                                       const std::vector<double> &area_weights,
                                       MPI_Comm comm)
    : m_source_field(field_name), m_area_weights(area_weights), m_comm(comm) {
  m_name = field_name + "_horiz_avg";
}

void HorizAvgDiagnostic::compute(const FieldDataProvider &fields,
                                 std::vector<double> &output) {
  const auto *field_data = fields.get_field(m_source_field);
  if (!field_data || field_data->empty()) {
    output.assign(1, 0.0);
    return;
  }

  int ncols = fields.get_ncols();
  int nlevs = fields.get_field_nlevs(m_source_field);
  output.resize(nlevs);

  // Check if we have area weights
  bool have_weights = (m_area_weights.size() == static_cast<size_t>(ncols));

  // Compute local weighted sum for each level
  std::vector<double> local_sum(nlevs, 0.0);
  double local_weight_sum = 0.0;

  for (int lev = 0; lev < nlevs; ++lev) {
    for (int col = 0; col < ncols; ++col) {
      int idx = lev * ncols + col; // Assuming [nlevs, ncols] layout
      if (idx >= static_cast<int>(field_data->size())) {
        // Fallback for [ncols] layout
        idx = col;
      }

      double weight = have_weights ? m_area_weights[col] : 1.0;
      local_sum[lev] += (*field_data)[idx] * weight;

      if (lev == 0) {
        local_weight_sum += weight;
      }
    }
  }

  // Global reduction
  std::vector<double> global_sum(nlevs);
  double global_weight_sum = 0.0;

  MPI_Allreduce(local_sum.data(), global_sum.data(), nlevs, MPI_DOUBLE, MPI_SUM,
                m_comm);
  MPI_Allreduce(&local_weight_sum, &global_weight_sum, 1, MPI_DOUBLE, MPI_SUM,
                m_comm);

  // Compute average
  if (global_weight_sum > 0.0) {
    for (int lev = 0; lev < nlevs; ++lev) {
      output[lev] = global_sum[lev] / global_weight_sum;
    }
  } else {
    std::fill(output.begin(), output.end(), 0.0);
  }
}

} // namespace emulator
