

# File vert\_slice\_diagnostic.cpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**diagnostics**](dir_8c59288332f7499e80b506ac6e5a3903.md) **>** [**vert\_slice\_diagnostic.cpp**](vert__slice__diagnostic_8cpp.md)

[Go to the documentation of this file](vert__slice__diagnostic_8cpp.md)


```C++


#include "vert_slice_diagnostic.hpp"

namespace emulator {

VertSliceDiagnostic::VertSliceDiagnostic(const std::string &field_name,
                                         int level_idx, int nlevs)
    : m_source_field(field_name), m_level_idx(level_idx), m_nlevs(nlevs) {
  m_name = field_name + "_at_lev" + std::to_string(level_idx);
}

void VertSliceDiagnostic::compute(const FieldDataProvider &fields,
                                  std::vector<double> &output) {
  const auto *field_data = fields.get_field(m_source_field);
  if (!field_data || field_data->empty()) {
    output.clear();
    return;
  }

  int ncols = fields.get_ncols();
  output.resize(ncols);

  // Check if source is 3D (stacked slices) or 2D
  int total_size = static_cast<int>(field_data->size());
  int detected_nlevs = total_size / ncols;

  if (detected_nlevs <= 1 || m_level_idx >= detected_nlevs) {
    // 2D field or invalid level - just copy the field
    for (int col = 0; col < ncols && col < total_size; ++col) {
      output[col] = (*field_data)[col];
    }
    return;
  }

  // Extract slice: assuming [nlevs, ncols] layout
  for (int col = 0; col < ncols; ++col) {
    int idx = m_level_idx * ncols + col;
    output[col] = (*field_data)[idx];
  }
}

} // namespace emulator
```


