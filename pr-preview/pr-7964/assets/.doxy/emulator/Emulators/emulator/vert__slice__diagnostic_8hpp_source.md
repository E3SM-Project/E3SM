

# File vert\_slice\_diagnostic.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**diagnostics**](dir_8c59288332f7499e80b506ac6e5a3903.md) **>** [**vert\_slice\_diagnostic.hpp**](vert__slice__diagnostic_8hpp.md)

[Go to the documentation of this file](vert__slice__diagnostic_8hpp.md)


```C++


#ifndef EMULATOR_VERT_SLICE_DIAGNOSTIC_HPP
#define EMULATOR_VERT_SLICE_DIAGNOSTIC_HPP

#include "derived_diagnostic.hpp"

namespace emulator {

class VertSliceDiagnostic : public DerivedDiagnostic {
public:
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
```


