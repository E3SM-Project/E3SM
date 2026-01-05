

# File horiz\_avg\_diagnostic.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**diagnostics**](dir_8c59288332f7499e80b506ac6e5a3903.md) **>** [**horiz\_avg\_diagnostic.hpp**](horiz__avg__diagnostic_8hpp.md)

[Go to the documentation of this file](horiz__avg__diagnostic_8hpp.md)


```C++


#ifndef EMULATOR_HORIZ_AVG_DIAGNOSTIC_HPP
#define EMULATOR_HORIZ_AVG_DIAGNOSTIC_HPP

#include "derived_diagnostic.hpp"
#include <mpi.h>

namespace emulator {

class HorizAvgDiagnostic : public DerivedDiagnostic {
public:
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
```


