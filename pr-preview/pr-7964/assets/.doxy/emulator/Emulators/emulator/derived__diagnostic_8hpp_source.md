

# File derived\_diagnostic.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**diagnostics**](dir_8c59288332f7499e80b506ac6e5a3903.md) **>** [**derived\_diagnostic.hpp**](derived__diagnostic_8hpp.md)

[Go to the documentation of this file](derived__diagnostic_8hpp.md)


```C++


#ifndef EMULATOR_DERIVED_DIAGNOSTIC_HPP
#define EMULATOR_DERIVED_DIAGNOSTIC_HPP

#include "../emulator_output_stream.hpp"
#include <memory>
#include <string>
#include <vector>

namespace emulator {

class DerivedDiagnostic {
public:
  virtual ~DerivedDiagnostic() = default;

  virtual std::string name() const = 0;

  virtual std::string source_field() const = 0;

  virtual void compute(const FieldDataProvider &fields,
                       std::vector<double> &output) = 0;

  virtual int output_size(int ncols, int nlevs) const = 0;
};

} // namespace emulator

#endif // EMULATOR_DERIVED_DIAGNOSTIC_HPP
```


