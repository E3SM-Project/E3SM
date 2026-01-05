

# File diagnostic\_factory.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**diagnostics**](dir_8c59288332f7499e80b506ac6e5a3903.md) **>** [**diagnostic\_factory.hpp**](diagnostic__factory_8hpp.md)

[Go to the documentation of this file](diagnostic__factory_8hpp.md)


```C++


#ifndef EMULATOR_DIAGNOSTIC_FACTORY_HPP
#define EMULATOR_DIAGNOSTIC_FACTORY_HPP

#include "derived_diagnostic.hpp"
#include <memory>
#include <mpi.h>

namespace emulator {

struct DiagnosticMetadata {
  std::vector<double> area_weights; 
  MPI_Comm comm = MPI_COMM_WORLD;   
  int nlevs = 1;                    
};

std::unique_ptr<DerivedDiagnostic>
create_diagnostic(const std::string &diag_name,
                  const DiagnosticMetadata &metadata);

bool is_derived_diagnostic(const std::string &name);

std::string get_base_field_name(const std::string &diag_name);

} // namespace emulator

#endif // EMULATOR_DIAGNOSTIC_FACTORY_HPP
```


