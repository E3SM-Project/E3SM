

# File atm\_io.hpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_io.hpp**](atm__io_8hpp.md)

[Go to the documentation of this file](atm__io_8hpp.md)


```C++


#ifndef ATM_IO_HPP
#define ATM_IO_HPP

#include "../../../common/src/emulator_logger.hpp"
#include "atm_field_manager.hpp"
#include <string>
#include <vector>

namespace emulator {
namespace impl {

bool read_atm_initial_conditions(const std::string &filename,
                                 int num_global_cols, int num_local_cols,
                                 const std::vector<int> &col_gids,
                                 const std::vector<double> &lat,
                                 AtmFieldManager &fields,
                                 const std::vector<std::string> &required_vars,
                                 Logger &logger, bool is_root);

} // namespace impl
} // namespace emulator

#endif // ATM_IO_HPP
```


