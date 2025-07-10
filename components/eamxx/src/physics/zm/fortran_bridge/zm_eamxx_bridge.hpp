#include "share/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

// Bridge functions to call fortran version of ZM functions from C++

namespace scream {
namespace zm {

// Glue functions to call fortran from from C++ with the Data struct
void zm_eamxx_bridge_init();

extern "C" { // _f function decls
}

}  // namespace zm
}  // namespace scream

