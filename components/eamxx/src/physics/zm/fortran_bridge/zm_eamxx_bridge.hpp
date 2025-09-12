#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

#include "zm_functions.hpp"

// Bridge functions to call fortran version of ZM functions from C++

namespace scream {
namespace zm {

using ZMF = zm::Functions<Real, DefaultDevice>;

// Glue functions to call fortran from from C++ with the Data struct
void zm_eamxx_bridge_init( Int pcols, Int pver );
void zm_eamxx_bridge_run( Int ncol, Int pver,
                          ZMF::zm_input_state& zm_input,
                          ZMF::zm_output_tend& zm_output,
                          ZMF::zm_runtime_opt& zm_opts
                        );

extern "C" { // _f function decls
}

}  // namespace zm
}  // namespace scream

