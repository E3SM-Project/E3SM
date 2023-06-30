#include "dp_functions_f90.hpp"

#include "dp_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include "share/util/scream_deep_copy.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C interface to DP fortran calls. The stubs below will link to fortran definitions in dp_iso_c.f90
//

extern "C" {

} // extern "C" : end _c decls

namespace scream {
namespace dp {

//
// Glue functions to call fortran from from C++ with the Data struct
//

// end _c impls

//
// _f function definitions. These expect data in C layout
//

} // namespace shoc
} // namespace scream
