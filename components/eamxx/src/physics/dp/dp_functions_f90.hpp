#ifndef SCREAM_DP_FUNCTIONS_F90_HPP
#define SCREAM_DP_FUNCTIONS_F90_HPP

#include "share/scream_types.hpp"
#include "physics/share/physics_test_data.hpp"

#include "dp_functions.hpp"
#include "physics_constants.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of dp functions from C++
//

namespace scream {
namespace dp {

// Glue functions to call fortran from from C++ with the Data struct

extern "C" { // _f function decls

} // end _f function decls

}  // namespace dp
}  // namespace scream

#endif // SCREAM_DP_FUNCTIONS_F90_HPP
