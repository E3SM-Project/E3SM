#ifndef SCREAM_SHOC_FUNCTIONS_F90_HPP
#define SCREAM_SHOC_FUNCTIONS_F90_HPP

#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

#include "shoc_functions.hpp"

#include <array>
#include <utility>

//
// Bridge functions to call fortran version of shoc functions from C++
//

namespace scream {
namespace shoc {

///////////////////////////////////////////////////////////////////////////////

// Converted subroutine helpers go here.

///////////////////////////////////////////////////////////////////////////////
// BFB math stuff
///////////////////////////////////////////////////////////////////////////////

extern "C" {

Real cxx_pow(Real base, Real exp);
Real cxx_sqrt(Real base);
Real cxx_cbrt(Real base);
Real cxx_log(Real input);
Real cxx_exp(Real input);

}

}  // namespace shoc
}  // namespace scream

#endif
