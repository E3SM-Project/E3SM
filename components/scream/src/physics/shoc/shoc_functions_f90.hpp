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

struct SHOCGridData
{
    static constexpr size_t NUM_ARRAYS = 14;

    //input
    int shcol, nlev, nlevi;
    Real* zt_grid, zi_grid; 

    //output
};
///////////////////////////////////////////////////////////////////////////////

// Converted subroutine helpers go here.

}  // namespace shoc
}  // namespace scream

#endif
