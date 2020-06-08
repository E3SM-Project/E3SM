#ifndef SCREAM_SHOC_FUNCTIONS_F90_HPP
#define SCREAM_SHOC_FUNCTIONS_F90_HPP

#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

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
struct SHOCGridData{
 // Inputs 
 Int shcol, nlev, nlevi;
 Real *zt_grid, *zi_grid, *pdel;

 // In/out
 Real *dz_zt, *dz_zi, *rho_zt;

};
void shoc_grid(Int nlev, SHOCGridData& d);

}  // namespace shoc
}  // namespace scream

#endif
