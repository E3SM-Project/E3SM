#ifndef SCREAM_SHOC_MAIN_WRAP_HPP
#define SCREAM_SHOC_MAIN_WRAP_HPP

#include "share/scream_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace shoc {

struct FortranData;

// Run SHOC subroutines, populating inout and out fields of d.
ekat::Int shoc_main(FortranData& d, bool use_fortran);


// Test SHOC by running initial conditions for a number of steps and comparing
// against reference data. If gen_plot_scripts is true, Python scripts are
// emitted that plot initial and final conditions.
int test_shoc_ic(bool use_fortran, bool gen_plot_scripts = false);

}  // namespace shoc
}  // namespace scream

#endif
