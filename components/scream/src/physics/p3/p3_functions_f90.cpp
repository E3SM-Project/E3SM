#include "p3_functions_f90.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"
#include "p3_f90.hpp"

using scream::Real;
using scream::Int;

extern "C" {

void p3_init_a_c(Real* itab, Real* itabcol);

}

namespace scream {
namespace p3 {

void p3_init_a(P3InitAFortranData& d)
{
  p3_init(); // need to initialize p3 first so that tables are loaded
  p3_init_a_c(d.itab.data(), d.itabcol.data());
}

} // namespace p3
} // namespace scream
