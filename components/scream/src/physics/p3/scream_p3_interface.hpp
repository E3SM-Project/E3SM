#ifndef SCREAM_P3_INTERFACE_HPP
#define SCREAM_P3_INTERFACE_HPP

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

// Put everything into a scream namespace
namespace scream {

using scream::Real;
using scream::Int;
extern "C"
{

// Fortran routines to be called from C
void p3_init_f90     (Real* q);
void p3_main_f90     (const Real& dtime, Real* q);
void p3_finalize_f90 ();

} // extern "C"

} // namespace scream

#endif // SCREAM_P3_INTERFACE_HPP
