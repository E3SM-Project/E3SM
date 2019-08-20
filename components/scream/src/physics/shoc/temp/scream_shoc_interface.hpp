#ifndef SCREAM_SHOC_INTERFACE_HPP
#define SCREAM_SHOC_INTERFACE_HPP

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

// Put everything into a scream namespace
namespace scream {

using scream::Real;
using scream::Int;
extern "C"
{

// Fortran routines to be called from C
void shoc_c_init     ();
void shoc_c_main     (const Real& dtime, Real* q);//, Real* FQ, const Real* qdp);
void shoc_c_finalize ();

} // extern "C"

} // namespace scream

#endif // SCREAM_SHOC_INTERFACE_HPP
