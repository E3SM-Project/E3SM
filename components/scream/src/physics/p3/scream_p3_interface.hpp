#ifndef SCREAM_P3_INTERFACE_HPP
#define SCREAM_P3_INTERFACE_HPP

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

// Put everything into a scream namespace
namespace scream {

extern "C"
{

// Fortran routines to be called from C
void p3_init_f90 ();
void p3_standalone_init_f90 (Real* q, Real* T, Real* zi, Real* pmid, Real* pdel,
                             Real* ast, Real* naai, Real* npccn );
void p3_main_f90 (const Real& dtime,
                  const Real* zi, const Real* pmid,
                  const Real* pdel, const Real* ast,
                  const Real* naai, const Real* npccn,
                  Real* q, Real* FQ, Real* T);
void p3_finalize_f90 ();

} // extern "C"

} // namespace scream

#endif // SCREAM_P3_INTERFACE_HPP
