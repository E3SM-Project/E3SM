#ifndef SCREAM_P3_INTERFACE_HPP
#define SCREAM_P3_INTERFACE_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_assert.hpp"

// Put everything into a scream namespace
namespace scream {

extern "C"
{

// Fortran routines to be called from C
void p3_init_f90 ();
void p3_standalone_init_f90 (Real* q, Real* T, Real* zi, Real* pmid, Real* dpres,
                             Real* ast, Real* ni_activated, Real* nc_nuceat_tend, Real* qv_prev, Real* T_prev );
void p3_main_f90 (const Real& dtime,
                  const Real* zi, const Real* pmid,
                  const Real* dpres, const Real* ast,
                  const Real* ni_activated, const Real* nc_nuceat_tend,
                  Real* q, Real* FQ, Real* T, Real* qv_prev, Real* T_prev);
void p3_finalize_f90 ();

} // extern "C"

} // namespace scream

#endif // SCREAM_P3_INTERFACE_HPP
