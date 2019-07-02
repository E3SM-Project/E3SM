#ifndef SCREAM_P3_INTERFACE_HPP
#define SCREAM_P3_INTERFACE_HPP

// Put everything into a scream namespace
namespace scream {

extern "C"
{

// Fortran routines to be called from C
void p3_init_f90     ();
void p3_main_f90     ();
void p3_finalize_f90 ();

} // extern "C"

} // namespace scream

#endif // SCREAM_P3_INTERFACE_HPP
