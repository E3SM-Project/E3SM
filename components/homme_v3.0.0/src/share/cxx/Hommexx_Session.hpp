/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_SESSION_HPP
#define HOMMEXX_SESSION_HPP

namespace Homme
{

struct Session {
  // Whether initialize_hommexx_session was called. If this is false,
  // initialize_hommexx_session does nothing.
  // This can be useful if Homme is handled externally by some other
  // library/executable, in which case we should let such lib/exe handle
  // initialization, so that Kokkos is initialized only once.
  static bool m_inited;

  // Whether this Session owns kokkos. If true, then initialize/finalize
  // hommexx session will init/finalize kokkos as well, otherwise it won't.
  // This can be useful if Homme is handled externally by some other
  // library/executable, in which case we should let such lib/exe handle
  // init/finalization of kokkos.
  static bool m_handle_kokkos;

  // If this is true, error routines in the Error namespace will
  // not call finalize_hommexx_session (and then MPI_Abort).
  // Instead, they will throw an exception. This allows a host
  // lib/exe that is calling homme to decide what to do.
  // In particular, this prevents hommexx from calling
  // kokkos finalization, which may be catastrophic for the host app.
  static bool m_throw_instead_of_abort;
};

extern "C" {

void initialize_hommexx_session ();
void finalize_hommexx_session ();

} // extern "C"

} // namespace Homme

#endif // HOMMEXX_SESSION_HPP
