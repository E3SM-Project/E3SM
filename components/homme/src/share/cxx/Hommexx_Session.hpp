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
  // finalize_hommexx_session does nothing. This can be useful if
  // Homme is handled externally by some other library/executable,
  // in which case we should let such lib/exe handle finalization.
  // This is because finalize_hommexx_session will also call
  // Kokkos::finalize, which can be catastrophic if the host
  // application still has allocated views and/or still needs
  // to use Kokkos.
  static bool m_inited;

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
