/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_SESSION_HPP
#define HOMMEXX_SESSION_HPP

namespace Homme
{

extern "C" {

void initialize_hommexx_session ();
void finalize_hommexx_session ();

} // extern "C"

} // namespace Homme

#endif // HOMMEXX_SESSION_HPP
