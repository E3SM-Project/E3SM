/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "FunctorsBuffersManager.hpp"
#include "ErrorDefs.hpp"

namespace Homme {

FunctorsBuffersManager::FunctorsBuffersManager ()
{
  m_size      = 0;
  m_allocated = false;
}

void FunctorsBuffersManager::request_size (const int num_doubles) {
  m_size = std::max(num_doubles, m_size);
}

void FunctorsBuffersManager:: allocate () {
  Errors::runtime_check(!m_allocated, "Error! Cannot call 'allocate' more than once.\n");

  m_buffer = ExecViewManaged<Real*>("",m_size);
  m_allocated = true;
}

} // namespace Homme
