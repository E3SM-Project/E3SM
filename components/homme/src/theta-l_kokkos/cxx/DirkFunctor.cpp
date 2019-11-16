/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "DirkFunctor.hpp"
#include "DirkFunctorImpl.hpp"
#include "Context.hpp"

#include "profiling.hpp"

#include <assert.h>
#include <type_traits>

namespace Homme {

DirkFunctor::DirkFunctor()
{
  m_dirk_impl.reset(new DirkFunctorImpl());
}

void DirkFunctor::run ()
{
  assert (m_dirk_impl);
  m_dirk_impl->run();
}

} // Namespace Homme
