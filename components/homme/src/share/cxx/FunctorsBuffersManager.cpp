/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "FunctorsBuffersManager.hpp"
#include "ErrorDefs.hpp"
#ifdef HOMMEXX_BFB_TESTING
#include "utilities/TestUtils.hpp"
#include <random>
#endif

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

  generate_random_data();

  m_allocated = true;
}

void FunctorsBuffersManager:: allocate (Real* data, const int new_size) {
  Errors::runtime_check(!m_allocated, "Error! Cannot call 'allocate' more than once.\n");

  m_size   = new_size;
  m_buffer = ExecViewManaged<Real*>(data, m_size);

  generate_random_data();

  m_allocated = true;
}

void FunctorsBuffersManager::generate_random_data ()
{
#ifdef HOMMEXX_BFB_TESTING
  // Here's the catch: when malloc allocates memory, it *may* get a fresh
  // new page from the OS. If that happens, the OS will zero-out the provided
  // memory, for security reasons (e.g., you may get other people passwords).
  // In functors unit tests, this could accidentally make things work.
  // E.g., if the functor does not initialize its temps to 0, and sums into
  // them, a unit test would pass with zero-ed memory.
  // Therefore, we randomly initialize the buffer, to protect against false
  // positives.
  using rngalg = std::mt19937_64;
  using rpdf = std::uniform_real_distribution<Real>;
  std::random_device rd;
  rngalg engine(rd());
  genRandArray(m_buffer,engine,rpdf(0.1,10.0));
#endif
}

} // namespace Homme
