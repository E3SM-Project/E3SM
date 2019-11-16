/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_DIRK_FUNCTOR_HPP
#define HOMMEXX_DIRK_FUNCTOR_HPP

#include "Types.hpp"
#include <memory>

namespace Homme {

struct DirkFunctorImpl;

class DirkFunctor {
public:
  DirkFunctor();
  DirkFunctor(const DirkFunctor &) = delete;
  DirkFunctor &operator=(const DirkFunctor &) = delete;

  void run();

private:
  std::unique_ptr<DirkFunctorImpl> m_dirk_impl;
};

} // Namespace Homme

#endif // HOMMEXX_DIRK_FUNCTOR_HPP
