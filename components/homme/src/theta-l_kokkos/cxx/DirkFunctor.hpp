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

class FunctorsBuffersManager;
class DirkFunctorImpl;
class Elements;
class HybridVCoord;

class DirkFunctor {
public:
  DirkFunctor(const int nelem);
  DirkFunctor(const DirkFunctor &) = delete;
  DirkFunctor &operator=(const DirkFunctor &) = delete;

  ~DirkFunctor();

  int requested_buffer_size() const;
  void init_buffers(const FunctorsBuffersManager& fbm);

  // Top-level interface, equivalent to compute_stage_value_dirk.
  void run(int nm1, Real alphadt_nm1, int n0, Real alphadt_n0, int np1, Real dt2,
           const Elements& elements, const HybridVCoord& hvcoord);

private:
  std::unique_ptr<DirkFunctorImpl> m_dirk_impl;
};

} // Namespace Homme

#endif // HOMMEXX_DIRK_FUNCTOR_HPP
