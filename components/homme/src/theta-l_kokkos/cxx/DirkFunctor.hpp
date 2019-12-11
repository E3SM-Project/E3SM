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

class DirkFunctorImpl;
class Elements;
class HybridVCoord;

class DirkFunctor {
public:
  DirkFunctor(const int nelem);
  DirkFunctor(const DirkFunctor &) = delete;
  DirkFunctor &operator=(const DirkFunctor &) = delete;

  ~DirkFunctor();

  // Top-level interface, equivalent to compute_stage_value_dirk_stripped.
  //   Set alphadt = 0 if the alphadt-dependent term is not to be added to the RHS.
  //   Set nm1 = -1 if the nm1-dependent term is not to be added to the RHS.
  void run(int n0, int nm1, int np1, Real alphadt, Real dt2,
           const Elements& elements, const HybridVCoord& hvcoord);

private:
  std::unique_ptr<DirkFunctorImpl> m_dirk_impl;
};

} // Namespace Homme

#endif // HOMMEXX_DIRK_FUNCTOR_HPP
