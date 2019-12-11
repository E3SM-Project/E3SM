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

DirkFunctor::DirkFunctor (int nelem) {
  m_dirk_impl.reset(new DirkFunctorImpl(nelem));
}

// Note: you cannot declare the default destructor in the header,
//       since its implementation requires a definition of whatever
//       the unique ptr is pointing to, defying the pimpl idiom purpose.
//       To fix this empasse, simply declare the destructor, and define
//       it in the cpp file.
DirkFunctor::~DirkFunctor () = default;

void DirkFunctor::run (int n0, int nm1, int np1, Real alphadt, Real dt2,
                       const Elements& elements, const HybridVCoord& hvcoord) {
  m_dirk_impl->run(n0, nm1, np1, alphadt, dt2, elements, hvcoord);
}

} // Namespace Homme
