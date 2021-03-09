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

int DirkFunctor::requested_buffer_size () const {
  return m_dirk_impl->requested_buffer_size();
}

void DirkFunctor::init_buffers (const FunctorsBuffersManager& fbm) {
  m_dirk_impl->init_buffers(fbm);
}

void DirkFunctor::run (int nm1, Real alphadt_nm1, int n0, Real alphadt_n0, int np1, Real dt2,
                       const Elements& elements, const HybridVCoord& hvcoord) {
  GPTLstart("compute_stage_value_dirk");
  m_dirk_impl->run(nm1, alphadt_nm1, n0, alphadt_n0, np1, dt2, elements, hvcoord);
  GPTLstop("compute_stage_value_dirk");
}

} // Namespace Homme
