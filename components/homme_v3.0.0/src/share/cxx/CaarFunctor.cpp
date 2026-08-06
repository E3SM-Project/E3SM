/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "CaarFunctor.hpp"
#include "CaarFunctorImpl.hpp"
#include "Context.hpp"
#include "Elements.hpp"
#include "ErrorDefs.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "ReferenceElement.hpp"
#include "SphereOperators.hpp"
#include "Tracers.hpp"
#include "mpi/MpiBuffersManager.hpp"

#include "profiling.hpp"

#include <assert.h>
#include <type_traits>

namespace Homme {

CaarFunctor::CaarFunctor()
  : is_setup(true)
{
  Elements&         elements   = Context::singleton().get<Elements>();
  Tracers&          tracers    = Context::singleton().get<Tracers>();
  ReferenceElement& ref_FE     = Context::singleton().get<ReferenceElement>();
  HybridVCoord&     hvcoord    = Context::singleton().get<HybridVCoord>();
  SphereOperators&  sphere_ops = Context::singleton().get<SphereOperators>();
  SimulationParams& params     = Context::singleton().get<SimulationParams>();

  // Build functor impl
  m_caar_impl.reset(new CaarFunctorImpl(elements,tracers,ref_FE,hvcoord,sphere_ops,params));
}

CaarFunctor::CaarFunctor(const Elements &elements, const Tracers &tracers,
                         const ReferenceElement &ref_FE,
                         const HybridVCoord &hvcoord,
                         const SphereOperators &sphere_ops,
                         const SimulationParams& params)
  : is_setup(true)
{
  // Build functor impl
  m_caar_impl.reset(new CaarFunctorImpl(elements, tracers, ref_FE, hvcoord,
                                        sphere_ops, params));
}

// This constructor is useful for using buffer functionality without
// having all other Functor information available.
// If this constructor is used, the setup() function must be called
// before using any other CaarFunctor functions.
CaarFunctor::CaarFunctor(const int num_elems, const SimulationParams& params)
  : is_setup(false)
{
  // Build functor impl
  m_caar_impl.reset(new CaarFunctorImpl(num_elems, params));
}

CaarFunctor::~CaarFunctor ()
{
  // This empty destructor (where CaarFunctorImpl type is completely known)
  // is necessary for pimpl idiom to work with unique_ptr. The issue is the
  // deleter, which needs to know the size of the stored type, and which
  // would be called from the implicitly declared default destructor, which
  // would be in the header file, where CaarFunctorImpl type is incomplete.
}

void CaarFunctor::setup(const Elements &elements, const Tracers &tracers,
                        const ReferenceElement &ref_FE, const HybridVCoord &hvcoord,
                        const SphereOperators &sphere_ops)
{
  assert (m_caar_impl);

  // Sanity check
  assert (!is_setup);

  m_caar_impl->setup(elements, tracers, ref_FE, hvcoord, sphere_ops);
  is_setup = true;
}

int CaarFunctor::requested_buffer_size () const {
  assert (m_caar_impl);
  return m_caar_impl->requested_buffer_size();
}

void CaarFunctor::init_buffers(const FunctorsBuffersManager& fbm) {
  assert (m_caar_impl);
  m_caar_impl->init_buffers(fbm);
}

void CaarFunctor::init_boundary_exchanges (const std::shared_ptr<MpiBuffersManager>& bm_exchange) {
  assert (m_caar_impl);

  // The Functor needs to be fully setup to use this function
  assert (is_setup);

  m_caar_impl->init_boundary_exchanges(bm_exchange);
}

void CaarFunctor::set_rk_stage_data (const RKStageData& data)
{
  // Sanity check (should NEVER happen)
  assert (m_caar_impl);

  // The Functor needs to be fully setup
  assert (is_setup);

  // Forward inputs to impl
  m_caar_impl->set_rk_stage_data(data);
}

void CaarFunctor::run (const RKStageData& data)
{
  // Sanity check (should NEVER happen)
  assert (m_caar_impl);

  // The Functor needs to be fully setup
  assert (is_setup);

  // Run functor
  m_caar_impl->run(data);
}

} // Namespace Homme
