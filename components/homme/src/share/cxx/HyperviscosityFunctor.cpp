/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "HyperviscosityFunctor.hpp"
#include "HyperviscosityFunctorImpl.hpp"

#include "Context.hpp"
#include "Elements.hpp"
#include "SimulationParams.hpp"

namespace Homme
{

HyperviscosityFunctor::HyperviscosityFunctor ()
{
  auto& c = Context::singleton();
  auto& params   = c.get<SimulationParams>();
  auto& elements = c.get<Elements>();

  m_hvf_impl.reset (new HyperviscosityFunctorImpl(params,elements));
}

HyperviscosityFunctor::~HyperviscosityFunctor ()
{
  // This empty destructor (where HyperviscosityFunctorImpl type is completely known)
  // is necessary for pimpl idiom to work with unique_ptr. The issue is the
  // deleter, which needs to know the size of the stored type, and which
  // would be called from the implicitly declared default destructor, which
  // would be in the header file, where HyperviscosityFunctorImpl type is incomplete.
}

void HyperviscosityFunctor::init_boundary_exchanges () {
  assert (m_hvf_impl);
  m_hvf_impl->init_boundary_exchanges();
}

void HyperviscosityFunctor::run (const int np1, const Real dt, const Real eta_ave_w)
{
  // Sanity check (this should NEVER happen by design)
  assert (m_hvf_impl);

  m_hvf_impl->run(np1,dt,eta_ave_w);
}

} // namespace Homme
