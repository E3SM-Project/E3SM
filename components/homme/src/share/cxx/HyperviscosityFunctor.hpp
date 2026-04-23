/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HYPERVISCOSITY_FUNCTOR_HPP
#define HOMMEXX_HYPERVISCOSITY_FUNCTOR_HPP

#include "Types.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ElementsDerivedState.hpp"
#include "SimulationParams.hpp"

#include <memory>

namespace Homme
{

class HyperviscosityFunctorImpl;
struct FunctorsBuffersManager;
class SimulationParams;

class HyperviscosityFunctor
{
public:

  HyperviscosityFunctor ();

  HyperviscosityFunctor (const int num_elems, const SimulationParams& params);

  ~HyperviscosityFunctor ();

  bool setup_needed() { return !is_setup; }
  void setup(const ElementsGeometry&     geometry,
             const ElementsState&        state,
             const ElementsDerivedState& derived);

  int requested_buffer_size () const;
  void init_buffers    (const FunctorsBuffersManager& fbm);

  void init_boundary_exchanges();

  // Forward nu_scale_top from Fortran to the impl
  void init_nu_scale_top_from_f90 (HostViewUnmanaged<const Real[NUM_PHYSICAL_LEV]> nu_scale_top_h,
                                   const int nlev_tom);

  void run (const int np1, const Real dt, const Real eta_ave_w);

private:

  std::unique_ptr<HyperviscosityFunctorImpl>  m_hvf_impl;
  bool is_setup;
};

} // namespace Homme

#endif // HOMMEXX_HYPERVISCOSITY_FUNCTOR_HPP
