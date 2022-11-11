/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CAAR_FUNCTOR_HPP
#define HOMMEXX_CAAR_FUNCTOR_HPP

#include "Types.hpp"
#include "RKStageData.hpp"
#include <memory>

namespace Homme {

struct CaarFunctorImpl;

class Elements;
class Tracers;
class ReferenceElement;
class SphereOperators;
class SimulationParams;
class HybridVCoord;

class MpiBuffersManager;
struct FunctorsBuffersManager;

class CaarFunctor {
public:
  CaarFunctor();
  CaarFunctor(const Elements &elements, const Tracers &tracers,
              const ReferenceElement &ref_FE, const HybridVCoord &hvcoord,
              const SphereOperators &sphere_ops, const SimulationParams& params);
  CaarFunctor(const int num_elems, const SimulationParams& params);
  CaarFunctor(const CaarFunctor &) = delete;

  ~CaarFunctor();

  CaarFunctor &operator=(const CaarFunctor &) = delete;

  bool setup_needed() { return !is_setup; }
  void setup(const Elements &elements, const Tracers &tracers,
             const ReferenceElement &ref_FE, const HybridVCoord &hvcoord,
             const SphereOperators &sphere_ops);

  int requested_buffer_size () const;
  void init_buffers(const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges(const std::shared_ptr<MpiBuffersManager> &bm_exchange);

  void set_rk_stage_data(const RKStageData& data);

  void run(const RKStageData& data);

private:
  std::unique_ptr<CaarFunctorImpl> m_caar_impl;
  bool is_setup;
};

} // Namespace Homme

#endif // HOMMEXX_CAAR_FUNCTOR_HPP
