/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CAAR_FUNCTOR_HPP
#define HOMMEXX_CAAR_FUNCTOR_HPP

#include "Derivative.hpp"
#include "HybridVCoord.hpp"
#include "SphereOperators.hpp"
#include "Types.hpp"
#include <memory>

namespace Homme {

struct BuffersManager;
struct CaarFunctorImpl;

class Elements;
class Tracers;

class CaarFunctor {
public:
  CaarFunctor();
  CaarFunctor(const Elements &elements, const Tracers &tracers,
              const Derivative &derivative, const HybridVCoord &hvcoord,
              const SphereOperators &sphere_ops, 
              const int rsplit);
  CaarFunctor(const CaarFunctor &) = delete;

  ~CaarFunctor();

  CaarFunctor &operator=(const CaarFunctor &) = delete;

  void
  init_boundary_exchanges(const std::shared_ptr<BuffersManager> &bm_exchange);

  void set_n0_qdp(const int n0_qdp);

  void set_rk_stage_data(const int nm1, const int n0, const int np1,
                         const Real dt, const Real eta_ave_w,
                         const bool compute_diagnostics);

  void run();

  void run(const int nm1, const int n0, const int np1, const Real dt,
           const Real eta_ave_w, const bool compute_diagnostics);

private:
  std::unique_ptr<CaarFunctorImpl> m_caar_impl;

  // Setup the policies
  Kokkos::TeamPolicy<ExecSpace, void> m_policy;
};

} // Namespace Homme

#endif // HOMMEXX_CAAR_FUNCTOR_HPP
