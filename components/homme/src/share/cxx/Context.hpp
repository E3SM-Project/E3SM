/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CONTEXT_HPP
#define HOMMEXX_CONTEXT_HPP

#include <string>
#include <map>
#include <memory>

namespace Homme {

class CaarFunctor;
class Derivative;
class Diagnostics;
class Elements;
class Tracers;
class HybridVCoord;
class HyperviscosityFunctor;
class SimulationParams;
class SphereOperators;
class TimeLevel;
class VerticalRemapManager;
class EulerStepFunctor;

/* A Context manages resources previously treated as singletons. Context is
 * meant to have two roles. First, a Context singleton is the only singleton in
 * the program. Second, a context need not be a singleton, and each Context
 * object can have different Elements, Derivative, etc., objects. (That
 * probably isn't needed, but Context immediately supports it.)
 *
 * Finally, Context has two singleton functions: singleton(), which returns
 * Context&, and finalize_singleton(). The second is called in a unit test exe
 * main before Kokkos::finalize().
 */
class Context {
private:
  // Note: using uniqe_ptr disables copy construction
  std::unique_ptr<CaarFunctor>            caar_functor_;
  std::unique_ptr<Elements>               elements_;
  std::unique_ptr<Tracers>                tracers_;
  std::unique_ptr<Derivative>             derivative_;
  std::unique_ptr<Diagnostics>            diagnostics_;
  std::unique_ptr<HybridVCoord>           hvcoord_;
  std::unique_ptr<HyperviscosityFunctor>  hyperviscosity_functor_;
  std::unique_ptr<SimulationParams>       simulation_params_;
  std::unique_ptr<TimeLevel>              time_level_;
  std::unique_ptr<VerticalRemapManager>   vertical_remap_mgr_;
  std::unique_ptr<SphereOperators>        sphere_operators_;
  std::unique_ptr<EulerStepFunctor>       euler_step_functor_;

  // Clear the objects Context manages.
  void clear();

public:
  Context();
  virtual ~Context();

  // Getters for each managed object.
  CaarFunctor& get_caar_functor();
  Diagnostics& get_diagnostics();
  Elements& get_elements();
  Tracers& get_tracers();
  Derivative& get_derivative();
  HybridVCoord& get_hvcoord();
  HyperviscosityFunctor& get_hyperviscosity_functor();
  SimulationParams& get_simulation_params();
  SphereOperators& get_sphere_operators();
  TimeLevel& get_time_level();
  EulerStepFunctor& get_euler_step_functor();
  VerticalRemapManager& get_vertical_remap_manager();

  // Exactly one singleton.
  static Context& singleton();

  static void finalize_singleton();
};

}

#endif // HOMMEXX_CONTEXT_HPP
