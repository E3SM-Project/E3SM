/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"

#include "mpi/BoundaryExchange.hpp"
#include "mpi/BuffersManager.hpp"
#include "mpi/Comm.hpp"
#include "mpi/Connectivity.hpp"

#include "CaarFunctor.hpp"
#include "Derivative.hpp"
#include "Diagnostics.hpp"
#include "Elements.hpp"
#include "Tracers.hpp"
#include "HybridVCoord.hpp"
#include "HyperviscosityFunctor.hpp"
#include "SphereOperators.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "VerticalRemapManager.hpp"
#include "EulerStepFunctor.hpp"

namespace Homme {

Context::Context() {}

Context::~Context() {}

CaarFunctor& Context::get_caar_functor() {
  if ( ! caar_functor_) {
    caar_functor_.reset(new CaarFunctor());
  }
  return *caar_functor_;
}

Comm& Context::get_comm() {
  if ( ! comm_) {
    comm_.reset(new Comm());
  }
  return *comm_;
}

void Context::create_comm(const int f_comm) {
  // You should NOT create a C MPI_Comm from F90 twice during the same execution
  assert (!comm_);

  MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
  comm_.reset(new Comm(c_comm));
}

Diagnostics& Context::get_diagnostics() {
  if ( ! diagnostics_) diagnostics_.reset(new Diagnostics());
  return *diagnostics_;
}

Elements& Context::get_elements() {
  if ( ! elements_) elements_.reset(new Elements());
  return *elements_;
}

Tracers& Context::get_tracers() {
  if ( ! tracers_) {
    auto &elems = get_elements();
    const auto &params = get_simulation_params();
    assert (params.params_set);
    tracers_.reset(new Tracers(elems.num_elems(), params.qsize));
  }
  return *tracers_;
}

HybridVCoord& Context::get_hvcoord() {
  if ( ! hvcoord_) hvcoord_.reset(new HybridVCoord());
  return *hvcoord_;
}

HyperviscosityFunctor& Context::get_hyperviscosity_functor() {
  if ( ! hyperviscosity_functor_) {
    Elements& e = get_elements();
    Derivative& d = get_derivative();
    SimulationParams& p = get_simulation_params();
    hyperviscosity_functor_.reset(new HyperviscosityFunctor(p,e,d));
  }
  return *hyperviscosity_functor_;
}

Derivative& Context::get_derivative() {
  //if ( ! derivative_) derivative_ = std::make_shared<Derivative>();
  if ( ! derivative_) derivative_.reset(new Derivative());
  return *derivative_;
}

SimulationParams& Context::get_simulation_params() {
  if ( ! simulation_params_) simulation_params_.reset(new SimulationParams());
  return *simulation_params_;
}

TimeLevel& Context::get_time_level() {
  if ( ! time_level_) time_level_.reset(new TimeLevel());
  return *time_level_;
}

VerticalRemapManager& Context::get_vertical_remap_manager() {
  if ( ! vertical_remap_mgr_) vertical_remap_mgr_.reset(new VerticalRemapManager());
  return *vertical_remap_mgr_;
}

std::shared_ptr<BuffersManager> Context::get_buffers_manager(short int exchange_type) {
  if ( ! buffers_managers_) {
    buffers_managers_.reset(new BMMap());
  }

  if (!(*buffers_managers_)[exchange_type]) {
    (*buffers_managers_)[exchange_type] = std::make_shared<BuffersManager>(get_connectivity());
  }
  return (*buffers_managers_)[exchange_type];
}

std::shared_ptr<Connectivity> Context::get_connectivity() {
  if ( ! connectivity_) connectivity_.reset(new Connectivity());
  return connectivity_;
}

SphereOperators& Context::get_sphere_operators() {
  if ( ! sphere_operators_) {
    Elements&   elements   = get_elements();
    Derivative& derivative = get_derivative();
    sphere_operators_.reset(new SphereOperators(elements,derivative));
  }
  return *sphere_operators_;
}

EulerStepFunctor& Context::get_euler_step_functor() {
  if ( ! euler_step_functor_) euler_step_functor_.reset(new EulerStepFunctor());
  return *euler_step_functor_;
}

void Context::clear() {
  comm_ = nullptr;
  elements_ = nullptr;
  tracers_ = nullptr;
  derivative_ = nullptr;
  diagnostics_ = nullptr;
  hvcoord_ = nullptr;
  hyperviscosity_functor_ = nullptr;
  connectivity_ = nullptr;
  buffers_managers_ = nullptr;
  simulation_params_ = nullptr;
  sphere_operators_ = nullptr;
  time_level_ = nullptr;
  vertical_remap_mgr_ = nullptr;
  caar_functor_ = nullptr;
  euler_step_functor_ = nullptr;
}

Context& Context::singleton() {
  static Context c;
  return c;
}

void Context::finalize_singleton() {
  singleton().clear();
}

} // namespace Homme
