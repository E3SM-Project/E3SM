/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#ifndef HOMMEXX_COMPOSE_TRANSPORT_HPP
#define HOMMEXX_COMPOSE_TRANSPORT_HPP

#include "Types.hpp"
#include <memory>

namespace Homme {

class FunctorsBuffersManager;
class ComposeTransportImpl;
class SimulationParams;
class TimeLevel;

class ComposeTransport {
public:
  ComposeTransport();
  ComposeTransport(const int num_elems);
  ComposeTransport(const ComposeTransport &) = delete;
  ComposeTransport &operator=(const ComposeTransport &) = delete;

  ~ComposeTransport();

  bool setup_needed() { return ! is_setup; }
  void setup();

  void reset(const SimulationParams& params);

  int requested_buffer_size() const;
  void init_buffers(const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges();

  void run(const TimeLevel& tl, const Real dt);
  void remap_q(const TimeLevel& tl);

  std::vector<std::pair<std::string, int> > run_unit_tests();

  typedef Kokkos::View<Real*****, Kokkos::LayoutRight> TestDepView;
  TestDepView::HostMirror test_trajectory(Real t0, Real t1, bool independent_time_steps);

  void test_2d(const bool bfb, const int nstep, std::vector<Real>& eval);

private:
  std::unique_ptr<ComposeTransportImpl> m_compose_impl;
  bool is_setup;
};

} // Namespace Homme

#endif // HOMMEXX_COMPOSE_TRANSPORT_HPP
#endif // HOMME_ENABLE_COMPOSE
