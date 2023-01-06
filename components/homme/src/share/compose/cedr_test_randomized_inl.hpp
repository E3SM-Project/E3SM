// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TEST_RANDOMIZED_INL_HPP
#define INCLUDE_CEDR_TEST_RANDOMIZED_INL_HPP

#include "cedr_test_randomized.hpp"

namespace cedr {
namespace test {

template <typename CDRT, typename ES>
Int TestRandomized::run (const Int nrepeat, const bool write) {
  const Int nt = tracers_.size(), nlclcells = gcis_.size();

  Values v(nt, nlclcells);
  generate_rho(v);
  for (const auto& t : tracers_) {
    generate_Q(t, v);
    perturb_Q(t, v);
  }

  if (write)
    for (const auto& t : tracers_)
      write_pre(t, v);

  const typename CDRT::DeviceOp cdr
    = static_cast<const typename CDRT::DeviceOp&>(get_cdr().get_device_op());
  ValuesDevice<ES> vd(v);
  vd.sync_device();

  {
    const auto rhom = vd.rhom();
    const auto set_rhom = KOKKOS_LAMBDA (const Int& i) {
      cdr.set_rhom(i, 0, rhom[i]);
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, nlclcells), set_rhom);
  }
  // repeat > 1 runs the same values repeatedly for performance
  // meaurement.
  for (Int trial = 0; trial <= nrepeat; ++trial) {
    const auto set_Qm = KOKKOS_LAMBDA (const Int& j) {
      const auto ti = j / nlclcells;
      const auto i = j % nlclcells;
      cdr.set_Qm(i, ti, vd.Qm(ti)[i], vd.Qm_min(ti)[i], vd.Qm_max(ti)[i],
                 vd.Qm_prev(ti)[i]);
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, nt*nlclcells), set_Qm);
    run_impl(trial);
  }
  {
    const auto get_Qm = KOKKOS_LAMBDA (const Int& j) {
      const auto ti = j / nlclcells;
      const auto i = j % nlclcells;
      vd.Qm(ti)[i] = cdr.get_Qm(i, ti);
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, nt*nlclcells), get_Qm);
  }
  vd.sync_host(); // => v contains computed values

  if (write)
    for (const auto& t : tracers_)
      write_post(t, v);
  return check(cdr_name_, *p_, tracers_, v);
}

} // namespace test
} // namespace cedr

#endif
