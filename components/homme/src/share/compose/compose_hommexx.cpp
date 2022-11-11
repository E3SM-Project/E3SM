#include "compose_hommexx.hpp"
#include "compose_homme.hpp"
#include "compose_slmm_islmpi.hpp"

namespace homme {

islmpi::IslMpi<>::Ptr get_isl_mpi_singleton();

bool cedr_should_run();
void cedr_sl_run_global();
void cedr_sl_run_local(const int limiter_option);
void cedr_sl_check();

void slmm_finalize();
void cedr_finalize();

namespace compose {

#ifdef COMPOSE_PORT
template <typename DataType>
using View = typename TracerArrays<ko::MachineTraits>::View<DataType>;
#endif

void set_views (const SetView<double***>& spheremp,
                const SetView<double****>& dp, const SetView<double*****>& dp3d,
                const SetView<double******>& qdp, const SetView<double*****>& q,
                const SetView<double*****>& dep_points) {
#ifdef COMPOSE_PORT
  auto& ta = *get_tracer_arrays();
  const auto nel = spheremp.extent_int(0);
  const auto np2 = spheremp.extent_int(1)*spheremp.extent_int(1);
  const auto nlev = dp.extent_int(3);
  ta.spheremp = View<Real**>(spheremp.data(), nel, np2);
  ta.dp = View<Real***>(dp.data(), nel, np2, nlev);
  ta.dp3d = View<Real****>(dp3d.data(), nel, dp3d.extent_int(1), np2, nlev);
  ta.qdp = View<Real*****>(qdp.data(), nel, qdp.extent_int(1), qdp.extent_int(2), np2, nlev);
  ta.q = View<Real****>(q.data(), nel, q.extent_int(1), np2, nlev);
  ta.dep_points = View<Real***[3]>(dep_points.data(), nel, dep_points.extent_int(1), np2);
#else
  slmm_throw_if(true, "Running a Hommexx code path with the non-Hommexx build"
                " is not supported.\n");
#endif
}

void advect (const int np1, const int n0_qdp, const int np1_qdp) {
  auto& cm = *get_isl_mpi_singleton();
  cm.tracer_arrays->np1 = np1;
  cm.tracer_arrays->n0_qdp = n0_qdp;
  cm.tracer_arrays->n1_qdp = np1_qdp;
  islmpi::step<>(cm, 0, cm.nelemd - 1, nullptr, nullptr, nullptr);
}

void set_dp3d_np1 (const int np1) {
  auto& cm = *get_isl_mpi_singleton();
  cm.tracer_arrays->np1 = np1;
}

bool property_preserve_global () {
  if ( ! cedr_should_run()) return false;
  homme::cedr_sl_run_global();
  return true;
}

bool property_preserve_local (const int limiter_option) {
  if ( ! cedr_should_run()) return false;
  homme::cedr_sl_run_local(limiter_option);
  return true;
}

void property_preserve_check () {
  homme::cedr_sl_check();
}

void finalize () {
  slmm_finalize();
  cedr_finalize();
}

} // namespace compose
} // namespace homme
