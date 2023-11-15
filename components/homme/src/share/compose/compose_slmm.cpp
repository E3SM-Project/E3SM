#include "compose.hpp"
#include "compose_homme.hpp"
#include "compose_slmm.hpp"
#include "compose_slmm_siqk.hpp"
#include "compose_slmm_advecter.hpp"
#include "compose_slmm_islmpi.hpp"

#include <sys/time.h>
#include <mpi.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <limits>
#include <algorithm>

namespace homme {
namespace islmpi {

#ifdef COMPOSE_PORT
template <typename ES> SLMM_KF
void print (const IslMpi<ko::MachineTraits>::ElemData<ES>& ed) {
  printf("me %p\n", ed.me);
  printf("nbrs %d:", ed.nbrs.size());
  for (int i = 0; i < ed.nbrs.size(); ++i)
    printf(" %d", ed.nbrs(i));
  printf("\nninhalo %d\n", ed.nin1halo);
  printf("own %d %d:", ed.own.capacity(), ed.own.size());
  for (int i = 0; i < ed.own.size(); ++i) {
    const auto& o = ed.own(i);
    printf(" %d %d,", (int) o.lev, (int) o.k);
  }
  printf("\nrmt %d %d:", ed.rmt.capacity(), ed.rmt.size());
  for (int i = 0; i < ed.rmt.size(); ++i) {
    const auto& r = ed.rmt(i);
    printf(" %d %d %d %d,", r.q_extrema_ptr, r.q_ptr, (int) r.lev, (int) r.k);
  }
  printf("\nsrc %d %d\n", ed.src.extent_int(0), ed.src.extent_int(1));
  printf("q_extrema %d %d\n", ed.q_extrema.extent_int(0), ed.q_extrema.extent_int(1));
}

template <typename MT, typename ESD, typename ESS>
void deep_copy (typename IslMpi<MT>::template ElemData<ESD>& d,
                const typename IslMpi<MT>::template ElemData<ESS>& s) {
  d.nbrs.copy(s.nbrs);
  const ptrdiff_t me_os = s.me - s.nbrs.data();
  d.me = d.nbrs.data() + me_os;
  d.nin1halo = s.nin1halo;
  //d.own.copy(s.own); unused in COMPOSE_PORT
  d.rmt.copy(s.rmt);
  siqk::resize_and_copy(d.src, s.src);
  d.q_extrema = typename IslMpi<MT>::template Array<Real**[2], ESD>(
    "q_extrema", s.q_extrema.extent_int(0), s.q_extrema.extent_int(1));
  ko::deep_copy(d.q_extrema, s.q_extrema);
  d.qdp = s.qdp;
  d.dp = s.dp;
  d.q = s.q;
}
#endif

template <typename MT>
void deep_copy (typename IslMpi<MT>::ElemDataListD& d,
                typename IslMpi<MT>::ElemDataListD::Mirror& m,
                const typename IslMpi<MT>::ElemDataListH& s) {
#ifdef COMPOSE_PORT
  const Int ned = s.size();
  // device view of device views
  d = typename IslMpi<MT>::ElemDataListD(ned);
  // host view of device views
  m = d.mirror();
  m.inc(ned);
  for (Int i = 0; i < ned; ++i)
    deep_copy<MT>(m(i), s(i));
  deep_copy(d, m);
# if 0
  for (int i = 0; i < ned; ++i) {
    printf(">> %d host\n", i);
    print(s(i));
    printf(">> %d device\n", i);
    ko::parallel_for(1, KOKKOS_LAMBDA (int) { print(d(i)); });
    ko::fence();
  }
# endif
#endif
}

template <typename MT> ko::EnableIfDiffSpace<MT>
sync_to_device (IslMpi<MT>& cm) { deep_copy<MT>(cm.ed_d, cm.ed_m, cm.ed_h); }

template <typename MT> ko::EnableIfSameSpace<MT>
sync_to_device (IslMpi<MT>& cm) { cm.ed_d = cm.ed_h; }

template <typename MT>
typename IslMpi<MT>::Ptr
init (const typename IslMpi<MT>::Advecter::ConstPtr& advecter,
      const mpi::Parallel::Ptr& p,
      Int np, Int nlev, Int qsize, Int qsized, Int nelemd,
      const Int* nbr_id_rank, const Int* nirptr,
      Int halo, Int traj_3d, Int traj_nsubstep) {
  slmm_throw_if(halo < 1, "halo must be 1 (default) or larger.");
  auto tracer_arrays = homme::init_tracer_arrays(nelemd, nlev, np, qsize, qsized);
  auto cm = std::make_shared<IslMpi<MT> >(p, advecter, tracer_arrays, np, nlev,
                                          qsize, qsized, nelemd, halo, traj_3d,
                                          traj_nsubstep);
  setup_comm_pattern(*cm, nbr_id_rank, nirptr);
  return cm;
}

// For const clarity, take the non-const advecter as an arg, even though cm
// already has a ref to the const'ed one.
template <typename MT>
void finalize_init_phase (IslMpi<MT>& cm, typename IslMpi<MT>::Advecter& advecter) {
  if (cm.halo > 1)
    extend_halo::extend_local_meshes<MT>(*cm.p, cm.ed_h, advecter);
  advecter.fill_nearest_points_if_needed();
  advecter.sync_to_device();
  sync_to_device(cm);
}

template <typename MT>
void set_hvcoord (IslMpi<MT>& cm, const Real etai_beg, const Real etai_end,
                  const Real* etam) {
  if (cm.etam.size() > 0) return;
#if defined COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    slmm_assert(cm.nlev > 0);
    cm.etai_beg = etai_beg;
    cm.etai_end = etai_end;
    cm.etam = typename IslMpi<MT>::template ArrayD<Real*>("etam", cm.nlev);
    const auto h = ko::create_mirror_view(cm.etam);
    for (int k = 0; k < cm.nlev; ++k) {
      h(k) = etam[k];
      slmm_assert(k == 0 || h(k) > h(k-1));
      slmm_assert(h(k) > 0 && h(k) < 1);
    }
    ko::deep_copy(cm.etam, h);
  }
#if defined COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

template void set_hvcoord(
  IslMpi<ko::MachineTraits>& cm, const Real etai_beg, const Real etai_end,
  const Real* etam);

// Set pointers to HOMME data arrays.
template <typename MT>
void set_elem_data (IslMpi<MT>& cm, const Int ie, Real* qdp, const Int n0_qdp,
                    const Real* dp, Real* q, const Int nelem_in_patch) {
  slmm_assert(ie < cm.ed_h.size());
  slmm_assert(cm.halo > 1 || cm.ed_h(ie).nbrs.size() == nelem_in_patch);
  auto& e = cm.ed_h(ie);
#if defined COMPOSE_PORT
  cm.tracer_arrays->pqdp.set_ie_ptr(ie, qdp);
  cm.tracer_arrays->n0_qdp = n0_qdp;
  cm.tracer_arrays->pdp.set_ie_ptr(ie, dp);
  cm.tracer_arrays->pq.set_ie_ptr(ie, q);
  e.qdp = e.dp = e.q = nullptr;
#else
  e.qdp = qdp + cm.nlev*cm.np2*cm.qsized*n0_qdp;
  e.dp = dp;
  e.q = q;
#endif
}

} // namespace islmpi

typedef ko::MachineTraits HommeMachineTraits;
typedef islmpi::IslMpi<HommeMachineTraits> HommeIslMpi;

static HommeIslMpi::Advecter::Ptr g_advecter;

void slmm_init (const Int np, const Int nelem, const Int nelemd,
                const Int transport_alg, const Int cubed_sphere_map,
                const slmm::Geometry::Type geometry, const Int sl_nearest_point_lev,
                const Int* lid2facenum) {
  g_advecter = std::make_shared<HommeIslMpi::Advecter>(
    np, nelemd, transport_alg, cubed_sphere_map, geometry, sl_nearest_point_lev);
  g_advecter->init_meta_data(nelem, lid2facenum);
}
} // namespace homme

namespace amb {
template <typename T> T strto(const char* s);
template <> inline int strto (const char* s) { return std::atoi(s); }
template <> inline bool strto (const char* s) { return std::atoi(s); }
template <> inline double strto (const char* s) { return std::atof(s); }
template <> inline std::string strto (const char* s) { return std::string(s); }

template <typename T>
bool getenv (const std::string& varname, T& var) {
  const char* var_s = std::getenv(varname.c_str());
  if ( ! var_s) return false;
  var = strto<T>(var_s);
  return true;
}

void dev_init_threads () {
#if defined COMPOSE_MIMIC_GPU
  static int nthr = -1;
  slmm_assert(omp_get_thread_num() == 0);
  if (nthr < 0) {
    nthr = 1;
    getenv("OMP_NUM_THREADS", nthr);
  }
  omp_set_num_threads(nthr);
  static_assert(std::is_same<ko::MachineTraits::DES, ko::OpenMP>::value,
                "in this dev code, should have OpenMP exe space on");
#endif
}

void dev_fin_threads () {
#if defined COMPOSE_MIMIC_GPU
  omp_set_num_threads(1);
#endif
}
} // namespace amb

namespace compose {
namespace test {
// Valid after slmm_init_local_mesh_ is called.
int slmm_unittest () {
  amb::dev_init_threads();
  int nerr = 0, ne;
  {
    ne = 0;
    const auto& p = homme::g_advecter->get_plane();
    const homme::Real length_scale
      = homme::g_advecter->is_sphere() ? 1 : std::max(p.Lx, p.Ly);
    for (int i = 0; i < homme::g_advecter->nelem(); ++i) {
      auto& m = homme::g_advecter->local_mesh_host(i);
      ne += slmm::unittest(m, m.tgt_elem, length_scale);
    }
    if (ne) fprintf(stderr, "FAIL slmm::unittest returned %d\n", ne);
    nerr += ne;
  }
  amb::dev_fin_threads();
  return nerr;
}
} // namespace test
} // namespace compose

#include <cstdlib>

namespace homme {
static homme::HommeIslMpi::Ptr g_csl_mpi;

HommeIslMpi::Ptr get_isl_mpi_singleton () { return g_csl_mpi; }

void slmm_finalize () {
#if defined COMPOSE_HORIZ_OPENMP
# pragma omp master
#endif
  {
    homme::g_csl_mpi = nullptr;
    homme::g_advecter = nullptr;
    homme::delete_tracer_arrays();
  }
}
} // namespace homme

static bool in_charge_of_kokkos = false;

static void initialize_kokkos () {
  if (Kokkos::is_initialized()) return;
  in_charge_of_kokkos = true;
  std::vector<char*> args;
#ifdef HOMMEXX_ENABLE_GPU
  int nd;
  const auto ret = cudaGetDeviceCount(&nd);
  if (ret != cudaSuccess) {
    // It isn't a big deal if we can't get the device count.
    nd = 1;
  }
  std::stringstream ss;
  ss << "--kokkos-ndevices=" << nd;
  const auto key = ss.str();
  std::vector<char> str(key.size()+1);
  std::copy(key.begin(), key.end(), str.begin());
  str.back() = 0;
  args.push_back(const_cast<char*>(str.data()));
#endif
  const char* silence = "--kokkos-disable-warnings";
  args.push_back(const_cast<char*>(silence));
  int narg = args.size();
  Kokkos::initialize(narg, args.data());
}

static void check_threading () {
#if defined COMPOSE_PORT && defined _OPENMP
  static bool first = true;
  if ( ! first) return;
  slmm_throw_if(omp_get_num_threads() > 1,
                "Compose API cannot be called in a threaded resion when "
                "COMPOSE_PORT is defined");
  first = false;
#endif
}

extern "C" {
// Interface for Homme, through compose_mod.F90.
void kokkos_init () {
  amb::dev_init_threads();
  initialize_kokkos();
  // Test these initialize correctly.
  Kokkos::View<int> v("hi");
  Kokkos::deep_copy(v, 0);
  homme::islmpi::FixedCapList<int,ko::MachineTraits::DES> fcl, fcl1(2);
  amb::dev_fin_threads();
}

void kokkos_finalize () {
  amb::dev_init_threads();
  if (in_charge_of_kokkos && Kokkos::is_initialized())
    Kokkos::finalize();
  amb::dev_fin_threads();
}

void slmm_init_impl (
  homme::Int fcomm, homme::Int transport_alg, homme::Int np,
  homme::Int nlev, homme::Int qsize, homme::Int qsized, homme::Int nelem,
  homme::Int nelemd, homme::Int cubed_sphere_map, homme::Int geometry,
  const homme::Int* lid2gid, const homme::Int* lid2facenum,
  const homme::Int* nbr_id_rank, const homme::Int* nirptr,
  homme::Int sl_halo, homme::Int sl_traj_3d, homme::Int sl_traj_nsubstep,
  homme::Int sl_nearest_point_lev,
  homme::Int, homme::Int, homme::Int, homme::Int)
{
  amb::dev_init_threads();
  homme::slmm_init(np, nelem, nelemd, transport_alg, cubed_sphere_map,
                   static_cast<slmm::Geometry::Type>(geometry),
                   sl_nearest_point_lev - 1, lid2facenum);
  slmm_throw_if(homme::g_advecter->is_cisl(), "CISL code was removed.");
  const auto p = homme::mpi::make_parallel(MPI_Comm_f2c(fcomm));
  homme::g_csl_mpi = homme::islmpi::init<homme::HommeMachineTraits>(
    homme::g_advecter, p, np, nlev, qsize, qsized, nelemd,
    nbr_id_rank, nirptr, sl_halo, sl_traj_3d, sl_traj_nsubstep);
  amb::dev_fin_threads();
}

void slmm_query_bufsz (homme::Int* sendsz, homme::Int* recvsz) {
  slmm_assert(homme::g_csl_mpi);
  if (ko::OnGpu<ko::MachineTraits::DES>::value)
    *sendsz = *recvsz = 0;
  else {
    homme::Int s = 0, r = 0;
    for (const auto e : homme::g_csl_mpi->sendsz) s += e;
    for (const auto e : homme::g_csl_mpi->recvsz) r += e;
    *sendsz = s;
    *recvsz = r;
  }
}

void slmm_init_plane (homme::Real Sx, homme::Real Sy, homme::Real Lx, homme::Real Ly) {
  slmm_assert(homme::g_advecter);
  homme::g_advecter->init_plane(Sx, Sy, Lx, Ly);
}

void slmm_set_bufs (homme::Real* sendbuf, homme::Real* recvbuf,
                    homme::Int, homme::Int) {
  slmm_assert(homme::g_csl_mpi);
#ifndef COMPOSE_WITH_HOMMEXX
  if (ko::OnGpu<ko::MachineTraits::DES>::value)
    sendbuf = recvbuf = nullptr;
#endif
  amb::dev_init_threads();
  homme::islmpi::alloc_mpi_buffers(*homme::g_csl_mpi, sendbuf, recvbuf);
  amb::dev_fin_threads();
}

void slmm_set_null_bufs () { slmm_set_bufs(nullptr, nullptr, 0, 0); }

void slmm_get_mpi_pattern (homme::Int* sl_mpi) {
  *sl_mpi = homme::g_csl_mpi ? 1 : 0;
}

void slmm_init_local_mesh (
  homme::Int ie, homme::Cartesian3D* neigh_corners, homme::Int nnc,
  homme::Cartesian3D* p_inside, homme::Int)
{
  amb::dev_init_threads();
  homme::g_advecter->init_local_mesh_if_needed(
    ie - 1, homme::FA3<const homme::Real>(
      reinterpret_cast<const homme::Real*>(neigh_corners), 3, 4, nnc),
    reinterpret_cast<const homme::Real*>(p_inside));
  amb::dev_fin_threads();
}

void slmm_init_finalize () {
  amb::dev_init_threads();
  if (homme::g_csl_mpi)
    homme::islmpi::finalize_init_phase(*homme::g_csl_mpi, *homme::g_advecter);
  amb::dev_fin_threads();
}

void slmm_check_ref2sphere (homme::Int ie, homme::Cartesian3D* p) {
  amb::dev_init_threads();
  homme::g_advecter->check_ref2sphere(
    ie - 1, reinterpret_cast<const homme::Real*>(p));
  amb::dev_fin_threads();
}

void slmm_set_hvcoord (const homme::Real etai_beg, const homme::Real etai_end,
                       const homme::Real* etam) {
  amb::dev_init_threads();
  slmm_assert(homme::g_csl_mpi);
  homme::islmpi::set_hvcoord(*homme::g_csl_mpi, etai_beg, etai_end, etam);
  amb::dev_fin_threads();
}

void slmm_calc_v_departure (
  homme::Int nets, homme::Int nete, homme::Int step, homme::Real dtsub,
  homme::Real* dep_points, homme::Int dep_points_ndim, homme::Real* vnode,
  homme::Real* vdep, homme::Int* info)
{
  amb::dev_init_threads();
  check_threading();
  slmm_assert(homme::g_csl_mpi);
  slmm_assert(homme::g_csl_mpi->sendsz.empty()); // alloc_mpi_buffers was called
  auto& cm = *homme::g_csl_mpi;
  slmm_assert(cm.dep_points_ndim == dep_points_ndim);
  {
    slmm::Timer timer("h2d");
    homme::sl_traj_h2d(*cm.tracer_arrays, dep_points, vnode, vdep,
                       cm.dep_points_ndim);
  }
  homme::islmpi::calc_v_departure(cm, nets - 1, nete - 1, step - 1,
                                  dtsub, dep_points, vnode, vdep);
  *info = 0;
  {
    slmm::Timer timer("d2h");
    homme::sl_traj_d2h(*cm.tracer_arrays, dep_points, vnode, vdep,
                       cm.dep_points_ndim);
  }
  amb::dev_fin_threads();
}

// Request extra data to be transferred for analysis.
static bool s_h2d, s_d2h;

void slmm_csl_set_elem_data (
  homme::Int ie, homme::Real* metdet, homme::Real* qdp, homme::Int n0_qdp,
  homme::Real* dp, homme::Real* q, homme::Int nelem_in_patch, bool h2d, bool d2h)
{
  amb::dev_init_threads();
  slmm_assert(homme::g_csl_mpi);
  homme::islmpi::set_elem_data(*homme::g_csl_mpi, ie - 1, qdp, n0_qdp - 1, dp, q,
                               nelem_in_patch);
  s_h2d = h2d;
  s_d2h = d2h;
  amb::dev_fin_threads();
}

void slmm_csl (homme::Int nets, homme::Int nete, homme::Real* dep_points,
               homme::Int dep_points_ndim, homme::Real* minq, homme::Real* maxq,
               homme::Int* info) {
  amb::dev_init_threads();
  check_threading();
  slmm_assert(homme::g_csl_mpi);
  slmm_assert(homme::g_csl_mpi->sendsz.empty()); // alloc_mpi_buffers was called
  auto& cm = *homme::g_csl_mpi;
  slmm_assert(cm.dep_points_ndim == dep_points_ndim);
  {
    slmm::Timer timer("h2d");
    homme::sl_h2d(*cm.tracer_arrays, s_h2d, dep_points, cm.dep_points_ndim);
  }
  *info = 0;
#if 1
  try {
    homme::islmpi::step(cm, nets - 1, nete - 1, dep_points, minq, maxq);
  } catch (const std::exception& e) {
    std::cerr << e.what();
    *info = -1;
  }
#else
  homme::islmpi::step(cm, nets - 1, nete - 1, dep_points, minq, maxq);
#endif
  {
    slmm::Timer timer("d2h");
    homme::sl_d2h(*cm.tracer_arrays, s_d2h, dep_points, cm.dep_points_ndim,
                  minq, maxq);
  }
  amb::dev_fin_threads();
}

void slmm_finalize () { homme::slmm_finalize(); }
} // extern "C"
