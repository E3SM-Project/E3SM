#ifndef SCREAM_KOKKOS_UTILS_HPP
#define SCREAM_KOKKOS_UTILS_HPP

#include "ekat/scream_kokkos_meta.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_std_meta.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/util/scream_arch.hpp"

#include "Kokkos_Random.hpp"

#include <chrono>
#include <cassert>
#include <type_traits>

// This file should not be merged with ekat/scream_kokkos_meta.hpp.
// That file contains functionalities that *ideally* should be in kokkos
// itself, and we can foresee being in kokkos in the future. Hence, that
// file may be gone at some point (in fact, it may be gone by the time
// you read this comment). This file, instead, contains functionalities
// that are probably not generic enough to appear in kokkos any time soon
// (or ever), and are more app-specific.

namespace scream {
namespace util {

// Retrieve the underlying scalar type from a MD array type
template<typename MDArrayType>
struct ValueType {
private:
  using no_ref = typename std::remove_reference<MDArrayType>::type;
  using no_arr = typename std::remove_all_extents<no_ref>::type;
public:
  using type = typename remove_all_pointers<no_arr>::type;
};

template<typename DataType>
struct GetRanks {
private:
  using dimension = typename Kokkos::View<DataType>::traits::dimension;
public:
  enum : int { rank = dimension::rank };
  enum : int { rank_dynamic = dimension::rank_dynamic };
};

template<typename DataTypeOut, typename DataTypeIn, typename... Props>
typename std::enable_if<GetRanks<DataTypeOut>::rank_dynamic==0,
                        typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>>>::type
reshape (Kokkos::View<DataTypeIn,Props...> view_in) {
  typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>> view_out(view_in.data());
  assert (view_in.size()==view_out.size());
  return view_out;
}

template<typename DataTypeOut, typename DataTypeIn, typename... Props>
typename std::enable_if<GetRanks<DataTypeOut>::rank_dynamic==1,
                        typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>>>::type
reshape (Kokkos::View<DataTypeIn,Props...> view_in,
         const int dim0) {
  typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>> view_out(view_in.data(),dim0);
  assert (view_in.size()==view_out.size());
  return view_out;
}

template<typename DataTypeOut, typename DataTypeIn, typename... Props>
typename std::enable_if<GetRanks<DataTypeOut>::rank_dynamic==2,
                        typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>>>::type
reshape (Kokkos::View<DataTypeIn,Props...> view_in,
         const int dim0, const int dim1) {
  typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>> view_out(view_in.data(),dim0,dim1);
  assert (view_in.size()==view_out.size());
  return view_out;
}

/*
 * ExeSpaceUtils is essentially a TeamPolicy factory. TeamPolicy objects
 * are what kokkos uses to define a thread layout (num teams, threads/team)
 * for a parallel kernel. On non-GPU archictures, we will generally have
 * thread teams of 1.
 */
template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
struct ExeSpaceUtils {
  using TeamPolicy = Kokkos::TeamPolicy<ExeSpace>;

  static TeamPolicy get_default_team_policy (Int ni, Int /* nk */) {
#ifdef SCREAM_MIMIC_GPU
    const int max_threads = ExeSpace::concurrency();
    const int team_size = max_threads < 7 ? max_threads : 7;
    return TeamPolicy(ni, team_size);
#else
    return TeamPolicy(ni, 1);
#endif
  }

  static TeamPolicy get_team_policy_force_team_size (Int ni, Int team_size) {
    return TeamPolicy(ni, team_size);
  }
};

/*
 * Specialization of above for Cuda execution space. Many GPU architectures can
 * support a great number of threads, so we'll need to expose additional
 * parallelism by having many threads per team.  This is due to having more
 * threads than the main kernel loop has indices.
 */
#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ExeSpaceUtils<Kokkos::Cuda> {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;

  static TeamPolicy get_default_team_policy (Int ni, Int nk) {
    return TeamPolicy(ni, std::min(128, 32*((nk + 31)/32)));
  }

  static TeamPolicy get_team_policy_force_team_size (Int ni, Int team_size) {
    return TeamPolicy(ni, team_size);
  }
};
#endif

/*
 * TeamUtils contains utilities for getting concurrency info for
 * thread teams. Don't use _TeamUtilsCommonBase directly, use
 * TeamUtils.
 */

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class _TeamUtilsCommonBase
{
 protected:
  int _team_size, _num_teams, _max_threads, _league_size;

 public:
  template <typename TeamPolicy>
  _TeamUtilsCommonBase(const TeamPolicy& policy)
  {
    _max_threads = ExeSpace::concurrency() / ( (!is_single_precision<Real>::value && OnGpu<ExeSpace>::value) ? 2 : 1);
    const int team_size = policy.team_size();
    _num_teams = _max_threads / team_size;
    _team_size = _max_threads / _num_teams;
    _league_size = policy.league_size();

    // We will never run more teams than the policy needs
    _num_teams = _num_teams > _league_size ? _league_size : _num_teams;

    scream_assert_msg(_num_teams > 0, "Should always be able to run at least 1 team. max_thrds=" << _max_threads << " team_size=" << team_size << " league_size=" << _league_size);
  }

  // How many thread teams can run concurrently
  int get_num_concurrent_teams() const { return _num_teams; }

  // How many threads can run concurrently
  int get_max_concurrent_threads() const { return _max_threads; }

  // How many ws slots are there
  int get_num_ws_slots() const { return _num_teams; }

  /*
   * Of the C concurrently running teams, which "slot" is open
   * for the given team.
   */
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& /*team_member*/) const
  { return 0; }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void release_workspace_idx(const MemberType& /*team_member*/, int /*ws_idx*/) const
  { }
};

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class TeamUtils : public _TeamUtilsCommonBase<ExeSpace>
{
 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy, const Real& = 1.0) :
    _TeamUtilsCommonBase<ExeSpace>(policy)
  { }
};

/*
 * Specialization for OpenMP execution space
 */
#ifdef KOKKOS_ENABLE_OPENMP
template <>
class TeamUtils<Kokkos::OpenMP> : public _TeamUtilsCommonBase<Kokkos::OpenMP>
{
 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy, const Real& = 1.0) :
    _TeamUtilsCommonBase<Kokkos::OpenMP>(policy)
  { }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& /*team_member*/) const
  { return omp_get_thread_num() / _team_size; }
};
#endif

/*
 * Specialization for Cuda execution space.
 */
#ifdef KOKKOS_ENABLE_CUDA
template <>
class TeamUtils<Kokkos::Cuda> : public _TeamUtilsCommonBase<Kokkos::Cuda>
{
  using Device = Kokkos::Device<Kokkos::Cuda, typename Kokkos::Cuda::memory_space>;
  using flag_type = int; // this appears to be the smallest type that correctly handles atomic operations
  using view_1d = typename KokkosTypes<Device>::view_1d<flag_type>;
  using RandomGenerator = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>;
  using rnd_type = typename RandomGenerator::generator_type;

  int             _num_ws_slots;    // how many workspace slots (potentially more than the num of concurrent teams due to overprovision factor)
  bool            _need_ws_sharing; // true if there are more teams in the policy than ws slots
  view_1d         _open_ws_slots;    // indexed by ws-idx, true if in current use, else false
  RandomGenerator _rand_pool;

 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy, const Real& overprov_factor = 1.0) :
    _TeamUtilsCommonBase<Kokkos::Cuda>(policy),
    _num_ws_slots(_league_size > _num_teams
                  ? (overprov_factor * _num_teams > _league_size ? _league_size : overprov_factor * _num_teams)
                  : _num_teams),
    _need_ws_sharing(_league_size > _num_ws_slots),
    _open_ws_slots("open_ws_slots", _need_ws_sharing ? _num_ws_slots : 0),
    _rand_pool()
  {
    if (_need_ws_sharing) {
      _rand_pool = RandomGenerator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }
  }

  // How many ws slots are there
  int get_num_ws_slots() const { return _num_ws_slots; }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
  {
    if (!_need_ws_sharing) {
      return team_member.league_rank();
    }
    else {
      int ws_idx = 0;
      Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
        ws_idx = team_member.league_rank() % _num_ws_slots;
        if (!Kokkos::atomic_compare_exchange_strong(&_open_ws_slots(ws_idx), (flag_type) 0, (flag_type)1)) {
          rnd_type rand_gen = _rand_pool.get_state(team_member.league_rank());
          ws_idx = Kokkos::rand<rnd_type, int>::draw(rand_gen) % _num_ws_slots;
          while (!Kokkos::atomic_compare_exchange_strong(&_open_ws_slots(ws_idx), (flag_type) 0, (flag_type)1)) {
            ws_idx = Kokkos::rand<rnd_type, int>::draw(rand_gen) % _num_ws_slots;
          }
        }
      });

      // broadcast the idx to the team with a simple reduce
      int ws_idx_max_reduce;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, 1), [&] (int, int& ws_idx_max) {
        ws_idx_max = ws_idx;
      }, Kokkos::Max<int>(ws_idx_max_reduce));
      team_member.team_barrier();
      return ws_idx_max_reduce;
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void release_workspace_idx(const MemberType& team_member, int ws_idx) const
  {
    if (_need_ws_sharing) {
      team_member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
        flag_type volatile* const e = &_open_ws_slots(ws_idx);
        *e = (flag_type)0;
      });
    }
  }
};
#endif

// Get a 1d subview of the i-th dimension of a 2d view
template <typename T, typename ...Parms> KOKKOS_FORCEINLINE_FUNCTION
ko::Unmanaged<Kokkos::View<T*, Parms...> >
subview (const Kokkos::View<T**, Parms...>& v_in, const int i) {
  scream_kassert(v_in.data() != nullptr);
  scream_kassert(i < v_in.extent_int(0));
  scream_kassert(i >= 0);
  return ko::Unmanaged<Kokkos::View<T*, Parms...> >(
    &v_in.impl_map().reference(i, 0), v_in.extent(1));
}

} // namespace util
} // namespace scream

#endif // SCREAM_KOKKOS_UTILS_HPP
