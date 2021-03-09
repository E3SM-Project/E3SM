#ifndef KOKKOS_UTILS_HPP
#define KOKKOS_UTILS_HPP

#include "Types.hpp"
#include "ExecSpaceDefs.hpp"

#include "Kokkos_Random.hpp"

#include <chrono>
#include <cassert>
#include <type_traits>

// That file contains functionalities that *ideally* should be in kokkos
// itself, and we can foresee being in kokkos in the future. Hence, that
// file may be gone at some point (in fact, it may be gone by the time
// you read this comment). This file, instead, contains functionalities
// that are probably not generic enough to appear in kokkos any time soon
// (or ever), and are more app-specific.

namespace Homme {

/*
 * TeamUtils contains utilities for getting concurrency info for
 * thread teams. Don't use _TeamUtilsCommonBase directly, use
 * TeamUtils.
 */

template <typename ExeSpace=ExecSpace>
class _TeamUtilsCommonBase
{
 protected:
  int _team_size, _num_teams, _max_threads, _league_size;

 public:
  template <typename TeamPolicy>
  _TeamUtilsCommonBase(const TeamPolicy& policy)
  {
    _max_threads = ExeSpace::concurrency() / ( OnGpu<ExeSpace>::value ? 2 : 1);
    const int team_size = policy.team_size();
    _num_teams = _max_threads / team_size;
    _team_size = _max_threads / _num_teams;
    _league_size = policy.league_size();

    // We will never run more teams than the policy needs
    _num_teams = _num_teams > _league_size ? _league_size : _num_teams;

    assert(_league_size==0 || _num_teams > 0);
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
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
  using Device = Kokkos::Device<Kokkos::Cuda, typename Kokkos::Cuda::memory_space>;
  using flag_type = int; // this appears to be the smallest type that correctly handles atomic operations
  using view_1d = ExecView<flag_type*>;
  using RandomGenerator = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>;
  using rnd_type = typename RandomGenerator::generator_type;

  int             _num_ws_slots;    // how many workspace slots (potentially more than the num of concurrent teams due to overprovision factor)
  bool            _need_ws_sharing; // true if there are more teams in the policy than ws slots
  view_1d         _open_ws_slots;    // indexed by ws-idx, true if in current use, else false
  RandomGenerator _rand_pool;
#endif

 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy, const Real& overprov_factor = 1.25) :
    _TeamUtilsCommonBase<Kokkos::Cuda>(policy)
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
    , _num_ws_slots(_league_size > _num_teams
                    ? (overprov_factor * _num_teams > _league_size ? _league_size : overprov_factor * _num_teams)
                    : _num_teams)
    , _need_ws_sharing(_league_size > _num_ws_slots)
    , _open_ws_slots("open_ws_slots", _need_ws_sharing ? _num_ws_slots : 0)
    , _rand_pool()
#endif
  {
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
    if (_need_ws_sharing) {
      _rand_pool = RandomGenerator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }
#endif
  }

  // How many ws slots are there
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
  int get_num_ws_slots() const { return _num_ws_slots; }
#else
  int get_num_ws_slots() const { return _league_size; }
#endif

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
  {
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
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
#else
    return team_member.league_rank();
#endif
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void release_workspace_idx(const MemberType& team_member, int ws_idx) const
  {
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
    if (_need_ws_sharing) {
      team_member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
        flag_type volatile* const e = &_open_ws_slots(ws_idx);
        *e = (flag_type)0;
      });
    }
#endif
  }
};
#endif

} // namespace Homme

#endif // KOKKOS_UTILS_HPP
