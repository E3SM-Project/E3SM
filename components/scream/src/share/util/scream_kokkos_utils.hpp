#ifndef SCREAM_KOKKOS_UTILS_HPP
#define SCREAM_KOKKOS_UTILS_HPP

#include "share/scream_kokkos_meta.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_std_meta.hpp"
#include "share/util/scream_utils.hpp"
#include "share/util/scream_arch.hpp"

#include <Kokkos_Core.hpp>

#include <cassert>
#include <type_traits>

// This file should not be merged with share/scream_kokkos_meta.hpp.
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
 * thread teams.
 */
template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class TeamUtils
{
  int _team_size, _num_teams, _max_threads;
#ifdef KOKKOS_ENABLE_CUDA
  using Device = Kokkos::Device<ExeSpace, typename ExeSpace::memory_space>;
  using flag_type = int;
  using view_1d = typename KokkosTypes<Device>::view_1d<flag_type>;
  using view_1d_slot = typename KokkosTypes<Device>::view_1d<int>;

  view_1d _open_ws_slots;
  view_1d_slot _ws_idxs;
#endif

 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy) : _team_size(0)
  {
    _max_threads = ExeSpace::concurrency() / ( (!is_single_precision<Real>::value && OnGpu<ExeSpace>::value) ? 2 : 1);
    const int team_size = policy.team_size();
    _num_teams = _max_threads / team_size;
    _team_size = _max_threads / _num_teams;

#ifdef KOKKOS_ENABLE_CUDA
    _open_ws_slots = view_1d("open_ws_slots", _num_teams);
    _ws_idxs       = view_1d_slot("ws_slots", policy.league_size());
#endif
  }

  // How many thread teams can run concurrently
  int get_num_concurrent_teams() const { return _num_teams; }

  // How many threads can run concurrently
  int get_max_concurrent_threads() const { return _max_threads; }

  /*
   * Of the C concurrently running teams, which "slot" is open
   * for the given team.
   */
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
  {
#ifdef KOKKOS_ENABLE_OPENMP
    return omp_get_thread_num() / _team_size;
#elif defined(KOKKOS_ENABLE_CUDA)
    // if (team_member.league_size() <= _num_teams) {
    if (false) {
      return team_member.league_rank();
    }
    else {
      Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
        //int my_ws_idx = team_member.league_rank() % _num_teams;
        int ws_idx = 0;
        while (!Kokkos::atomic_compare_exchange_strong(&_open_ws_slots(ws_idx), (flag_type) 0/*false*/, (flag_type)1/*true*/)) {
          // or random?
          ws_idx = (ws_idx+1) % _num_teams;
        }
        _ws_idxs(team_member.league_rank()) = ws_idx;
      });
      team_member.team_barrier();
      return _ws_idxs(team_member.league_rank());
    }
#else
    return 0;
#endif
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void release_workspace_idx(const MemberType& team_member, int ws_idx) const
  {
#ifdef KOKKOS_ENABLE_CUDA
    //if (team_member.league_size() > _num_teams) {
    if (true) {
      team_member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
        Kokkos::atomic_exchange(&_open_ws_slots(ws_idx), (flag_type) 0);
          // flag_type volatile* const e = &_open_ws_slots(ws_idx);
          // *e = (flag_type)0;
      });
      team_member.team_barrier();
    }
#endif
  }
};

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
