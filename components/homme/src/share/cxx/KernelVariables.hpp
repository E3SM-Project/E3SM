/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef KERNEL_VARIABLES_HPP
#define KERNEL_VARIABLES_HPP

#include "Types.hpp"
#include "kokkos_utils.hpp"

namespace Homme {

struct KernelVariables {
private:
  struct TeamInfo {

    template<typename ExecSpaceType>
    static
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!std::is_same<ExecSpaceType,Hommexx_Cuda>::value &&
                            !std::is_same<ExecSpaceType,Hommexx_OpenMP>::value,
                            int
                           >::type
    get_team_idx (const int /*team_size*/, const int /*league_rank*/)
    {
      return 0;
    }

#ifdef KOKKOS_ENABLE_CUDA
#ifdef __CUDA_ARCH__
    template <typename ExecSpaceType>
    static KOKKOS_INLINE_FUNCTION typename std::enable_if<
        OnGpu<ExecSpaceType>::value, int>::type
    get_team_idx(const int /*team_size*/, const int league_rank) {
      return league_rank;
    }
#else
    template <typename ExecSpaceType>
    static KOKKOS_INLINE_FUNCTION typename std::enable_if<
        OnGpu<ExecSpaceType>::value, int>::type
    get_team_idx(const int /*team_size*/, const int /*league_rank*/) {
      assert(false); // should never happen
      return -1;
    }
#endif // __CUDA_ARCH__
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_OPENMP
    template<typename ExecSpaceType>
    static
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if<std::is_same<ExecSpaceType,Kokkos::OpenMP>::value,int>::type
    get_team_idx (const int team_size, const int /*league_rank*/)
    {
      // Kokkos puts consecutive threads into the same team.
      return omp_get_thread_num() / team_size;
    }
#endif
  };

public:

  KOKKOS_INLINE_FUNCTION
  KernelVariables(const TeamMember &team_in)
      : team(team_in)
      , ie(team_in.league_rank())
      , iq(-1)
      , team_idx(TeamInfo::get_team_idx<ExecSpace>(team_in.team_size(),team_in.league_rank()))
      , team_utils(nullptr)
  {
    // Nothing to be done here
  }

  KOKKOS_INLINE_FUNCTION
  KernelVariables(const TeamMember &team_in, const TeamUtils<ExecSpace>& utils)
      : team(team_in)
      , ie(team_in.league_rank())
      , iq(-1)
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
      , team_idx(utils.get_workspace_idx(team_in))
      , team_utils(&utils)
#else
      , team_idx(TeamInfo::get_team_idx<ExecSpace>(team_in.team_size(),team_in.league_rank()))
      , team_utils(nullptr)
#endif
  {
    // Nothing to be done here
  }

  KOKKOS_INLINE_FUNCTION
  KernelVariables(const TeamMember &team_in, const int qsize)
      : team(team_in)
      , ie(team_in.league_rank() / qsize)
      , iq(team_in.league_rank() % qsize)
      , team_idx(TeamInfo::get_team_idx<ExecSpace>(team_in.team_size(),team_in.league_rank()))
      , team_utils(nullptr)
  {
    // Nothing to be done here
  }

  KOKKOS_INLINE_FUNCTION
  KernelVariables(const TeamMember &team_in, const int qsize, const TeamUtils<ExecSpace>& utils)
      : team(team_in)
      , ie(team_in.league_rank() / qsize)
      , iq(team_in.league_rank() % qsize)
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
      , team_idx(utils.get_workspace_idx(team_in))
      , team_utils(&utils)
#else
      , team_idx(TeamInfo::get_team_idx<ExecSpace>(team_in.team_size(),team_in.league_rank()))
      , team_utils(nullptr)
#endif
  {
    // Nothing to be done here
  }

#ifdef HOMMEXX_CUDA_SHARE_BUFFER
  KOKKOS_INLINE_FUNCTION
  ~KernelVariables() {
    if (team_utils != nullptr) {
      team_utils->release_workspace_idx(team, team_idx);
    }
  }
#endif

  template <typename Primitive, typename Data>
  KOKKOS_INLINE_FUNCTION Primitive *allocate_team() const {
    ScratchView<Data> view(team.team_scratch(0));
    return view.data();
  }

  template <typename Primitive, typename Data>
  KOKKOS_INLINE_FUNCTION Primitive *allocate_thread() const {
    ScratchView<Data> view(team.thread_scratch(0));
    return view.data();
  }

  KOKKOS_INLINE_FUNCTION
  static size_t shmem_size(int team_size) {
    size_t mem_size = 0 * team_size;
    return mem_size;
  }

  const TeamMember &team;

  KOKKOS_FORCEINLINE_FUNCTION void team_barrier() const {
    team.team_barrier();
  }

  int ie, iq;
  const int team_idx;
  const TeamUtils<ExecSpace>* team_utils;
}; // KernelVariables

} // Homme

#endif // KERNEL_VARIABLES_HPP
