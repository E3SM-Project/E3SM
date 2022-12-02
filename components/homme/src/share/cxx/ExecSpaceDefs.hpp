/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_EXEC_SPACE_DEFS_HPP
#define HOMMEXX_EXEC_SPACE_DEFS_HPP

#include <cassert>

#include <Kokkos_Core.hpp>

#include "Config.hpp"
#include "Dimensions.hpp"
#include "vector/vector_pragmas.hpp"

namespace Homme
{

// Some in-house names for Kokkos exec spaces, which are
// always defined, possibly as alias of void

#ifdef HOMMEXX_ENABLE_GPU

#ifdef KOKKOS_ENABLE_CUDA
using HommexxGPU = Kokkos::Cuda;
#endif

#ifdef KOKKOS_ENABLE_HIP
using HommexxGPU = Kokkos::Experimental::HIP;
#endif

#else
using HommexxGPU = void;
#endif



#ifdef KOKKOS_ENABLE_OPENMP
using Hommexx_OpenMP = Kokkos::OpenMP;
#else
using Hommexx_OpenMP = void;
#endif

#ifdef KOKKOS_ENABLE_PTHREADS
using Hommexx_Threads = Kokkos::Threads;
#else
using Hommexx_Threads = void;
#endif

#ifdef KOKKOS_ENABLE_SERIAL
using Hommexx_Serial = Kokkos::Serial;
#else
using Hommexx_Serial = void;
#endif

#ifdef HOMMEXX_ENABLE_GPU
# define HOMMEXX_STATIC
#else
# define HOMMEXX_STATIC static
#endif

// Selecting the execution space. If no specific request, use Kokkos default
// exec space
#ifdef HOMMEXX_ENABLE_GPU
using ExecSpace = HommexxGPU;
#elif defined(HOMMEXX_OPENMP_SPACE)
using ExecSpace = Hommexx_OpenMP;
#elif defined(HOMMEXX_THREADS_SPACE)
using ExecSpace = Hommexx_Threads;
#elif defined(HOMMEXX_SERIAL_SPACE)
using ExecSpace = Hommexx_Serial;
#elif defined(HOMMEXX_DEFAULT_SPACE)
using ExecSpace = Kokkos::DefaultExecutionSpace::execution_space;
#else
#error "No valid execution space choice"
#endif // HOMMEXX_EXEC_SPACE

static_assert (!std::is_same<ExecSpace,void>::value,
               "Error! You are trying to use an ExecutionSpace not enabled in Kokkos.\n");

template <typename ExeSpace>
struct OnGpu { enum : bool { value = false }; };

template <>
struct OnGpu<HommexxGPU> { enum : bool { value = true }; };

// Call this instead of Kokkos::initialize.
void initialize_kokkos();

// What follows provides utilities to parameterize the parallel machine (CPU/KNL
// cores within a rank, GPU attached to a rank) optimally. The parameterization
// is a nontrivial function of available resources, number of parallel
// iterations to be performed, and kernel-specific preferences regarding team
// and vector dimensions of parallelization.
//   So far, we are able to hide the details inside a call like this:
//     Homme::get_default_team_policy<ExecSpace>(data.num_elems * data.qsize);
// thus, all that follows except the function get_default_team_policy may not
// need to be used except in the implementation of get_default_team_policy.
//   If that remains true forever, we can move all of this code to
// ExecSpaceDefs.cpp.

// Preferences to guide distribution of physical threads among team and vector
// dimensions. Default values are sensible.
struct ThreadPreferences {
  // Max number of threads a kernel can use. Default: NP*NP.
  int max_threads_usable;
  // Max number of vectors a kernel can use. Default: NUM_PHYSICAL_LEV.
  int max_vectors_usable;
  // Prefer threads to vectors? Default: true.
  bool prefer_threads;
  // Prefer larger teams even when league_size is large.
  bool prefer_larger_team;

  ThreadPreferences();
};

namespace Parallel {
// Determine (#threads, #vectors) as a function of a pool of threads provided to
// the process and the number of parallel iterations to perform.
std::pair<int, int>
team_num_threads_vectors_from_pool(
  const int pool_size, const int num_parallel_iterations,
  const ThreadPreferences tp = ThreadPreferences());

// Determine (#threads, #vectors) as a function of a pool of warps provided to
// the process, the number of threads per warp, the maximum number of warps a
// team can use, and the number of parallel iterations to perform.
std::pair<int, int>
team_num_threads_vectors_for_gpu(
  const int num_warps_total, const int num_threads_per_warp,
  const int min_num_warps, const int max_num_warps,
  const int num_parallel_iterations,
  const ThreadPreferences tp = ThreadPreferences());
} // namespace Parallel

// Device-dependent distribution of physical threads over teams and vectors. The
// general case is for a machine with a pool of threads, like KNL and CPU.
template <typename ExecSpaceType>
struct DefaultThreadsDistribution {
  static std::pair<int, int>
  team_num_threads_vectors(const int num_parallel_iterations,
                           const ThreadPreferences tp = ThreadPreferences()) {
    return Parallel::team_num_threads_vectors_from_pool(
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      ExecSpaceType::thread_pool_size()
#else
      ExecSpaceType::impl_thread_pool_size()
#endif
      , num_parallel_iterations, tp);
  }
};

// Specialization for a GPU, where threads can't be viewed as existing simply in
// a pool.
template <>
struct DefaultThreadsDistribution<HommexxGPU> {
  static std::pair<int, int>
  team_num_threads_vectors(const int num_parallel_iterations,
                           const ThreadPreferences tp = ThreadPreferences());
};

// Return a TeamPolicy using defaults that, so far, have been good for all use
// cases. Use of this function means you don't have to use
// DefaultThreadsDistribution.
template <typename ExecSpace, typename... Tags>
Kokkos::TeamPolicy<ExecSpace, Tags...>
get_default_team_policy(const int num_parallel_iterations,
                        const ThreadPreferences tp = ThreadPreferences()) {
  const auto threads_vectors =
    DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
      num_parallel_iterations, tp);
  auto policy = Kokkos::TeamPolicy<ExecSpace, Tags...>(num_parallel_iterations,
                                                   threads_vectors.first,
                                                   threads_vectors.second);
  policy.set_chunk_size(1);
  return policy;
}

template<typename ExecSpaceType, typename... Tags>
static
typename std::enable_if<!OnGpu<ExecSpaceType>::value,int>::type
get_num_concurrent_teams (const Kokkos::TeamPolicy<ExecSpaceType,Tags...>& policy) {
  const int team_size = policy.team_size();
  const int concurrency = ExecSpaceType::concurrency();
  return (concurrency + team_size - 1) / team_size;
}

template<typename ExecSpaceType, typename... Tags>
static
typename std::enable_if<OnGpu<ExecSpaceType>::value,int>::type
get_num_concurrent_teams (const Kokkos::TeamPolicy<ExecSpaceType,Tags...>& policy) {
  // const int team_size = policy.team_size() * policy.vector_length();
  // const int concurrency = ExecSpaceType::concurrency();
  // return (concurrency + team_size - 1) / team_size;
  return policy.league_size();
}

template <typename ExecSpace, typename... Tags>
static int get_num_concurrent_teams(const int num_parallel_iterations) {
  return get_num_concurrent_teams(
      get_default_team_policy<ExecSpace, Tags...>(num_parallel_iterations));
}

// A templated typedef for MD range policy (used in RK stages)
template<typename ExecutionSpace, int Rank>
using MDRangePolicy = Kokkos::MDRangePolicy
                          < ExecutionSpace,
                            Kokkos::Rank
                              < Rank,
                                Kokkos::Iterate::Right,
                                Kokkos::Iterate::Right
                              >,
                            Kokkos::IndexType<int>
                          >;

template <typename ExeSpace>
struct Memory {
  enum : bool { on_gpu = OnGpu<ExeSpace>::value };

  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION static
  Scalar* get_shmem (const typename Kokkos::TeamPolicy<ExeSpace>::member_type&,
                     const size_t sz = 0) {
    return nullptr;
  }

  template <typename Scalar, int N>
  class AutoArray {
    Scalar data_[N];
  public:
    KOKKOS_INLINE_FUNCTION AutoArray (Scalar*) {}
    KOKKOS_INLINE_FUNCTION Scalar& operator[] (const int& i) {
      assert(i >= 0);
      assert(i < N);
      return data_[i];
    }
    KOKKOS_INLINE_FUNCTION Scalar* data () { return data_; }
  };
};

template <>
struct Memory<HommexxGPU> {
  enum : bool { on_gpu = OnGpu<HommexxGPU>::value };

  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION static
  Scalar* get_shmem (const Kokkos::TeamPolicy<HommexxGPU>::member_type& team,
                     const size_t n = 0) {
    return static_cast<Scalar*>(team.team_shmem().get_shmem(n*sizeof(Scalar)));
  }

  template <typename Scalar, int N>
  class AutoArray {
    Scalar* data_;
  public:
    KOKKOS_INLINE_FUNCTION AutoArray (Scalar* data) : data_(data) {}
    KOKKOS_INLINE_FUNCTION Scalar& operator[] (const int& i) {
      assert(i >= 0);
      assert(i < N);
      return data_[i];
    }
    KOKKOS_INLINE_FUNCTION Scalar* data () { return data_; }
  };
};

// A templated namespace to hold variations on intra-team parallel
// dispatches.
template <typename ExeSpace=ExecSpace>
struct Dispatch {
  // Match the HOMMEXX_BFB_TESTING function in the Cuda
  // specialization.
  template<typename LoopBdyType, class Lambda, typename ValueType>
  static KOKKOS_FORCEINLINE_FUNCTION
  void parallel_reduce (
    const Kokkos::TeamPolicy<ExecSpace>::member_type& team,
    const LoopBdyType& loop_boundaries,
    const Lambda& lambda, ValueType& result)
  {
    Kokkos::parallel_reduce(loop_boundaries, lambda, result);
  }

  // Improve performance in two ways: NP*NP is a compile-time range
  // limit, and use pragma simd. This substantially speeds up the
  // limiters, for example.
  template<class Lambda>
  static KOKKOS_FORCEINLINE_FUNCTION
  void parallel_for_NP2 (
    const typename Kokkos::TeamPolicy<ExeSpace>::member_type& team,
    const Lambda& lambda)
  {
VECTOR_SIMD_LOOP
    for (int k = 0; k < NP*NP; ++k)
      lambda(k);
  }

  template<class Lambda, typename ValueType>
  static KOKKOS_FORCEINLINE_FUNCTION
  void parallel_reduce_NP2 (
    const typename Kokkos::TeamPolicy<ExeSpace>::member_type& team,
    const Lambda& lambda, ValueType& result)
  {
    result = Kokkos::reduction_identity<ValueType>::sum();
    for (int k = 0; k < NP*NP; ++k)
      lambda(k, result);
  }

  template<class Lambda>
  static KOKKOS_FORCEINLINE_FUNCTION
  void parallel_scan (
    const typename Kokkos::TeamPolicy<ExeSpace>::member_type& team,
    const int num_iters,
    const Lambda& lambda)
  {
    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,num_iters), lambda);
  }
};

#ifdef HOMMEXX_ENABLE_GPU
template <>
struct Dispatch<HommexxGPU> {
  using ExeSpace = HommexxGPU;

  template<typename LoopBdyType, class Lambda, typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION
  static void parallel_reduce (
    const Kokkos::TeamPolicy<ExeSpace>::member_type& team,
    const LoopBdyType& loop_boundaries,
    const Lambda& lambda, ValueType& result)
  {
#ifdef HOMMEXX_BFB_TESTING
    // We want to get C++ on GPU to match F90 on CPU. Thus, need to
    // serialize parallel reductions.

    // All threads init result.
    // NOTE: we *need* an automatic temporary, since we do not know
    //       where 'result' comes from. If result itself is an automatic
    //       variable, using it would be fine (only one vector lane would
    //       actually have a nonzero value after the single). But if
    //       result is taken from a view, then all vector lanes would
    //       see the updated value before the vector_reduce call,
    //       which will cause the final answer to be multiplied by the
    //       size of the warp.
    auto local_tmp = Kokkos::reduction_identity<ValueType>::sum();
    // One thread sums.
    Kokkos::single(Kokkos::PerThread(team), [&] () {
        for (auto i = loop_boundaries.start; i < loop_boundaries.end; ++i)
          lambda(i, local_tmp);
      });
    // Broadcast result to all threads by doing sum of one thread's
    // non-0 value and the rest of the 0s.
    Kokkos::Impl::CudaTeamMember::vector_reduce(
      Kokkos::Sum<ValueType>(local_tmp));
    result = local_tmp;
#else
    Kokkos::parallel_reduce(loop_boundaries, lambda, result);
#endif
  }

  // Match the performance-improving impls in the non-GPU impl.
  template<class Lambda>
  static KOKKOS_FORCEINLINE_FUNCTION
  void parallel_for_NP2 (
    const Kokkos::TeamPolicy<HommexxGPU>::member_type& team,
    const Lambda& lambda)
  {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, NP*NP), lambda);
  }

  template<class Lambda, typename ValueType>
  static KOKKOS_FORCEINLINE_FUNCTION
  void parallel_reduce_NP2 (
    const Kokkos::TeamPolicy<HommexxGPU>::member_type& team,
    const Lambda& lambda, ValueType& result)
  {
    parallel_reduce(team, Kokkos::ThreadVectorRange(team, NP*NP),
                    lambda, result);
  }

  template<class Lambda>
  static KOKKOS_FORCEINLINE_FUNCTION
  void parallel_scan (
    const typename Kokkos::TeamPolicy<ExeSpace>::member_type& team,
    const int num_iters,
    const Lambda& lambda)
  {
#ifdef HOMMEXX_BFB_TESTING
    // We want to get C++ on GPU to match F90 on CPU. Thus, need to
    // serialize parallel scans.

    // Detect the value type
    using value_type =
      typename Kokkos::Impl::FunctorAnalysis
        < Kokkos::Impl::FunctorPatternInterface::SCAN
        , void
        , Lambda >::value_type ;

    // All threads init result.
    value_type accumulator = Kokkos::reduction_identity<value_type>::sum();
    // Only one thread does the work, i.e., only one sweeps (so last arg to lambda is true)
    Kokkos::single(Kokkos::PerThread(team), [&] () {
      for (int i = 0; i < num_iters; ++i) {
        lambda(i, accumulator, true);
      }
    });
#else
    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,num_iters), lambda);
#endif
  }
};
#endif

#ifdef KOKKOS_ENABLE_CUDA
// Cuda-provided GPU-safe replacements for std functions.
using ::isnan;
#else
using std::isnan;
#endif

} // namespace Homme

#endif // HOMMEXX_EXEC_SPACE_DEFS_HPP
