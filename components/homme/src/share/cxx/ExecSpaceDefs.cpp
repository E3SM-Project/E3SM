/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include <cassert>

#include <sstream>
#include <vector>

#include "ExecSpaceDefs.hpp"
#include "Dimensions.hpp"
#include "utilities/MathUtils.hpp"

namespace Homme {

// Since we're initializing from inside a Fortran code and don't have access to
// char** args to pass to Kokkos::initialize, we need to create and set the kokkos
// initialization settings manually.
void initialize_kokkos () {
  auto const settings = Kokkos::InitializationSettings()
    .set_map_device_id_by("mpi_rank")
    .set_disable_warnings(true);
  Kokkos::initialize(settings);
}

ThreadPreferences::ThreadPreferences ()
  : max_threads_usable(NP*NP),
    max_vectors_usable(NUM_PHYSICAL_LEV),
    prefer_threads(true),
    prefer_larger_team(true)
{}

namespace Parallel {

std::pair<int, int>
team_num_threads_vectors_from_pool (
  const int pool_size, const int num_parallel_iterations,
  const ThreadPreferences tp)
{
  assert(pool_size >= 1);
  assert(num_parallel_iterations >= 0);
  assert(tp.max_threads_usable >= 1 && tp.max_vectors_usable >= 1);

  const int num_threads_per_element = ( (num_parallel_iterations > 0 &&
                                         pool_size > num_parallel_iterations) ?
                                        pool_size / num_parallel_iterations :
                                        1 );

  const int num_threads = ( (tp.max_threads_usable > num_threads_per_element) ?
                            num_threads_per_element :
                            tp.max_threads_usable );

  const int p1 = num_threads, p2 = num_threads_per_element / num_threads;
  return tp.prefer_threads ? std::make_pair(p1, p2) : std::make_pair(p2, p1);
}

std::pair<int, int>
team_num_threads_vectors_for_gpu (
  const int num_warps_total, const int num_threads_per_warp,
  const int min_num_warps, const int max_num_warps,
  const int num_parallel_iterations, const ThreadPreferences tp)
{
  using Homme::nextpow2;
  using Homme::prevpow2;

  assert(num_warps_total > 0 && num_threads_per_warp > 0 && min_num_warps > 0);
  assert(num_parallel_iterations >= 0);
  assert(min_num_warps <= max_num_warps);
  assert(num_warps_total >= max_num_warps);
  assert(tp.max_threads_usable >= 1 && tp.max_vectors_usable >= 1);

#ifndef KOKKOS_ENABLE_SYCL
  int num_warps;
  if (tp.prefer_larger_team) {
    const int num_warps_usable =
      (tp.max_threads_usable *
       // Vector size is limited to threads/warp.
       std::min(num_threads_per_warp, tp.max_vectors_usable) +
       num_threads_per_warp - 1) /
      num_threads_per_warp;
    num_warps =
      std::max<int>( min_num_warps,
                     nextpow2(
                       std::min( max_num_warps,
                                 std::max( 1, num_warps_usable ))));
  } else {
    num_warps =
      // Min and max keep num_warps in bounds.
      std::max<int>( min_num_warps,
                     // We want num_warps to divide 32, the number of cores per SM, so
                     // apply nextpow2.
                     nextpow2(
                       std::min( max_num_warps,
                                 ( (num_parallel_iterations > 0 &&
                                    num_warps_total > num_parallel_iterations) ?
                                   (num_warps_total / num_parallel_iterations) :
                                   1 ))));
  }

  const int num_device_threads = num_warps * num_threads_per_warp;

  if (tp.prefer_threads) {
    const int num_threads = ( (tp.max_threads_usable > num_device_threads) ?
                              num_device_threads :
                              tp.max_threads_usable );

    return std::make_pair( num_threads,
                           prevpow2(num_device_threads / num_threads) );
  } else {
    const int num_vectors = prevpow2( (tp.max_vectors_usable > num_device_threads) ?
                                      num_device_threads :
                                      tp.max_vectors_usable );

    return std::make_pair( num_device_threads / num_vectors,
                           num_vectors );
  }
#else
  return std::make_pair(16,8);
#endif
}

} // namespace Parallel

std::pair<int, int>
DefaultThreadsDistribution<HommexxGPU>::
team_num_threads_vectors (const int num_parallel_iterations,
                          const ThreadPreferences tp) {
  // It appears we can't use Kokkos to tell us this. On current devices, using
  // fewer than 4 warps/thread block limits the thread occupancy to that
  // number/4. That seems to be in Cuda specs, but I don't know of a function
  // that provides this number. Use a configuration option that defaults to 4.
  int min_num_warps = HOMMEXX_CUDA_MIN_WARP_PER_TEAM;

  int max_num_warps = HOMMEXX_CUDA_MAX_WARP_PER_TEAM;
#ifdef KOKKOS_ENABLE_DEBUG
  // In debug builds, team size must be smaller because of Kokkos-side data.
  max_num_warps = std::min(max_num_warps, 8);
#endif

#ifdef KOKKOS_ENABLE_CUDA
  // No such thing sometime after Kokkos version 4.2:
  //    Kokkos::Impl::cuda_internal_maximum_concurrent_block_count();
  // This number is not super important as long as it's large. 3456 is the A100
  // spec, so we'll use that.
  const int num_warps_device = 3456;
  const int num_threads_warp = Kokkos::Impl::CudaTraits::WarpSize;
#elif defined(KOKKOS_ENABLE_HIP)
  // Use 64 wavefronts per CU and 120 CUs.
  const int num_warps_device = 120*64; // no such thing Kokkos::Impl::hip_internal_maximum_warp_count();
  const int num_threads_warp = Kokkos::Impl::HIPTraits::WarpSize;
#else
  // I want thread-distribution rules to be unit-testable even when GPU spaces
  // are off. Thus, make up a GPU-like machine:
  const int num_warps_device = 1792;
  const int num_threads_warp = 32;
  max_num_warps = 16;
#endif

  min_num_warps = std::min(min_num_warps, max_num_warps);

  return Parallel::team_num_threads_vectors_for_gpu(
    num_warps_device, num_threads_warp,
    min_num_warps, max_num_warps,
    num_parallel_iterations, tp);
}

} // namespace Homme
