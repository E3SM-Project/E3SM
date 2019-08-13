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

#ifdef KOKKOS_ENABLE_CUDA
# include <cuda.h>
#endif

namespace Homme {

// Since we're initializing from inside a Fortran code and don't have access to
// char** args to pass to Kokkos::initialize, we need to do some work on our
// own. As a side benefit, we'll end up running on GPU platforms optimally
// without having to specify --kokkos-ndevices on the command line.
void initialize_kokkos () {
  // This is in fact const char*, but Kokkos::initialize requires char*.
  std::vector<char*> args;

  //   This is the only way to get the round-robin rank assignment Kokkos
  // provides, as that algorithm is hardcoded in Kokkos::initialize(int& narg,
  // char* arg[]). Once the behavior is exposed in the InitArguments version of
  // initialize, we can remove this string code.
  //   If for some reason we're running on a GPU platform, have Cuda enabled,
  // but are using a different execution space, this initialization is still
  // OK. The rank gets a GPU assigned and simply will ignore it.
#ifdef KOKKOS_ENABLE_CUDA
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

ThreadPreferences::ThreadPreferences ()
  : max_threads_usable(NP*NP),
    max_vectors_usable(NUM_PHYSICAL_LEV),
    prefer_threads(true)
{}

namespace Parallel {

// Use a well-known trick for nextpow2, with a few tweaks for prevpow2;
// see, e.g.,
//   https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
unsigned short prevpow2 (unsigned short n) {
  if (n == 0) return 0;
  n >>= 1;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  ++n;
  return n;
}

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
  assert(num_warps_total > 0 && num_threads_per_warp > 0 && min_num_warps > 0);
  assert(num_parallel_iterations >= 0);
  assert(min_num_warps <= max_num_warps);
  assert(num_warps_total >= max_num_warps);
  assert(tp.max_threads_usable >= 1 && tp.max_vectors_usable >= 1);

  const int num_warps =
    std::max( min_num_warps,
              std::min( max_num_warps,
                        ( (num_parallel_iterations > 0 &&
                           num_warps_total > num_parallel_iterations) ?
                          (num_warps_total / num_parallel_iterations) :
                          1 )));

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
}

} // namespace Parallel

std::pair<int, int>
DefaultThreadsDistribution<Hommexx_Cuda>::
team_num_threads_vectors (const int num_parallel_iterations,
                          const ThreadPreferences tp) {
  // It appears we can't use Kokkos to tell us this. On current devices, using
  // fewer than 4 warps/thread block limits the thread occupancy to that
  // number/4. That seems to be in Cuda specs, but I don't know of a function
  // that provides this number. Use a configuration option that defaults to 4.
  const int min_num_warps = HOMMEXX_CUDA_MIN_WARP_PER_TEAM;
#ifdef KOKKOS_ENABLE_CUDA
  const int num_warps_device = Kokkos::Impl::cuda_internal_maximum_concurrent_block_count();
  const int num_threads_warp = Kokkos::Impl::CudaTraits::WarpSize;
  // The number on P100 is 16, but for some reason
  // Kokkos::Impl::cuda_internal_maximum_grid_count() returns 8. I may be
  // misusing the function. I have an open issue with the Kokkos team to resolve
  // this. For now:
# ifdef KOKKOS_ENABLE_DEBUG
  // Work around spurious "terminate called after throwing an instance
  // of 'std::runtime_error'" message.
  const int max_num_warps = 8;
# else
  const int max_num_warps = HOMMEXX_CUDA_MAX_WARP_PER_TEAM; //Kokkos::Impl::cuda_internal_maximum_grid_count();
# endif
#else
  // I want thread-distribution rules to be unit-testable even when Cuda is
  // off. Thus, make up a P100-like machine:
  const int num_warps_device = 1792;
  const int num_threads_warp = 32;
  const int max_num_warps = 16;
#endif

  return Parallel::team_num_threads_vectors_for_gpu(
    num_warps_device, num_threads_warp,
    min_num_warps, max_num_warps,
    num_parallel_iterations, tp);
}

} // namespace Homme
