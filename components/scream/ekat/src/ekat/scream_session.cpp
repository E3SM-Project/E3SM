#include "ekat/scream_session.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_arch.hpp"

#include <vector>

#ifdef SCREAM_FPE
# include <xmmintrin.h>
#endif

namespace scream_impl {
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

} // anonymous namespace

namespace scream {

void initialize_scream_session () {
  enable_default_fpes ();

  if (!Kokkos::is_initialized()) {
    printf("Calling initialize_kokkos\n");
    scream_impl::initialize_kokkos();
  }

  std::cout << util::config_string() << "\n";
}

void initialize_scream_session (int argc, char **argv) {
  enable_default_fpes ();
  Kokkos::initialize(argc, argv);
  std::cout << util::config_string() << "\n";
}
extern "C" {
void finalize_scream_session () {
  Kokkos::finalize();
}
} // extern "C"

} // namespace scream
