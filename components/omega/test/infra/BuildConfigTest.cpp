//===-- Test driver for OMEGA build configuration ----------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA build configuration
///
/// This driver checks that the architecture Omega was configured with agrees
/// with the Kokkos backend it was actually built against. OMEGA_ARCH is chosen
/// by the build system - from the CIME machine settings in an E3SM build, or
/// from the standalone configuration - and it drives both the OMEGA_ENABLE_*
/// and OMEGA_TARGET_DEVICE macros and the Kokkos backend selection. Nothing at
/// runtime re-derives it, so a mismatch between the two is silent: Omega would
/// believe it is running on one architecture while Kokkos executes on another.
/// That is the failure this test exists to catch, in particular a host-only
/// Kokkos paired with a GPU OMEGA_ARCH or the reverse.
///
//
//===-----------------------------------------------------------------------===/

#include "DataTypes.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "mpi.h"

#include <string>
#include <type_traits>

using namespace OMEGA;

// OMEGA_ARCH is defined as a bare token (e.g. -DOMEGA_ARCH=CUDA) so it must be
// stringified before it can be compared.
#define OMEGA_STR_HELPER(Arg) #Arg
#define OMEGA_STR(Arg)        OMEGA_STR_HELPER(Arg)

//------------------------------------------------------------------------------
// Returns the Kokkos execution space name expected for a given OMEGA_ARCH.
// These are the names Kokkos itself reports, so the comparison stays valid
// across Kokkos versions that rename the underlying types.

std::string expectedKokkosSpace(const std::string &Arch // [in] OMEGA_ARCH value
) {

   if (Arch == "CUDA")
      return "Cuda";
   if (Arch == "HIP")
      return "HIP";
   if (Arch == "SYCL")
      return "SYCL";
   if (Arch == "OPENMP")
      return "OpenMP";
   if (Arch == "SERIAL")
      return "Serial";

   return "";

} // end expectedKokkosSpace

//------------------------------------------------------------------------------
// The test driver for the build configuration.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   // These are needed to set up logging for output
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   initLogging(DefEnv);
   LOG_INFO("----- Build Configuration Unit Testing -----");

   {
      const std::string ArchName = OMEGA_STR(OMEGA_ARCH);
      const std::string KokkosName{Kokkos::DefaultExecutionSpace::name()};

      LOG_INFO("BuildConfigTest: OMEGA_ARCH = {}", ArchName);
      LOG_INFO("BuildConfigTest: Kokkos default execution space = {}",
               KokkosName);

      // OMEGA_ARCH must be one of the architectures the build system knows how
      // to configure. An empty or unrecognized value means the build system
      // failed to determine one.
      const std::string ExpectedSpace = expectedKokkosSpace(ArchName);
      if (ExpectedSpace.empty())
         ABORT_ERROR("BuildConfigTest: FAIL OMEGA_ARCH '{}' is not one of "
                     "CUDA, HIP, SYCL, OPENMP, SERIAL",
                     ArchName);

      // The matching OMEGA_ENABLE_<ARCH> macro must have been defined, since
      // the source uses it to select architecture-specific code paths.
      std::string EnabledArch;
#ifdef OMEGA_ENABLE_CUDA
      EnabledArch = "CUDA";
#endif
#ifdef OMEGA_ENABLE_HIP
      EnabledArch = "HIP";
#endif
#ifdef OMEGA_ENABLE_SYCL
      EnabledArch = "SYCL";
#endif
#ifdef OMEGA_ENABLE_OPENMP
      EnabledArch = "OPENMP";
#endif
#ifdef OMEGA_ENABLE_SERIAL
      EnabledArch = "SERIAL";
#endif

      if (EnabledArch != ArchName)
         ABORT_ERROR("BuildConfigTest: FAIL OMEGA_ARCH is '{}' but the macro "
                     "OMEGA_ENABLE_{} was not the one defined (got '{}')",
                     ArchName, ArchName, EnabledArch);

      // Kokkos must actually be running on the architecture Omega was told it
      // would use. This is the check that catches a host-only Kokkos in a build
      // that believes it is on a GPU, or the reverse.
      if (KokkosName != ExpectedSpace)
         ABORT_ERROR("BuildConfigTest: FAIL OMEGA_ARCH '{}' expects Kokkos "
                     "execution space '{}' but Kokkos is using '{}'",
                     ArchName, ExpectedSpace, KokkosName);

      // OMEGA_TARGET_DEVICE must be set for, and only for, a build whose
      // default execution space has its own memory space. Omega uses it to
      // decide whether device arrays need explicit host mirrors.
      constexpr bool KokkosIsDevice =
          !std::is_same_v<Kokkos::DefaultExecutionSpace::memory_space,
                          Kokkos::HostSpace>;
#ifdef OMEGA_TARGET_DEVICE
      constexpr bool OmegaIsDevice = true;
#else
      constexpr bool OmegaIsDevice = false;
#endif

      if (OmegaIsDevice != KokkosIsDevice)
         ABORT_ERROR("BuildConfigTest: FAIL OMEGA_TARGET_DEVICE is {} but the "
                     "Kokkos default execution space '{}' is {} device space",
                     OmegaIsDevice ? "defined" : "not defined", KokkosName,
                     KokkosIsDevice ? "a" : "not a");

      LOG_INFO("BuildConfigTest: OMEGA_TARGET_DEVICE = {}",
               OmegaIsDevice ? "defined" : "not defined");
   }

   LOG_INFO("----- Build Configuration Unit Tests Successful -----");
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0; // if we made it here, return successfully

} // end of main
//===-----------------------------------------------------------------------===/
