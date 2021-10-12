#include "catch2/catch.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "spa_unit_tests_common.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace spa {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::InitPIO {

  template <typename S>
  using view_1d = typename KokkosTypes<D>::template view_1d<S>;

  static void run()
  {
    // Set up the mpi communicator and init the pio subsystem
    ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
    MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
    scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler
   
    REQUIRE(true);
 
    // All Done 
//    scorpio::eam_pio_finalize();
  }

}; // struct InitPIO

} // namespace unit_test
} // namespace spa
} //namespace scream

namespace {

TEST_CASE("spa_init_pio", "spa")
{
  using TestStruct = scream::spa::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::InitPIO;

  TestStruct::run();
}

}
