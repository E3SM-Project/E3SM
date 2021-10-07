#include "catch2/catch.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "spa_unit_tests_common.hpp"
#include "physics/spa/spa_functions.hpp"

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
struct UnitWrap::UnitTest<D>::TestReadDataFile {

  template <typename S>
  using view_1d = typename KokkosTypes<D>::template view_1d<S>;

  static void run_property()
  {
    // Set up the mpi communicator and init the pio subsystem
    ekat::Comm spa_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
    MPI_Fint fcomm = MPI_Comm_c2f(spa_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
    scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

    // Establish the SPA function object
    using SPAFunc         = spa::SPAFunctions<Real, DefaultDevice>;

    std::string spa_data_file = "spa_data_for_testing.nc";
    Int time_index = 1;
    Int ncols    = 2;
    Int nlevs    = 4;
    Int nswbands = 2;
    Int nlwbands = 3;
    SPAFunc::SPAHorizInterp spa_horiz_interp(10);
    SPAFunc::SPAData data_beg(ncols, nlevs, nswbands, nlwbands);
    SPAFunc::SPAData data_end(ncols, nlevs, nswbands, nlwbands);

    // Construct the horiz_interp data on the fly
    spa_horiz_interp.source_grid_ncols = 1;
    
    for (int ii=0;ii<spa_horiz_interp.length;ii++)
    {
      spa_horiz_interp.weights(ii) = Real(ii+1)/Real(ii+2);
      spa_horiz_interp.source_grid_loc(ii) = ii % spa_horiz_interp.source_grid_ncols + 1;
      spa_horiz_interp.target_grid_loc(ii) = ii % ncols + 1;
    }

    SPAFunc::update_spa_data_from_file(spa_comm,spa_data_file, time_index, nswbands, nlwbands,
                                       spa_horiz_interp, data_beg);
    SPAFunc::update_spa_data_from_file(spa_comm,spa_data_file, time_index+1, nswbands, nlwbands,
                                       spa_horiz_interp, data_end);
    
  }

}; // struct TestReadDataFile

} // namespace unit_test
} // namespace spa
} //namespace scream

namespace {

TEST_CASE("spa_read_data_from_file_test", "spa")
{
  using TestStruct = scream::spa::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestReadDataFile;

  TestStruct::run_property();
}

}
