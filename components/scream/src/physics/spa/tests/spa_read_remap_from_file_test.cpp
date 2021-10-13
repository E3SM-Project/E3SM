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
struct UnitWrap::UnitTest<D>::TestReadRemapData {

  template <typename S>
  using view_1d = typename KokkosTypes<D>::template view_1d<S>;

  static void run_property()
  {
//    // Set up the mpi communicator and init the pio subsystem
//    ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
//    MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
//    scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

    // Establish the SPA function object
    using SPAFunc         = spa::SPAFunctions<Real, DefaultDevice>;
    SPAFunc::SPAHorizInterp spa_horiz_interp;

//    const std::string remap_file_name = "map_ne2np4_to_ne4np4_mono.nc";
    const std::string remap_file_name = "spa_data_for_testing.nc";

    Int target_grid_ncols = 2;
    view_1d<int> dofs_gids("",target_grid_ncols);
    for (int ii=0;ii<ncols;ii++) {
      dofs_gids(ii) = ii;
    }
    SPAFunc::get_remap_weights_from_file(remap_file_name,target_grid_ncols,dofs_gids,spa_horiz_interp);

    REQUIRE(spa_horiz_interp.length==6); // The test file has 17,498 points in it
    REQUIRE(spa_horiz_interp.source_grid_ncols==3); // The test file maps ne2 (218 cols) to ne4 (866 cols)

    // We have a few metrics to ensure that the data read from file matches the data in the file.
    Real tol = 1e5*std::numeric_limits<Real>::epsilon();
    Int col_sum = 0;
    Int row_sum = 0;
    Real wgt_sum = 0.0;
    view_1d<Real> wgts("",target_grid_ncols);
    Kokkos::deep_copy(wgts,0.0);
    for (int i=0; i<spa_horiz_interp.length; i++) {
      col_sum += spa_horiz_interp.target_grid_loc[i];
      row_sum += spa_horiz_interp.source_grid_loc[i];
      wgt_sum += spa_horiz_interp.weights[i];
      wgts(spa_horiz_interp.target_grid_loc[i]) += spa_horiz_interp.weights[i];
    }
    REQUIRE(col_sum == 9);
    REQUIRE(row_sum == 12);
    REQUIRE(std::abs(wgt_sum - 2.0) < tol);
    // The sum of remap weights should always be 1.0
    for (int i=0; i<target_grid_ncols; i++) {
      REQUIRE(wgts[i]==1.0);
    }
    
//    // All Done 
//    scorpio::eam_pio_finalize();
  }

}; // struct TestReadRemapData

} // namespace unit_test
} // namespace spa
} //namespace scream

namespace {

TEST_CASE("spa_read_remap_test", "spa")
{
  using TestStruct = scream::spa::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestReadRemapData;

  TestStruct::run_property();
}

}
