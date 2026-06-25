#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_config.hpp"

namespace scream {

TEST_CASE ("exceptions")
{
  // Test correctness of some exception handling inside the DataInterpolation source code
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);
  auto grid = create_point_grid("pg",data_ngcols,data_nlevs,comm,1);

  auto fields = create_fields(grid,false);

  REQUIRE_THROWS (create_interp(nullptr,fields)); // Invalid grid pointer

  auto interp = create_interp(grid,fields);

  strvec_t files = {"/etc/shadow"};
  REQUIRE_THROWS (interp->setup_linear_time_database(files)); // Input file not readable

  interp->setup_linear_time_database({"./data_interpolation_0.nc"});
  util::TimeStamp t0 ({2000,1,1},{0,0,0});
  REQUIRE_THROWS (interp->init_time_interpolation(t0,DataInterpolation::Linear)); // linear timeline, but t0<first_slice

  util::TimeStamp t1 ({2020,1,1},{0,0,0});
  REQUIRE_THROWS (interp->init_time_interpolation(t1,DataInterpolation::Linear)); // linear timeline, but t0>last_slice

  scorpio::finalize_subsystem();
}

} // namespace scream
