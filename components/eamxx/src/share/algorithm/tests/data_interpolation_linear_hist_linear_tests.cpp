#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_config.hpp"

namespace scream {

TEST_CASE ("interpolation_linear_hist_linear")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  // Regardless of how EAMxx is configured, ignore leap years for this test
  set_use_leap_year(false);

  scorpio::init_subsystem(comm);

  auto data_grid   = create_point_grid("pg",   data_ngcols, data_nlevs, comm, 1);
  auto hfine_grid  = create_point_grid("pg_h", fine_ngcols, data_nlevs, comm, 1);
  auto vfine_grid  = create_point_grid("pg_v", data_ngcols, fine_nlevs, comm, 1);
  auto hvfine_grid = create_point_grid("pg_hv",fine_ngcols, fine_nlevs, comm, 1);

  auto time_interp_type = DataInterpolation::Linear;
  auto timeline = util::TimeLine::Linear;

  // We assume we start with t0 sometime between the first and second input slice.
  auto t_beg = get_first_slice_time();
  strvec_t files         = {"data_interpolation_0.nc",        "data_interpolation_1.nc",        "data_interpolation_2.nc"        };
  strvec_t files_no_ilev = {"data_interpolation_0_no_ilev.nc","data_interpolation_1_no_ilev.nc","data_interpolation_2_no_ilev.nc"};

  SECTION ("no-horiz") {
    SECTION ("no-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=NO ..........\n");
      run_tests(data_grid,  files,         t_beg, timeline, time_interp_type);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=NO .......... PASS\n");
    }
    SECTION ("p1d-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=p1d .........\n");
      run_tests(vfine_grid, files_no_ilev, t_beg, timeline, time_interp_type, P1D);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=p1d ......... PASS\n");
    }
    SECTION ("p2d-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=p2d .........\n");
      run_tests(vfine_grid, files_no_ilev, t_beg, timeline, time_interp_type, P2D);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=p2d ......... PASS\n");
    }
    SECTION ("p3d-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=p3d .........\n");
      run_tests(vfine_grid, files_no_ilev, t_beg, timeline, time_interp_type, P3D);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=NO,  vert_remap=p3d ......... PASS\n");
    }
  }

  SECTION ("yes-horiz") {
    SECTION ("no-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=NO ..........\n");
      run_tests(hfine_grid, files,         t_beg, timeline, time_interp_type);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=NO .......... PASS\n");
    }
    SECTION ("p1d-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=p1d .........\n");
      run_tests(hvfine_grid,files_no_ilev, t_beg, timeline, time_interp_type, P1D);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=p1d ......... PASS\n");
    }
    SECTION ("p2d-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=p2d .........\n");
      run_tests(hvfine_grid,files_no_ilev, t_beg, timeline, time_interp_type, P2D);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=p2d ......... PASS\n");
    }
    SECTION ("p3d-vert") {
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=p3d .........\n");
      run_tests(hvfine_grid,files_no_ilev, t_beg, timeline, time_interp_type, P3D);
      root_print(comm,"  interp=LINEAR, timeline=LINEAR,   horiz_remap=YES, vert_remap=p3d ......... PASS\n");
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
