#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_config.hpp"

namespace scream {

// data_interpolation_static.nc      – all fields not time-dep, delta = 0
// data_interpolation_static_no_ilev.nc – same but no ilev coord (for vremap)
// data_interpolation_0_no_ilev.nc   – time-dep file; we pick time_index=3 (Oct),
//   delta = delta_data[9] = 90.

TEST_CASE ("static_data")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  // Regardless of how EAMxx is configured, ignore leap years for this test
  set_use_leap_year(false);

  scorpio::init_subsystem(comm);

  auto data_grid   = create_point_grid("pg",   data_ngcols, data_nlevs, comm, 1);
  auto hfine_grid  = create_point_grid("pg_h", fine_ngcols, data_nlevs, comm, 1);
  auto vfine_grid  = create_point_grid("pg_v", data_ngcols, fine_nlevs, comm, 1);
  auto hvfine_grid = create_point_grid("pg_hv",fine_ngcols, fine_nlevs, comm, 1);

  SECTION ("no-time-dim") {
    // File has no time dimension; all fields are constant (delta = 0)
    SECTION ("no-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=NO ..........\n");
        run_static_tests(data_grid,  "data_interpolation_static.nc",        -1, 0.0);
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=p1d .........\n");
        run_static_tests(vfine_grid, "data_interpolation_static_no_ilev.nc", -1, 0.0, P1D);
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=p2d .........\n");
        run_static_tests(vfine_grid, "data_interpolation_static_no_ilev.nc", -1, 0.0, P2D);
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=p3d .........\n");
        run_static_tests(vfine_grid, "data_interpolation_static_no_ilev.nc", -1, 0.0, P3D);
        root_print(comm,"  static, no-time-dim, horiz_remap=NO,  vert_remap=p3d ......... PASS\n");
      }
    }
    SECTION ("yes-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=NO ..........\n");
        run_static_tests(hfine_grid, "data_interpolation_static.nc",        -1, 0.0);
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=p1d .........\n");
        run_static_tests(hvfine_grid,"data_interpolation_static_no_ilev.nc", -1, 0.0, P1D);
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=p2d .........\n");
        run_static_tests(hvfine_grid,"data_interpolation_static_no_ilev.nc", -1, 0.0, P2D);
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=p3d .........\n");
        run_static_tests(hvfine_grid,"data_interpolation_static_no_ilev.nc", -1, 0.0, P3D);
        root_print(comm,"  static, no-time-dim, horiz_remap=YES, vert_remap=p3d ......... PASS\n");
      }
    }
  }

  SECTION ("specific-snapshot") {
    // File has a time dimension; we select time_index=3 (October slice).
    // delta_data at that slice: mm=3 in data_interpolation_0, mm_index=9 → delta_data[9]=90
    constexpr double oct_delta = 90.0;
    SECTION ("no-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=NO ..........\n");
        run_static_tests(data_grid,  "data_interpolation_0.nc",        3, oct_delta);
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=p1d .........\n");
        run_static_tests(vfine_grid, "data_interpolation_0_no_ilev.nc", 3, oct_delta, P1D);
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=p2d .........\n");
        run_static_tests(vfine_grid, "data_interpolation_0_no_ilev.nc", 3, oct_delta, P2D);
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=p3d .........\n");
        run_static_tests(vfine_grid, "data_interpolation_0_no_ilev.nc", 3, oct_delta, P3D);
        root_print(comm,"  static, snapshot[3], horiz_remap=NO,  vert_remap=p3d ......... PASS\n");
      }
    }
    SECTION ("yes-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=NO ..........\n");
        run_static_tests(hfine_grid, "data_interpolation_0.nc",        3, oct_delta);
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=p1d .........\n");
        run_static_tests(hvfine_grid,"data_interpolation_0_no_ilev.nc", 3, oct_delta, P1D);
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=p2d .........\n");
        run_static_tests(hvfine_grid,"data_interpolation_0_no_ilev.nc", 3, oct_delta, P2D);
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=p3d .........\n");
        run_static_tests(hvfine_grid,"data_interpolation_0_no_ilev.nc", 3, oct_delta, P3D);
        root_print(comm,"  static, snapshot[3], horiz_remap=YES, vert_remap=p3d ......... PASS\n");
      }
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
