#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_config.hpp"

namespace scream {

// Test that setup_yearly_periodic_time_database with a start_ts correctly rotates
// the slice array, so that the 12-months period begins at start_ts.
//
// The input data are 4 files, each with 6 months. We do two runs:
//  - run 1: start_ts=Sep of 1st year, which should pick data Sep01-Aug02
//  - run 2: start_ts=Sep of 2nd year, which should pick data Sep02-Dec02-Jan02-Aug02

TEST_CASE ("start_ts_rotation")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  // Regardless of how EAMxx is configured, ignore leap years for this test
  set_use_leap_year(false);

  scorpio::init_subsystem(comm);

  auto data_grid = create_point_grid("pg", data_ngcols, data_nlevs, comm, 1);

  strvec_t files = {
    "data_interpolation_0.nc",
    "data_interpolation_1.nc",
    "data_interpolation_2.nc",
    "data_interpolation_3.nc"
  };

  auto first_slice = get_first_slice_time();
  auto sep01 = first_slice;
  for (int mm=0; mm<8; ++mm) {
    sep01 += spd * sep01.days_in_curr_month();
  }
  auto sep02 = sep01 + spd*365;

  int dt = 20*spd;
  for (auto start_ts : {sep01,sep02}) {
    std::cout << "--- start_ts: " << start_ts.to_string() << " ---\n";
    // t0 = 5 days after start_ts, so the interval is [Sep,Oct]
    // auto t0 = start_ts + 5*spd;
    auto t0 = start_ts + 0*spd;

    auto fields = create_fields(data_grid, false);
    fields.pop_back(); // Don't interpolate p1d

    auto interp = create_interp(data_grid, fields);
    interp->setup_periodic_time_database(files, start_ts);     // start_ts: rotate so Oct is first slice
    interp->create_horiz_remappers("");  // no hremap (data_grid == src_grid)
    interp->create_vert_remapper();      // no vremap (VRemapType::None)
    interp->init_time_interpolation(t0,DataInterpolation::Linear);

    auto base     = create_fields(data_grid, true);
    auto ones     = create_fields(data_grid, false);
    auto diff     = create_fields(data_grid, false);
    auto expected = create_fields(data_grid, false);
    for (auto& f : ones) { f.deep_copy(1); }

    int nfields = static_cast<int>(fields.size());

    auto t_beg = start_ts;
    auto t_end = t_beg + spd*t_beg.days_in_curr_month();
    auto time = t0;
    for (int istep=0; istep<4; ++istep) {

      interp->run(time);

      // Check
      auto alpha = static_cast<double>(time-t_beg) / (t_end-t_beg);

      auto mm_beg = t_beg.get_month()-1;
      auto mm_end = t_end.get_month()-1;
      if (t_beg.get_year()!=first_slice.get_year())
        mm_beg += 12;
      if (t_end.get_year()!=first_slice.get_year())
        mm_end += 12;
      double delta = delta_data[mm_beg]*(1-alpha) + delta_data[mm_end]*alpha;

      for (int i=0; i<nfields; ++i) {
        expected[i].update(base[i], 1, 0);
        expected[i].update(ones[i], delta, 1.0);

        diff[i].update(fields[i],   1, 1);
        diff[i].update(expected[i], -1, 1);
        diff[i].scale(1.0 / frobenius_norm(expected[i]).as<Real>());
        REQUIRE(frobenius_norm(diff[i]).as<Real>() < tol);
      }

      // Udpate time and, if neede, the time interval for manual calc
      time += dt;
      if (t_end<time) {
        t_beg = t_end;
        t_end += spd*t_end.days_in_curr_month();
      }
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
