#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_config.hpp"

namespace scream {

// Test that setup_yearly_periodic_time_database with a start_ts correctly rotates
// the slice array, so that the logical "first interval" starts at start_ts.
//
// Data: three files, where the first two cover 12 months starting from July
// (slices 0..11 = Jul..Jun), while the third adds 6 extra months with a large
// offset in the values. If start_ts truncation is broken, the wrap-around test
// below will pick those offset slices and fail.
// We set start_ts = get_first_slice_time() + 3 months (= October 15th).
// After rotation, slices[0]=Oct, slices[1]=Nov, ..., slices[8]=Jun, slices[9]=Jul, ...
//
// We then verify that running at a time between Oct and Nov produces the expected
// linear interpolation between delta_data[9]=90 (Oct) and delta_data[10]=60 (Nov).

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
    "data_interpolation_2.nc"
  };

  // Build start_ts = October 15th (index 3 in the original slice array)
  // get_first_slice_time() = July 15th
  auto start_ts = get_first_slice_time();
  for (int mm=0; mm<3; ++mm) {
    start_ts += spd * start_ts.days_in_curr_month();
  }
  // start_ts is now October 15th (same year as first slice)

  // t0 = October 30th: 15 days after start_ts, within [Oct15, Nov15]
  auto t0 = start_ts + 15*spd;

  auto fields = create_fields(data_grid, false);
  fields.pop_back(); // Don't interpolate p1d

  auto interp = create_interp(data_grid, fields);
  interp->setup_periodic_time_database(files,
      DataInterpolation::Linear,
      start_ts);     // start_ts: rotate so Oct is first slice
  interp->create_horiz_remappers("");  // no hremap (data_grid == src_grid)
  interp->create_vert_remapper();      // no vremap (VRemapType::None)
  interp->init_data_interval(t0);

  auto base     = create_fields(data_grid, true);
  auto ones     = create_fields(data_grid, false);
  auto diff     = create_fields(data_grid, false);
  auto expected = create_fields(data_grid, false);
  for (auto& f : ones) { f.deep_copy(1); }

  int nfields = static_cast<int>(fields.size());

  // Helper: verify fields at time ts with expected interpolation coeff alpha
  auto check = [&](const util::TimeStamp& ts, double alpha) {
    // delta_data indices for Oct=9, Nov=10
    double delta = delta_data[9]*(1-alpha) + delta_data[10]*alpha;

    interp->run(ts);

    for (int i=0; i<nfields; ++i) {
      expected[i].update(base[i], 1, 0);
      expected[i].update(ones[i], delta, 1.0);

      diff[i].update(fields[i],   1, 1);
      diff[i].update(expected[i], -1, 1);
      diff[i].scale(1.0 / frobenius_norm(expected[i]).as<Real>());
      REQUIRE(frobenius_norm(diff[i]).as<Real>() < tol);
    }
  };

  // t0 is 15 days into October's 31-day interval
  auto timeline = util::TimeLine::YearlyPeriodic;
  util::TimeInterval from_oct(start_ts, t0, timeline);
  double alpha_t0 = from_oct.length / start_ts.days_in_curr_month();

  SECTION ("within-oct-nov-interval") {
    // alpha = 15/31 ≈ 0.484
    check(t0, alpha_t0);
  }

  SECTION ("near-end-of-oct-nov-interval") {
    // 28 days into October (just before Nov 15)
    auto ts = start_ts + 28*spd;
    util::TimeInterval from_oct2(start_ts, ts, timeline);
    double alpha = from_oct2.length / start_ts.days_in_curr_month();
    check(ts, alpha);
  }

  SECTION ("wrap-around-sep-to-oct") {
    // Test that the wrap-around at the end of the rotated array works:
    // the last interval in the rotated array is Sep→Oct.
    // September 15th (slices[11] in the rotated array) to October 15th (slices[0]).
    // t_sep = start_ts + 11 months
    auto t_sep = start_ts;
    for (int mm=0; mm<11; ++mm) {
      t_sep += spd * t_sep.days_in_curr_month();
    }
    // t_sep is now Sep 15th.  We test a point 10 days after Sep 15th,
    // in the interval [Sep15, Oct15].  Expected: interp between
    // delta_data[8]=120 (Sep) and delta_data[9]=90 (Oct).
    auto ts = t_sep + 10*spd;
    util::TimeInterval from_sep(t_sep, ts, timeline);
    double alpha = from_sep.length / t_sep.days_in_curr_month();
    double delta = delta_data[8]*(1-alpha) + delta_data[9]*alpha;

    interp->run(ts);

    for (int i=0; i<nfields; ++i) {
      expected[i].update(base[i], 1, 0);
      expected[i].update(ones[i], delta, 1.0);

      diff[i].update(fields[i],   1, 1);
      diff[i].update(expected[i], -1, 1);
      diff[i].scale(1.0 / frobenius_norm(expected[i]).as<Real>());
      REQUIRE(frobenius_norm(diff[i]).as<Real>() < tol);
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
