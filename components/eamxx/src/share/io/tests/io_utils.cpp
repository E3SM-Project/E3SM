#include <catch2/catch.hpp>

#include <share/io/eamxx_io_utils.hpp>
#include <share/io/eamxx_io_control.hpp>
#include <share/util/eamxx_time_stamp.hpp>

#include <fstream>

TEST_CASE ("find_filename_in_rpointer") {
  using namespace scream;

  constexpr auto AVG = OutputAvgType::Average;
  constexpr auto INST = OutputAvgType::Instant;

  ekat::Comm comm(MPI_COMM_WORLD);

  util::TimeStamp t0({2023,9,7},{12,0,0});
  util::TimeStamp t1({2023,9,7},{13,0,0});

  // Create a dummy rpointer
  std::ofstream rpointer ("rpointer.atm");

  IOControl foo_c, bar_c, bar2_c;
  foo_c.frequency  = 3; foo_c.frequency_units  = "nsteps";
  bar_c.frequency  = 1; bar_c.frequency_units  = "ndays";
  bar2_c.frequency = 6; bar2_c.frequency_units = "nhours";

  std::string suffix = ".np" + std::to_string(comm.size()) + "." + t0.to_string() + ".nc";
  std::string foo_fname  = "foo.r.INSTANT.nsteps_x3"     + suffix;
  std::string bar_fname  = "bar.rhist.AVERAGE.ndays_x1"  + suffix;
  std::string bar2_fname = "bar.rhist.AVERAGE.nhours_x6" + suffix;

  rpointer << foo_fname<< "\n";
  rpointer << bar_fname<< "\n";
  rpointer << bar2_fname << "\n";
  rpointer.close();

  // Now test find_filename_in_rpointer with different inputs
  REQUIRE_THROWS (find_filename_in_rpointer("baz",false,comm,t0,false,AVG)); // missing control (needed for rhist files)
  REQUIRE_THROWS (find_filename_in_rpointer("baz",false,comm,t0,false,AVG,bar_c)); // wrong prefix
  REQUIRE_THROWS (find_filename_in_rpointer("bar",false,comm,t1,false,AVG,bar_c)); // wrong timestamp
  REQUIRE_THROWS (find_filename_in_rpointer("bar",true, comm,t0,false,AVG,bar_c)); // bar is not model restart
  REQUIRE_THROWS (find_filename_in_rpointer("bar",false,comm,t0,false,INST,bar_c)); // wrong avg type
  REQUIRE_THROWS (find_filename_in_rpointer("bar",false,comm,t0,false,INST,bar2_c)); // wrong freq specs
  REQUIRE_THROWS (find_filename_in_rpointer("foo",false,comm,t0,false,INST,foo_c)); // foo is model restart
  REQUIRE_THROWS (find_filename_in_rpointer("foo",true, comm,t0,false,AVG)); // model restart MUST be INSTANT
  auto not_found = find_filename_in_rpointer("bar",false,comm,t0,true,INST,bar2_c); // Allowed to not find it
  REQUIRE (not_found=="");

  REQUIRE (find_filename_in_rpointer("bar",false,comm,t0,false,AVG,bar_c)==bar_fname);
  REQUIRE (find_filename_in_rpointer("bar",false,comm,t0,false,AVG,bar2_c)==bar2_fname);
  REQUIRE (find_filename_in_rpointer("foo",true, comm,t0)==foo_fname);
}

TEST_CASE ("io_control") {
  using namespace scream;

  util::TimeStamp t0({2023,9,7},{12,0,0});

  IOControl control;
  control.frequency = 2;
  control.last_write_ts = t0;

  SECTION ("none") {
    control.frequency_units = "none";
    REQUIRE (not control.output_enabled());
  }

  SECTION ("never") {
    control.frequency_units = "never";
    REQUIRE (not control.output_enabled());
  }

  SECTION ("nsteps") {
    control.frequency_units = "nsteps";
    control.compute_next_write_ts();
    auto t1 = t0 + 1;
    auto t2 = t1 + 1;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nsecs") {
    control.frequency_units = "nsecs";
    control.compute_next_write_ts();
    auto t1 = t0 + 1;
    auto t2 = t1 + 1;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nmins") {
    control.frequency_units = "nmins";
    control.compute_next_write_ts();
    auto t1 = t0 + 60;
    auto t2 = t1 + 60;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nhours") {
    control.frequency_units = "nhours";
    control.compute_next_write_ts();
    auto t1 = t0 + 3600;
    auto t2 = t1 + 3600;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("ndays") {
    control.frequency_units = "ndays";
    control.compute_next_write_ts();
    auto t1 = t0 + 86400;
    auto t2 = t1 + 86400;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nmonths") {
    control.frequency_units = "nmonths";
    control.compute_next_write_ts();
    util::TimeStamp t1({2023,10,7},{12,0,0});
    util::TimeStamp t2({2023,11,7},{12,0,0});
    util::TimeStamp t3({2023,11,7},{13,0,0});
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
    REQUIRE (not control.is_write_step(t3));
  }

  SECTION ("nyears") {
    control.frequency_units = "nyears";
    control.compute_next_write_ts();
    util::TimeStamp t1({2024,9,7},{12,0,0});
    util::TimeStamp t2({2025,9,7},{12,0,0});
    util::TimeStamp t3({2025,9,7},{13,0,0});
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
    REQUIRE (not control.is_write_step(t3));
  }
}

TEST_CASE ("parse_cf_time_units") {
  using namespace scream;

  SECTION ("seconds since with full timestamp") {
    auto ts_units = "seconds since 1970-01-01 00:00:00";
    auto ts = parse_cf_time_units(ts_units);
    REQUIRE (ts.get_year() == 1970);
    REQUIRE (ts.get_month() == 1);
    REQUIRE (ts.get_day() == 1);
    REQUIRE (ts.get_hours() == 0);
    REQUIRE (ts.get_minutes() == 0);
    REQUIRE (ts.get_seconds() == 0);
  }

  SECTION ("days since with full timestamp") {
    auto ts_units = "days since 2000-01-01 00:00:00";
    auto ts = parse_cf_time_units(ts_units);
    REQUIRE (ts.get_year() == 2000);
    REQUIRE (ts.get_month() == 1);
    REQUIRE (ts.get_day() == 1);
    REQUIRE (ts.get_hours() == 0);
    REQUIRE (ts.get_minutes() == 0);
    REQUIRE (ts.get_seconds() == 0);
  }

  SECTION ("hours since with full timestamp") {
    auto ts_units = "hours since 2023-09-07 12:30:45";
    auto ts = parse_cf_time_units(ts_units);
    REQUIRE (ts.get_year() == 2023);
    REQUIRE (ts.get_month() == 9);
    REQUIRE (ts.get_day() == 7);
    REQUIRE (ts.get_hours() == 12);
    REQUIRE (ts.get_minutes() == 30);
    REQUIRE (ts.get_seconds() == 45);
  }

  SECTION ("minutes since with full timestamp") {
    auto ts_units = "minutes since 2015-06-15 08:15:30";
    auto ts = parse_cf_time_units(ts_units);
    REQUIRE (ts.get_year() == 2015);
    REQUIRE (ts.get_month() == 6);
    REQUIRE (ts.get_day() == 15);
    REQUIRE (ts.get_hours() == 8);
    REQUIRE (ts.get_minutes() == 15);
    REQUIRE (ts.get_seconds() == 30);
  }

  SECTION ("date only (no time)") {
    auto ts_units = "days since 2010-12-25";
    auto ts = parse_cf_time_units(ts_units);
    REQUIRE (ts.get_year() == 2010);
    REQUIRE (ts.get_month() == 12);
    REQUIRE (ts.get_day() == 25);
    REQUIRE (ts.get_hours() == 0);
    REQUIRE (ts.get_minutes() == 0);
    REQUIRE (ts.get_seconds() == 0);
  }

  SECTION ("with timezone (should be ignored)") {
    auto ts_units = "seconds since 1980-03-15 10:20:30 UTC";
    auto ts = parse_cf_time_units(ts_units);
    REQUIRE (ts.get_year() == 1980);
    REQUIRE (ts.get_month() == 3);
    REQUIRE (ts.get_day() == 15);
    REQUIRE (ts.get_hours() == 10);
    REQUIRE (ts.get_minutes() == 20);
    REQUIRE (ts.get_seconds() == 30);
  }

  SECTION ("with fractional seconds (should be truncated)") {
    auto ts_units = "seconds since 2020-05-10 14:25:36.789";
    auto ts = parse_cf_time_units(ts_units);
    REQUIRE (ts.get_year() == 2020);
    REQUIRE (ts.get_month() == 5);
    REQUIRE (ts.get_day() == 10);
    REQUIRE (ts.get_hours() == 14);
    REQUIRE (ts.get_minutes() == 25);
    REQUIRE (ts.get_seconds() == 36); // fractional part is ignored
  }

  SECTION ("invalid format (no 'since')") {
    auto ts_units = "seconds from 1970-01-01 00:00:00";
    REQUIRE_THROWS (parse_cf_time_units(ts_units));
  }
}
