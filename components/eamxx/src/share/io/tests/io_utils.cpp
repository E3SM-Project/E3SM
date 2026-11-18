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
  // File naming: eamxx.{r|rh}[.{stream_name}].{output_spec}.{date_spec}.nc (no case name for unit tests)
  std::string foo_fname  = "eamxx.r.INSTANT.nsteps_x3"    + suffix;
  std::string bar_fname  = "eamxx.rh.bar.AVERAGE.ndays_x1"  + suffix;
  std::string bar2_fname = "eamxx.rh.bar.AVERAGE.nhours_x6" + suffix;

  std::string dot_fname  = "eamxx.rh.my.stream-v1.AVERAGE.ndays_x1" + suffix;
  rpointer << foo_fname<< "\n";
  rpointer << dot_fname<< "\n";
  rpointer << bar_fname<< "\n";
  rpointer << bar2_fname << "\n";
  rpointer.close();

  // Now test find_filename_in_rpointer with different inputs
  // Note: prefix is [<case>.]eamxx; stream_name is only needed for history restart files
  auto f = [&](const std::string& stream, bool model_rest, const util::TimeStamp& t,
               bool allow_nf, OutputAvgType avg, const IOControl& c) {
    return find_filename_in_rpointer("eamxx",model_rest,comm,t,allow_nf,avg,c,stream);
  };
  REQUIRE_THROWS (find_filename_in_rpointer("eamxx",false,comm,t0,false,AVG,IOControl(),"bar")); // missing control (needed for rhist files)
  REQUIRE_THROWS (find_filename_in_rpointer("eamxx",false,comm,t0,false,AVG,bar_c)); // missing stream name (needed for rhist files)
  REQUIRE_THROWS (f("baz",false,t0,false,AVG,bar_c));  // wrong stream
  REQUIRE_THROWS (f("ba",false,t0,false,AVG,bar_c));   // stream is only a prefix of the real one
  REQUIRE_THROWS (f("bar",false,t1,false,AVG,bar_c));  // wrong timestamp
  REQUIRE_THROWS (f("bar",true, t0,false,AVG,bar_c));  // bar is not model restart
  REQUIRE_THROWS (f("bar",false,t0,false,INST,bar_c)); // wrong avg type
  REQUIRE_THROWS (f("bar",false,t0,false,INST,bar2_c));// wrong freq specs
  REQUIRE_THROWS (f("bar",false,t0,false,INST,foo_c)); // foo is model restart
  REQUIRE_THROWS (f("",   true, t0,false,AVG,IOControl())); // model restart MUST be INSTANT
  REQUIRE (f("bar",false,t0,true,INST,bar2_c)==""); // Allowed to not find it

  REQUIRE (f("bar",false,t0,false,AVG,bar_c)==bar_fname);
  REQUIRE (f("bar",false,t0,false,AVG,bar2_c)==bar2_fname);
  REQUIRE (f("my.stream-v1",false,t0,false,AVG,bar_c)==dot_fname); // dots/dashes in stream name
  REQUIRE (f("my",false,t0,true,AVG,bar_c)==""); // 'my' must not match 'my.stream-v1'
  REQUIRE (find_filename_in_rpointer("eamxx",true,comm,t0)==foo_fname);
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
