#include <catch2/catch.hpp>

#include <share/io/scream_io_utils.hpp>
#include <share/util/scream_time_stamp.hpp>

#include <fstream>

TEST_CASE ("find_filename_in_rpointer") {
  using namespace scream;

  ekat::Comm comm(MPI_COMM_WORLD);

  util::TimeStamp t0({2023,9,7},{12,0,0});
  util::TimeStamp t1({2023,9,7},{13,0,0});

  // Create a dummy rpointer
  std::ofstream rpointer ("rpointer.atm");

  rpointer << "foo.r." + t0.to_string() + ".nc\n";
  rpointer << "bar2.rhist." + t0.to_string() + ".nc\n";
  rpointer << "bar.rhist." + t0.to_string() + ".nc\n";
  rpointer.close();

  // Now test find_filename_in_rpointer with different inputs

  REQUIRE_THROWS (find_filename_in_rpointer("baz",false,comm,t0)); // wrong prefix
  REQUIRE_THROWS (find_filename_in_rpointer("bar",false,comm,t1)); // wrong timestamp
  REQUIRE_THROWS (find_filename_in_rpointer("bar",true, comm,t0)); // bar is not model restart
  REQUIRE_THROWS (find_filename_in_rpointer("foo",false,comm,t0)); // foo is model restart

  REQUIRE (find_filename_in_rpointer("bar", false,comm,t0)==("bar.rhist."+t0.to_string()+".nc"));
  REQUIRE (find_filename_in_rpointer("bar2",false,comm,t0)==("bar2.rhist."+t0.to_string()+".nc"));
  REQUIRE (find_filename_in_rpointer("foo", true, comm,t0)==("foo.r."+t0.to_string()+".nc"));
}

TEST_CASE ("io_control") {
  using namespace scream;

  util::TimeStamp t0({2023,9,7},{12,0,0});

  IOControl control;
  control.frequency = 2;
  control.timestamp_of_last_write = t0;

  SECTION ("none") {
    control.frequency_units = "none";
    REQUIRE (not control.output_enabled());
    REQUIRE (not control.is_write_step(t0));
  }

  SECTION ("never") {
    control.frequency_units = "never";
    REQUIRE (not control.output_enabled());
    REQUIRE (not control.is_write_step(t0));
  }

  SECTION ("nsteps") {
    control.frequency_units = "nsteps";
    auto t1 = t0 + 1;
    auto t2 = t1 + 1;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nsecs") {
    control.frequency_units = "nsecs";
    auto t1 = t0 + 1;
    auto t2 = t1 + 1;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nmins") {
    control.frequency_units = "nmins";
    auto t1 = t0 + 60;
    auto t2 = t1 + 60;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nhours") {
    control.frequency_units = "nhours";
    auto t1 = t0 + 3600;
    auto t2 = t1 + 3600;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("ndays") {
    control.frequency_units = "ndays";
    auto t1 = t0 + 86400;
    auto t2 = t1 + 86400;
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
  }

  SECTION ("nmonths") {
    control.frequency_units = "nmonths";
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
    util::TimeStamp t1({2024,9,7},{12,0,0});
    util::TimeStamp t2({2025,9,7},{12,0,0});
    util::TimeStamp t3({2025,9,7},{13,0,0});
    REQUIRE (control.output_enabled());
    REQUIRE (not control.is_write_step(t1));
    REQUIRE (control.is_write_step(t2));
    REQUIRE (not control.is_write_step(t3));
  }
}
