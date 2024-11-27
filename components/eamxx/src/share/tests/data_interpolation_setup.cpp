#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/io/scream_io_utils.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"

namespace scream {

TEST_CASE ("data_interpolation_setup")
{
  // NOTE: ensure these match what's used in data_interpolation_tests.cpp
  constexpr int ncols = 12;
  constexpr int nlevs = 32;

  auto t_ref = get_t_ref();

  // Init test session
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  EKAT_REQUIRE_MSG (comm.size()==1,
      "Error! You should run the data_interpolation_setup test with ONE rank.\n");

  // Create grid
  std::shared_ptr<const AbstractGrid> grid = create_point_grid("pg",ncols,nlevs,comm);

  // Create and setup two files, so we can test both YearlyPeriodic and LinearHistory
  std::vector<std::string> files = {
    "data_interpolation_0.nc",
    "data_interpolation_1.nc",
    "data_interpolation_2.nc",
    "data_interpolation_3.nc"
  };
  for (const std::string& fname : files) {
    scorpio::register_file(fname,scorpio::Write);
    scorpio::define_dim(fname,"ncol",ncols);
    scorpio::define_dim(fname,"lev",nlevs);
    scorpio::define_dim(fname,"ilev",nlevs+1);
    scorpio::define_dim(fname,"dim2",ncmps);
    scorpio::define_time(fname,"days since " + t_ref.to_string());

    scorpio::define_var(fname,"s2d",  {"ncol"},              "real", true);
    scorpio::define_var(fname,"v2d",  {"ncol","dim2"},       "real", true);
    scorpio::define_var(fname,"s3d_m",{"ncol","lev"},        "real", true);
    scorpio::define_var(fname,"v3d_m",{"ncol","dim2","lev"}, "real", true);
    scorpio::define_var(fname,"s3d_i",{"ncol","ilev"},       "real", true);
    scorpio::define_var(fname,"v3d_i",{"ncol","dim2","ilev"},"real", true);

    scorpio::enddef(fname);
  }

  // Fields and some helper fields (for later)
  auto base_fields = create_fields (grid,true);
  auto fields = create_fields(grid,false);
  auto ones = create_fields(grid,false);
  for (const auto& f : ones) {
    f.deep_copy(1);
  }
  int nfields = fields.size();

  // Loop over time, and add 30 to the value for the first 6 months,
  // and subtract 30 for the last 6 months. This guarantees that the data
  // is indeed periodic. We'll write at the 15th of each month
  // Generate three files:
  //   - one to be used for yearly-periodic interp
  //   - two to be used for linear-hystory interp
  util::TimeStamp time = get_first_slice_time ();
  for (int mm=0; mm<24; ++mm) {
    std::string file_name = "data_interpolation_" + std::to_string(mm/6) + ".nc";

    scorpio::update_time(file_name,time.days_from(t_ref));
    for (int i=0; i<nfields; ++i) {
      auto& f = fields[i];
      f.deep_copy(base_fields[i]);
      f.update(ones[i],delta_data[mm % 12],1.0);
      scorpio::write_var(file_name,f.name(),f.get_internal_view_data<Real,Host>());
    }
    time += 86400*time.days_in_curr_month();
  }

  for (const std::string& fname : files) {
    write_timestamp(fname,"reference_time_stamp",t_ref);
    scorpio::release_file(fname);
  }

  scorpio::finalize_subsystem();
}

} // anonymous namespace
