#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/io/scream_io_utils.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"

namespace scream {

TEST_CASE ("data_interpolation_setup")
{
  // NOTE: ensure these match what's used in data_interpolation_tests.cpp
  constexpr int ngcols = data_ngcols;
  constexpr int nlevs  = data_nlevs;

  auto t_ref = get_t_ref();

  // Init test session
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  EKAT_REQUIRE_MSG (comm.size()==1,
      "Error! You should run the data_interpolation_setup test with ONE rank.\n");

  // Create grid
  std::shared_ptr<const AbstractGrid> grid = create_point_grid("pg",ngcols,nlevs,comm);

  // Create and setup two files, so we can test both YearlyPeriodic and LinearHistory
  std::vector<std::string> files = {
    "data_interpolation_0",
    "data_interpolation_1"
  };

  for (auto with_pressure : {true, false}) {
    auto suffix = with_pressure ? "_with_p.nc" : ".nc";
    for (const std::string& fname : files) {
      scorpio::register_file(fname+suffix,scorpio::Write);

      scorpio::define_dim (fname+suffix,"ncol",ngcols);
      scorpio::define_dim (fname+suffix,"lev",nlevs);
      if (not with_pressure) {
        scorpio::define_dim (fname+suffix,"ilev",nlevs+1);
      }
      scorpio::define_dim (fname+suffix,"dim2",ncmps);
      scorpio::define_time(fname+suffix,"days since " + t_ref.to_string());

      std::string ilev_tag = with_pressure ? "lev" : "ilev";

      scorpio::define_var(fname+suffix,"s2d",  {"ncol"},                 "real", true);
      scorpio::define_var(fname+suffix,"s2d",  {"ncol"},                 "real", true);
      scorpio::define_var(fname+suffix,"v2d",  {"ncol","dim2"},          "real", true);
      scorpio::define_var(fname+suffix,"s3d_m",{"ncol","lev"},           "real", true);
      scorpio::define_var(fname+suffix,"v3d_m",{"ncol","dim2","lev"},    "real", true);
      scorpio::define_var(fname+suffix,"s3d_i",{"ncol",ilev_tag},        "real", true);
      scorpio::define_var(fname+suffix,"v3d_i",{"ncol","dim2",ilev_tag}, "real", true);

      if (with_pressure) {
        scorpio::define_var(fname+suffix,"p1d",{"lev"},"real", false);
        scorpio::define_var(fname+suffix,"p3d",{"ncol","lev"},"real", true);
      }

      scorpio::enddef(fname+suffix);
    }

    // Fields and some helper fields (for later)
    // NOTE: if we save a pressure field, there is not distinction
    //       between interfaces and midpoints in the file
    auto base_fields = create_fields (grid,true,with_pressure);
    auto fields      = create_fields(grid,false,with_pressure);
    auto ones        = create_fields(grid,false,with_pressure);
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

    if (with_pressure) {
      // Create p1d as slice of p3d, and ensure it's the same on all ranks, then write it.
      auto p1d = fields.back().subfield(0,0).clone("p1d");
      auto comm = grid->get_comm();
      comm.broadcast(p1d.get_internal_view_data<Real,Host>(),nlevs,0);
      p1d.sync_to_dev();
      for (const std::string& fname : files) {
        scorpio::write_var(fname+suffix,p1d.name(),p1d.get_internal_view_data<Real,Host>());
      }
    }

    for (int mm=0; mm<12; ++mm) {
      std::string file_name = "data_interpolation_" + std::to_string(mm/6) + suffix;

      // We start the files with July
      int mm_index = mm+6;
      scorpio::update_time(file_name,time.days_from(t_ref));
      for (int i=0; i<nfields; ++i) {
        auto& f = fields[i];
        f.deep_copy(base_fields[i]);
        f.update(ones[i],delta_data[ mm_index % 12],1.0);
        scorpio::write_var(file_name,f.name(),f.get_internal_view_data<Real,Host>());
      }
      time += 86400*time.days_in_curr_month();
    }

    for (const std::string& fname : files) {
      write_timestamp(fname+suffix,"reference_time_stamp",t_ref);
      scorpio::release_file(fname+suffix);
    }
  }

  // Now write a map file for horiz remap, that splits each dof interval in two
  const int ngdofs_src = data_ngcols;
  const int ngdofs_tgt = fine_ngcols;

  // Existing dofs are "copied", added dofs are averaged from neighbors
  const int nnz = ngdofs_src + 2*(ngdofs_src-1);
  
  std::string filename = map_file_name;
  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::define_dim(filename, "n_a", ngdofs_src);
  scorpio::define_dim(filename, "n_b", ngdofs_tgt);
  scorpio::define_dim(filename, "n_s", nnz);

  scorpio::define_var(filename, "col", {"n_s"}, "int");
  scorpio::define_var(filename, "row", {"n_s"}, "int");
  scorpio::define_var(filename, "S",   {"n_s"}, "double");

  scorpio::enddef(filename);

  std::vector<int> col(nnz), row(nnz);
  std::vector<double> S(nnz);
  for (int i=0,nnz=0; i<ngdofs_tgt; ++i) {
    if (i % 2 == 0) {
      // Fine grid point also in the src grid: only diag entry
      col[nnz] = i / 2;
      row[nnz] = i;
        S[nnz] = 1.0;

      ++nnz;
    } else {
      // Add 0.5 times the left src dof and 0.5 times the right src dof
      col[nnz] = i / 2;
      row[nnz] = i;
        S[nnz] = 0.5;
      ++nnz;

      col[nnz] = (i / 2) + 1;
      row[nnz] = i;
        S[nnz] = 0.5;
      ++nnz;
    }
  }

  scorpio::write_var(filename,"row",row.data());
  scorpio::write_var(filename,"col",col.data());
  scorpio::write_var(filename,"S",  S.data());

  scorpio::release_file(filename);

  scorpio::finalize_subsystem();
}

} // anonymous namespace
