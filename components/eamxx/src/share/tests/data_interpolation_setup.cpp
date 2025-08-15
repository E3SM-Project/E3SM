#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/io/eamxx_io_utils.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
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

  // We use raw scorpio calls without decomp, so ensure we're in serial case
  EKAT_REQUIRE_MSG (comm.size()==1,
      "Error! You should run the data_interpolation_setup test with ONE rank.\n");

  // Create grid
  std::shared_ptr<const AbstractGrid> grid = create_point_grid("pg",ngcols,nlevs,comm);

  // Create and setup two files, so we can test both YearlyPeriodic and LinearHistory
  std::vector<std::string> files = {
    "data_interpolation_0",
    "data_interpolation_1"
  };

  for (auto int_same_as_mid : {true, false}) {
    auto suffix = int_same_as_mid ? "_no_ilev.nc" : ".nc";
    for (const std::string& fname : files) {
      scorpio::register_file(fname+suffix,scorpio::Write);

      scorpio::define_dim (fname+suffix,"ncol",ngcols);
      scorpio::define_dim (fname+suffix,"lev",nlevs);
      scorpio::define_dim (fname+suffix,"dim2",ncmps);
      scorpio::define_time(fname+suffix,"days since " + t_ref.to_string());
      if (not int_same_as_mid) {
        scorpio::define_dim (fname+suffix,"ilev",nlevs+1);
      }

      std::string ilev = int_same_as_mid ? "lev" : "ilev";

      scorpio::define_var(fname+suffix,"s2d",  {"ncol"},             "real", true);
      scorpio::define_var(fname+suffix,"s2d",  {"ncol"},             "real", true);
      scorpio::define_var(fname+suffix,"v2d",  {"ncol","dim2"},      "real", true);
      scorpio::define_var(fname+suffix,"s3d_m",{"ncol","lev"},       "real", true);
      scorpio::define_var(fname+suffix,"v3d_m",{"ncol","dim2","lev"},"real", true);
      scorpio::define_var(fname+suffix,"s3d_i",{"ncol",ilev},        "real", true);
      scorpio::define_var(fname+suffix,"v3d_i",{"ncol","dim2",ilev}, "real", true);

      // We keep p1d, p2d, p3d, hyam, and hybm as NOT time-dep
      scorpio::define_var(fname+suffix,"p1d",  {"lev"},"real", false);
      scorpio::define_var(fname+suffix,"p2d",  {"ncol"},"real", false);
      scorpio::define_var(fname+suffix,"p3d",  {"ncol","lev"},"real", false);
      scorpio::define_var(fname+suffix,"hyam", {"lev"},"real", false);
      scorpio::define_var(fname+suffix,"hybm", {"lev"},"real", false);

      scorpio::enddef(fname+suffix);
    }

    // Fields and some helper fields (for later)
    // NOTE: if we save a pressure field, there is not distinction
    //       between interfaces and midpoints in the file
    // NOTE: do not pad, so that we can grab pointers and pass them to scorpio
    auto base_fields = create_fields(grid,true, int_same_as_mid,false);
    auto fields      = create_fields(grid,false,int_same_as_mid,false);
    auto ones        = create_fields(grid,false,int_same_as_mid,false);
    for (auto& f : ones) {
      f.deep_copy(1);
    }
    // Loop over time, and add 30 to the value for the first 6 months,
    // and subtract 30 for the last 6 months. This guarantees that the data
    // is indeed periodic. We'll write at the 15th of each month
    // Generate three files:
    //   - one to be used for yearly-periodic interp
    //   - two to be used for linear-hystory interp
    util::TimeStamp time = get_first_slice_time ();

    // We keep pressures fields NOT time-dep, so we write outside the loop. Also write hyam/hybm here
    auto p1d = base_fields[6];
    auto p3d = base_fields[2].alias("p3d");
    auto p2d = base_fields[0].alias("p2d");
    auto hybm = p1d.alias("hybm");
    auto hyam = hybm.clone("hyam");
    p1d.sync_to_host();
    p2d.sync_to_host();
    p3d.sync_to_host();
    hyam.deep_copy(0); hyam.sync_to_host();
    for (const std::string& fname : files) {
      scorpio::write_var(fname+suffix,p1d.name(),p1d.get_internal_view_data<Real,Host>());
      scorpio::write_var(fname+suffix,p2d.name(),p2d.get_internal_view_data<Real,Host>());
      scorpio::write_var(fname+suffix,p3d.name(),p3d.get_internal_view_data<Real,Host>());
      scorpio::write_var(fname+suffix,hybm.name(),hybm.get_internal_view_data<Real,Host>());
      scorpio::write_var(fname+suffix,hyam.name(),hyam.get_internal_view_data<Real,Host>());
    }

    int nfields = fields.size() - 1; // Don't handle p1d, since it's done above
    for (int mm=0; mm<12; ++mm) {
      std::string file_name = "data_interpolation_" + std::to_string(mm/6) + suffix;

      // We start the files with July
      int mm_index = mm+6;
      scorpio::update_time(file_name,time.days_from(t_ref));
      for (int i=0; i<nfields; ++i) {
        auto& f = fields[i];
        f.deep_copy(base_fields[i]);
        f.update(ones[i],delta_data[ mm_index % 12],1.0);
        f.sync_to_host();
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
