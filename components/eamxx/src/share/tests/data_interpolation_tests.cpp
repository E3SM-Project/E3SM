#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/eamxx_data_interpolation.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/scream_config.hpp"

// NOTE: ensure these are the same used in data_interpolation_setup.cpp
constexpr int data_ncols = 12;
constexpr int data_nlevs = 32;

namespace scream {

using strvec_t = std::vector<std::string>;

std::shared_ptr<DataInterpolation>
create_interp (const std::shared_ptr<const AbstractGrid>& grid,
               const std::vector<Field>& fields)
{
  return std::make_shared<DataInterpolation>(grid,fields);
}

TEST_CASE ("exceptions")
{
  // Test correctness of some exception handling inside the DataInterpolation source code
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);
  auto grid = create_point_grid("pg",data_ncols,data_nlevs,comm);

  auto fields = create_fields(grid,false);

  REQUIRE_THROWS (create_interp(nullptr,fields)); // Invalid grid pointer

  auto interp = create_interp(grid,fields);

  strvec_t files = {"/etc/shadow"};
  REQUIRE_THROWS (interp->setup_time_database(files,util::TimeLine::Linear)); // Input file not readable
  
  interp->setup_time_database({"./data_interpolation_0.nc"},util::TimeLine::Linear);
  util::TimeStamp t0 ({2000,1,1},{0,0,0});
  REQUIRE_THROWS (interp->init_data_interval(t0)); // linear timeline, but t0<first_slice

  util::TimeStamp t1 ({2020,1,1},{0,0,0});
  REQUIRE_THROWS (interp->init_data_interval(t1)); // linear timeline, but t0>last_slice

  scorpio::finalize_subsystem();
}

TEST_CASE ("interpolation")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  // Regardless of how EAMxx is configured, ignore leap years for this test
  set_use_leap_year(false);

  scorpio::init_subsystem(comm);

  Real tol = std::numeric_limits<Real>::epsilon()*10;
  SECTION ("only-time")
  {
    auto grid = create_point_grid("pg",data_ncols,data_nlevs,comm);

    // Create a bunch of copies of the fields. All but the first set of them
    // are used in the test phase to check the result against what's expected
    auto fields   = create_fields(grid,false);
    auto base_f   = create_fields(grid,true);
    auto ones     = create_fields (grid,false);
    auto diff     = create_fields (grid,false);
    auto expected = create_fields (grid,false);
    for (auto& f : ones) {
      f.deep_copy(1);
    }
    int nfields = fields.size();

    auto interp = create_interp(grid,fields);

    SECTION ("periodic")
    {
      strvec_t files = {"data_interpolation_0.nc","data_interpolation_1.nc"};

      util::TimeStamp t0 ({2020,1,1},{0,0,0});

      // t_beg/t_end keep track of the interval [beg,end] where we currently are,
      // where beg/end are timestamps at which we have data.
      // We assume we start with t0 after the last input slice of the 1st year.
      // NOTE: -365*spd is since we actually need to *rewind* t0, if we want beg<t0<end
      util::TimeStamp t_beg = t0 + get_last_slice_time().frac_of_year_in_days()*spd - 365*spd;
      util::TimeStamp t_end = t_beg + t_beg.days_in_curr_month()*spd;

      interp->setup_time_database(files,util::TimeLine::YearlyPeriodic);
      interp->init_data_interval(t0);

      // Loop for two year at a 20 day increment
      int dt = 20*spd;
      for (auto time = t0+dt; time.days_from(t0)<365; time+=dt) {
        if (t_end<time) {
          // update t_beg/t_end
          t_beg = t_end;
          t_end += t_end.days_in_curr_month()*spd;
        }

        // Compute the delta to add to field base value for mm_beg and mm_end
        int mm_beg = t_beg.get_month()-1;
        int mm_end = t_end.get_month()-1;

        // Since input data is at the 15th of the month, we need to compute
        // the distance between current time and the beg slice, then use it
        // to do a convex interpolation between f(t=t_beg) and f(t=t_end).
        util::TimeInterval time_from_beg(t_beg,time,util::TimeLine::YearlyPeriodic);;
        double alpha = time_from_beg.length / t_beg.days_in_curr_month();
        double delta = delta_data[mm_beg]*(1-alpha) + delta_data[mm_end]*alpha;

        if (alpha<0 or alpha>1) {
          std::cout << "TEST ERROR:\n"
                    << " t beg: " << t_beg.to_string() << "\n"
                    << " t end: " << t_end.to_string() << "\n"
                    << " time : " << time.to_string() << "\n"
                    << " t-beg: " << time_from_beg.length << "\n"
                    << " days in mm_beg: " << t_beg.days_in_curr_month() << "\n"
                    << " alpha: " << alpha << "\n"
                    << " delta_beg: " << delta_data[mm_beg] << "\n"
                    << " delta_end: " << delta_data[mm_end] << "\n"
                    << " delta: " << delta << "\n";
        }
        // Compute expected difference from base value
        interp->run(time);  
        for (int i=0; i<nfields; ++i) {
          // Compute expected, then subtract computed field
          expected[i].update(base_f[i],1,0);
          expected[i].update(ones[i],delta,1.0);
          diff[i].update(fields[i],1,1);
          diff[i].update(expected[i],-1,1);
          diff[i].scale(1.0 / frobenius_norm<Real>(expected[i]));
          if (frobenius_norm<Real>(diff[i])>=tol) {
            auto n = fields[i].name();
            print_field_hyperslab(fields[i].alias(n+"_computed"));
            print_field_hyperslab(expected[i].alias(n+"_expected"));
            print_field_hyperslab(diff[i].alias(n+"_diff"));
          }
          REQUIRE (frobenius_norm<Real>(diff[i])<tol);
        }
      }
    }

    SECTION ("linear")
    {
      // t_beg/t_end keep track of the interval [beg,end] where we currently are,
      // where beg/end are timestamps at which we have data.
      // We assume we start with t0 sometime between the first and second input slice.
      auto t_beg = get_first_slice_time();
      for (int mm=0; mm<6; ++mm) {
        t_beg += spd*t_beg.days_in_curr_month();
      }
      auto t_end = t_beg + spd*t_beg.days_in_curr_month();
      auto t0    = t_beg + spd*10;

      strvec_t files = {"data_interpolation_1.nc","data_interpolation_2.nc"};
      interp->setup_time_database(files,util::TimeLine::Linear);
      interp->init_data_interval(t0);

      // Loop for two year at a 20 day increment
      int dt = 20*spd;
      for (auto time = t0+dt; time.days_from(t0)<200; time+=dt) {
        if (t_end<time) {
          // update t_beg/t_end
          t_beg = t_end;
          t_end += t_end.days_in_curr_month()*spd;
        }

        // Compute the delta to add to field base value for mm_beg and mm_end
        int mm_beg = t_beg.get_month()-1;
        int mm_end = t_end.get_month()-1;

        // Since input data is at the 15th of the month, we need to compute
        // the distance between current time and the beg slice, then use it
        // to do a convex interpolation between f(t=t_beg) and f(t=t_end).
        util::TimeInterval time_from_beg(t_beg,time,util::TimeLine::Linear);;
        double alpha = time_from_beg.length / t_beg.days_in_curr_month();
        double delta = delta_data[mm_beg]*(1-alpha) + delta_data[mm_end]*alpha;

        if (alpha<0 or alpha>1) {
          std::cout << "TEST ERROR:\n"
                    << " t beg: " << t_beg.to_string() << "\n"
                    << " t end: " << t_end.to_string() << "\n"
                    << " time : " << time.to_string() << "\n"
                    << " t-beg: " << time_from_beg.length << "\n"
                    << " days in mm_beg: " << t_beg.days_in_curr_month() << "\n"
                    << " alpha: " << alpha << "\n"
                    << " delta_beg: " << delta_data[mm_beg] << "\n"
                    << " delta_end: " << delta_data[mm_end] << "\n"
                    << " delta: " << delta << "\n";
        }
        // Compute expected difference from base value
        interp->run(time);  
        for (int i=0; i<nfields; ++i) {
          // Compute expected, then subtract computed field
          expected[i].update(base_f[i],1,0);
          expected[i].update(ones[i],delta,1.0);
          diff[i].update(fields[i],1,1);
          diff[i].update(expected[i],-1,1);
          diff[i].scale(1.0 / frobenius_norm<Real>(expected[i]));
          if (frobenius_norm<Real>(diff[i])>=tol) {
            print_field_hyperslab(fields[i]);
            print_field_hyperslab(expected[i]);
            print_field_hyperslab(diff[i]);
          }
          REQUIRE (frobenius_norm<Real>(diff[i])<tol);
        }
      }
    }
  }

  scorpio::finalize_subsystem();
}

} // anonymous namespace
