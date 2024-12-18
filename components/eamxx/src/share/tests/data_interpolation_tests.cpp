#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/eamxx_data_interpolation.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/scream_config.hpp"

namespace scream {

constexpr auto tol = std::numeric_limits<Real>::epsilon()*10;
constexpr auto P1D = DataInterpolation::Static1D;
constexpr auto P3D = DataInterpolation::Dynamic3D;
using strvec_t = std::vector<std::string>;

util::TimeStamp reset_year (const util::TimeStamp& t_in, int yy)
{
  auto date = t_in.get_date();
  auto time = t_in.get_time();
  date[0] = yy;
  return util::TimeStamp(date,time);
}

std::shared_ptr<DataInterpolation>
create_interp (const std::shared_ptr<const AbstractGrid>& grid,
               const std::vector<Field>& fields)
{
  return std::make_shared<DataInterpolation>(grid,fields);
}

void root_print (const ekat::Comm& comm,
                 const std::string& msg)
{
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

// Run the data interpolation to the input grid, and check against expected values
void run_tests (const std::shared_ptr<const AbstractGrid>& grid,
                const strvec_t& input_files, util::TimeStamp t_beg,
                const util::TimeLine timeline,
                const DataInterpolation::VRemapType vr_type = DataInterpolation::None)
{
  auto t_end = t_beg + t_beg.days_in_curr_month()*spd;
  auto t0 = t_beg + (t_end-t_beg)/2;

  // These are the fields we will compute
  auto fields = create_fields(grid,false,false);

  std::string map_file = grid->get_num_global_dofs()==data_ngcols ? "" : map_file_name;

  // These are used to check the answer
  auto base_f   = create_fields(grid,true);
  auto ones     = create_fields(grid,false);
  auto diff     = create_fields(grid,false);
  auto expected = create_fields(grid,false);
  for (auto& f : ones) {
    f.deep_copy(1);
  }
  int nfields = fields.size();

  std::string data_pname = vr_type==P1D ? "p1d" : "p3d";  // if vr_type==None, it's not used anyways
  auto model_pmid = base_f[2].clone("pmid"); // ensure the 2nd field is s3d_m
  auto model_pint = base_f[4].clone("pint"); // ensure the 4th field is s3d_i

  auto interp = create_interp(grid,fields);
  interp->setup_time_database(input_files,util::TimeLine::YearlyPeriodic);
  interp->setup_remappers (map_file,vr_type,data_pname,model_pmid,model_pint);
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
    util::TimeInterval time_from_beg(t_beg,time,timeline);
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

TEST_CASE ("exceptions")
{
  // Test correctness of some exception handling inside the DataInterpolation source code
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);
  auto grid = create_point_grid("pg",data_ngcols,data_nlevs,comm);

  auto fields = create_fields(grid,false,false);

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

  auto data_grid   = create_point_grid("pg",data_ngcols,data_nlevs,comm);
  auto hfine_grid  = create_point_grid("pg",fine_ngcols,data_nlevs,comm);
  auto vfine_grid  = create_point_grid("pg",data_ngcols,fine_nlevs,comm);
  auto hvfine_grid = create_point_grid("pg",fine_ngcols,fine_nlevs,comm);

  SECTION ("periodic") {
    // We assume we start with t0 sometime between the first and second input slice.
    auto t_beg = reset_year(get_last_slice_time(),2019);
    auto timeline = util::TimeLine::YearlyPeriodic;
    strvec_t files = {"data_interpolation_0.nc","data_interpolation_1.nc"};
    strvec_t files_with_p = {"data_interpolation_0_with_p.nc","data_interpolation_1_with_p.nc"};

    SECTION ("no-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=NO ..........\n");
        run_tests (data_grid,files,t_beg,timeline);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p1d .........\n");
        run_tests (data_grid,files_with_p,t_beg,timeline,P1D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p3d .........\n");
        run_tests (data_grid,files_with_p,t_beg,timeline,P3D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p3d ......... PASS\n");
      }
    }

    SECTION ("yes-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=NO ..........\n");
        run_tests (hfine_grid,files,t_beg,timeline);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p1d .........\n");
        run_tests (data_grid,files_with_p,t_beg,timeline,P1D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p3d .........\n");
        run_tests (data_grid,files_with_p,t_beg,timeline,P3D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p3d ......... PASS\n");
      }
    }
  }

  SECTION ("linear") {
    // We assume we start with t0 sometime between the first and second input slice.
    auto t_beg = get_first_slice_time();
    auto timeline = util::TimeLine::Linear;
    strvec_t files = {"data_interpolation_0.nc","data_interpolation_1.nc"};

    SECTION ("no-horiz-no-vert") {
      root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=NO ..........\n");
      run_tests (data_grid,files,t_beg,timeline);
      root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=NO .......... PASS\n");
    }

    SECTION ("yes-horiz-no-vert") {
      root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=NO ..........\n");
      run_tests (hfine_grid,files,t_beg,timeline);
      root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=NO .......... PASS\n");
    }
  }

  scorpio::finalize_subsystem();
}

} // anonymous namespace
