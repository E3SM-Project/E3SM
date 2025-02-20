#include <catch2/catch.hpp>

#include "data_interpolation_tests.hpp"

#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_data_interpolation.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/eamxx_config.hpp"

namespace scream {

// Give ourselves some room for roundoff errors, since our manual
// evaluation may be different (in finite prec) from the one in the class.
constexpr auto tol = std::numeric_limits<Real>::epsilon()*10;

// NOTE: due to how we build p3d, it is true that p3d is the tensor product of
//       p1d and p2d. Since hybm=p1d and hyam=0, the P2D case is the same as
//       the P3D one, once the 3d pressure has been reconstructed
constexpr auto NOP = DataInterpolation::None;
constexpr auto P1D = DataInterpolation::Static1D;
constexpr auto P2D = DataInterpolation::Dynamic3DRef;
constexpr auto P3D = DataInterpolation::Dynamic3D;

using strvec_t = std::vector<std::string>;
using namespace ShortFieldTagsNames;

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

  auto ncols = grid->get_num_local_dofs();
  auto nlevs = grid->get_num_vertical_levels();

  auto vcoarse_grid = grid->clone("vcoarse",true);
  vcoarse_grid->reset_num_vertical_lev(data_nlevs);

  // These are the fields we will compute
  auto fields = create_fields(grid,false,false);
  fields.pop_back(); // We don't interpolate p1d...

  std::string map_file = grid->get_num_global_dofs()==data_ngcols ? "" : map_file_name;

  // These are used to check the answer
  auto base     = create_fields(grid,true);
  auto ones     = create_fields(grid,false);
  auto diff     = create_fields(grid,false);
  auto expected = create_fields(grid,false);
  for (auto& f : ones) {
    f.deep_copy(1);
  }

  auto model_pmid = base[2].clone("pmid");
  auto model_pint = base[4].clone("pint");
  if (vr_type==P1D) {
    // It's complicated to test the static profile, since we'd have to really do
    // a manual interpolation. But setting all model pressure equal to the 1st col
    // of the pressure makes things doable
    auto comm = grid->get_comm();
    for (const Field& p : {model_pmid, model_pint}) {
      auto col_0 = p.subfield(0,0);
      auto len = col_0.get_header().get_identifier().get_layout().size();
      comm.broadcast(col_0.get_internal_view_data<Real,Host>(),len,0);
      col_0.sync_to_dev();

      for (int icol=0; icol<ncols; ++icol) {
        auto col_i = p.subfield(0,icol);
        col_i.deep_copy(col_0);
      }
    }
  }

  std::vector<Field> expected_vcoarse, base_vcoarse, ones_vcoarse;
  if (vr_type!=NOP) {
    // If we do remap, there is some P0 extrapolation,
    // for which we need to know the data at the top/bot
    // NOTE: the data has NO ilev coord in this case
    expected_vcoarse = create_fields(vcoarse_grid,false,true);
    base_vcoarse = create_fields(vcoarse_grid,true,true);
    ones_vcoarse = create_fields(vcoarse_grid,false,true);
    for (auto& f : ones_vcoarse) {
      f.deep_copy(1);
    }
  }

  DataInterpolation::RemapData remap_data;
  remap_data.hremap_file = map_file;
  remap_data.vr_type = vr_type;
  remap_data.pname =  // if vr_type==None, it's not used anyways
    vr_type==P1D ? "p1d" :
                   (vr_type==P3D ? "p3d" : "p2d");
  remap_data.pmid = model_pmid;
  remap_data.pint = model_pint;

  int nfields = fields.size();
  auto interp = create_interp(grid,fields);
  interp->setup_time_database(input_files,util::TimeLine::YearlyPeriodic);
  interp->setup_remappers (remap_data);
  interp->init_data_interval(t0);

  // We jump ahead by 2 months, but the shift interval logic cannot keep up with
  // a dt that long, so we should get an error due to the interpolation param being
  // outside the [0,1] interval.
  REQUIRE_THROWS (interp->run(t0+60*spd));

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

    // Just in case our testing logic is buggy, the run call below should print
    // similar information, so we can more easily debug.
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
      const auto& layout = fields[i].get_header().get_identifier().get_layout();
      bool midpoints  = layout.has_tag(LEV);
      bool interfaces = layout.has_tag(ILEV);
      bool has_levels = midpoints or interfaces;

      // Compute expected field, possibly clipping the top/bot value, if vr_type!=None
      // Since we run with 2x the number of levels as in the data, we get ONE level
      // out-of-bounds at both top and bot. The easiest way to convince yourself is to
      // make a quick diagram. Recall that in case of vremap, the input data is at MIDPOINTS
      // for both lev/ilev fields. The vfine grid has interfaces corresponding to the union
      // of midpoints and interfaces on the data grid.
      expected[i].update(base[i],1,0);
      expected[i].update(ones[i],delta,1.0);
      if ( vr_type!=NOP and has_levels) {
        auto vlen = midpoints ? nlevs : nlevs+1;
        auto scalar = not layout.is_vector_layout();

        expected_vcoarse[i].update(base_vcoarse[i],1,0);
        expected_vcoarse[i].update(ones_vcoarse[i],delta,1.0);

        expected[i].sync_to_host();
        expected_vcoarse[i].sync_to_host();
        if (scalar) {
          auto e = expected[i].get_view<Real**,Host>();
          auto e_vcoarse = expected_vcoarse[i].get_view<const Real**,Host>();

          for (int icol=0; icol<ncols; ++icol) {
            e(icol,0) = e_vcoarse(icol,0);
            e(icol,vlen-1) = e_vcoarse(icol,data_nlevs-1);
          }
        } else {
          auto e = expected[i].get_view<Real***,Host>();
          auto e_vcoarse = expected_vcoarse[i].get_view<const Real***,Host>();

          for (int icol=0; icol<ncols; ++icol) {
            for (int idim=0; idim<layout.get_vector_dim(); ++idim) {
              e(icol,idim,0) = e_vcoarse(icol,idim,0);
              e(icol,idim,vlen-1) = e_vcoarse(icol,idim,data_nlevs-1);
            }
          }
        }
        expected[i].sync_to_dev();
      }

      // Compute the rel error
      diff[i].update(fields[i],1,1);
      diff[i].update(expected[i],-1,1);
      diff[i].scale(1.0 / frobenius_norm<Real>(expected[i]));
      if (frobenius_norm<Real>(diff[i])>=tol) {
        print_field_hyperslab(expected[i]);
        print_field_hyperslab(fields[i]);
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
  auto hfine_grid  = create_point_grid("pg_h",fine_ngcols,data_nlevs,comm);
  auto vfine_grid  = create_point_grid("pg_v",data_ngcols,fine_nlevs,comm);
  auto hvfine_grid = create_point_grid("pg_hv",fine_ngcols,fine_nlevs,comm);

  SECTION ("periodic") {
    // We assume we start with t0 sometime between the first and second input slice.
    auto t_beg = reset_year(get_last_slice_time(),2019);
    auto timeline = util::TimeLine::YearlyPeriodic;
    strvec_t files = {"data_interpolation_0.nc","data_interpolation_1.nc"};
    strvec_t files_no_ilev = {"data_interpolation_0_no_ilev.nc","data_interpolation_1_no_ilev.nc"};

    SECTION ("no-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=NO ..........\n");
        run_tests (data_grid,files,t_beg,timeline);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p1d .........\n");
        run_tests (vfine_grid,files_no_ilev,t_beg,timeline,P1D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p2d .........\n");
        run_tests (vfine_grid,files_no_ilev,t_beg,timeline,P2D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=NO,  vert_remap=p3d .........\n");
        run_tests (vfine_grid,files_no_ilev,t_beg,timeline,P3D);
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
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=p1d .........\n");
        run_tests (hvfine_grid,files_no_ilev,t_beg,timeline,P1D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=p2d .........\n");
        run_tests (hvfine_grid,files_no_ilev,t_beg,timeline,P2D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=p3d .........\n");
        run_tests (hvfine_grid,files_no_ilev,t_beg,timeline,P3D);
        root_print(comm,"  timeline=PERIODIC, horiz_remap=YES, vert_remap=p3d ......... PASS\n");
      }
    }
  }

  SECTION ("linear") {
    // We assume we start with t0 sometime between the first and second input slice.
    auto t_beg = get_first_slice_time();
    auto timeline = util::TimeLine::Linear;
    strvec_t files = {"data_interpolation_0.nc","data_interpolation_1.nc"};
    strvec_t files_no_ilev = {"data_interpolation_0_no_ilev.nc","data_interpolation_1_no_ilev.nc"};

    SECTION ("no-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=NO ..........\n");
        run_tests (data_grid,files,t_beg,timeline);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=p1d .........\n");
        run_tests (vfine_grid,files_no_ilev,t_beg,timeline,P1D);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=p2d .........\n");
        run_tests (vfine_grid,files_no_ilev,t_beg,timeline,P2D);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=p3d .........\n");
        run_tests (vfine_grid,files_no_ilev,t_beg,timeline,P3D);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=NO,  vert_remap=p3d ......... PASS\n");
      }
    }

    SECTION ("yes-horiz") {
      SECTION ("no-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=NO ..........\n");
        run_tests (hfine_grid,files,t_beg,timeline);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=NO .......... PASS\n");
      }
      SECTION ("p1d-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=p1d .........\n");
        run_tests (hvfine_grid,files_no_ilev,t_beg,timeline,P1D);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=p1d ......... PASS\n");
      }
      SECTION ("p2d-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=p2d .........\n");
        run_tests (hvfine_grid,files_no_ilev,t_beg,timeline,P2D);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=p2d ......... PASS\n");
      }
      SECTION ("p3d-vert") {
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=p3d .........\n");
        run_tests (hvfine_grid,files_no_ilev,t_beg,timeline,P3D);
        root_print(comm,"  timeline=LINEAR,   horiz_remap=YES, vert_remap=p3d ......... PASS\n");
      }
    }
  }

  scorpio::finalize_subsystem();
}

} // anonymous namespace
