#include <catch2/catch.hpp>

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util//eamxx_setup_random_test.hpp"

#include "share/eamxx_types.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>

namespace scream {

constexpr Real FillValue = constants::DefaultFillValue<float>().value;

std::shared_ptr<FieldManager>
get_test_fm(const std::shared_ptr<const AbstractGrid>& grid);

std::shared_ptr<FieldManager>
clone_fm (const std::shared_ptr<const FieldManager>& fm);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& comm, const Int num_gcols, const Int num_levs);

template<typename Engine>
void randomize_fields (const FieldManager& fm, Engine& engine);

void time_advance (const FieldManager& fm,
                   const std::list<ekat::CaseInsensitiveString>& fnames,
                   const int dt);

TEST_CASE("output_restart","io")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  // If running with 2+ ranks, this will check that restart works correctly
  // even if some ranks own no dofs
  int num_gcols = std::max(comm.size()-1,1);
  int num_levs = 3;
  int dt = 1;

  auto engine = setup_random_test(&comm);

  // First set up a field manager and grids manager to interact with the output functions
  auto gm = get_test_gm(comm,num_gcols,num_levs);
  auto grid = gm->get_grid("Point Grid");

  // The the IC field manager
  auto fm0 = get_test_fm(grid);
  randomize_fields(*fm0,engine);

  const auto& out_fields = fm0->get_groups_info().at("output")->m_fields_names;

  // Initialize the pio_subsystem for this test:
  scorpio::init_subsystem(comm);

  // Timestamp of the simulation initial time
  util::TimeStamp t0 ({2000,1,1},{0,0,0});

  // Create output params (some options are set below, depending on the run type
  ekat::ParameterList output_params;
  output_params.set<std::string>("Floating Point Precision","real");
  output_params.set<std::vector<std::string>>("Field Names",{"field_1", "field_2", "field_3", "field_4","field_5"});
  output_params.set<double>("fill_value",FillValue);
  output_params.set<int>("flush_frequency",1);
  output_params.sublist("output_control").set<std::string>("frequency_units","nsteps");
  output_params.sublist("output_control").set<int>("Frequency",10);
  output_params.sublist("Checkpoint Control").set<int>("Frequency",5);
  // This skips a test that only matters for AD runs
  output_params.sublist("Checkpoint Control").set<bool>("is_unit_testing","true");

  // Creates and runs an OM from output_params and given inputs
  auto run = [&](std::shared_ptr<FieldManager> fm,
                 const util::TimeStamp& case_t0,
                 const util::TimeStamp& run_t0,
                 const int nsteps)
  {
    OutputManager output_manager;
    output_manager.initialize(comm, output_params, run_t0, case_t0, false);
    output_manager.setup(fm,gm);

    // We advance the fields, by adding dt to each entry of the fields at each time step
    // The output restart data is written every 5 time steps, while the output freq is 10.
    auto time = run_t0;
    for (int i=0; i<nsteps; ++i) {
      output_manager.init_timestep(time,dt);
      time_advance(*fm,out_fields,dt);
      time += dt;
      output_manager.run(time);
    }
    output_manager.finalize();
  };

  auto print = [&] (const std::string& s, int line_len = -1) {
    if (comm.am_i_root()) {
      if (line_len<0) {
        std::cout << s;
      } else {
        std::cout << std::left << std::setw(line_len) << std::setfill('.') << s;
      }
    }
  };
  // Run test for different avg type choices
  for (const std::string avg_type : {"INSTANT","AVERAGE"}) {
    {
      // In normal runs, the OM for the model restart takes care of nuking rpointer.atm,
      // and re-creating a new one. Here, we don't have that, so we must nuke it manually
      std::ofstream ofs;
      ofs.open("rpointer.atm", std::ofstream::out | std::ofstream::trunc);
    }
    print("   -> Averaging type: " + avg_type + " ", 40);
    output_params.set<std::string>("Averaging Type",avg_type);

    // 1. Run for full 20 days, no restarts needed
    auto fm_mono = clone_fm(fm0);
    output_params.set<std::string>("filename_prefix","monolithic");
    output_params.sublist("Checkpoint Control").set<std::string>("frequency_units","never");
    run(fm_mono,t0,t0,20);

    // 2. Run for 15 days on fm0, write restart every 5 steps
    auto fm_rest = clone_fm(fm0);
    output_params.set<std::string>("filename_prefix","restarted");
    output_params.sublist("Checkpoint Control").set<std::string>("frequency_units","nsteps");
    run(fm_rest,t0,t0,15);

    // 3. Restart the second run at step=15, and do 5 more steps
    // NOTE: keep fm_rest FM, since we are not testing the restart of the state, just the history.
    //       Here, we proceed as if the AD already restarted the state correctly.
    output_params.sublist("Checkpoint Control").set<std::string>("frequency_units","never");

    // Ensure nsteps is equal to 15 upon restart
    auto run_t0 = (t0+15*dt).clone(15);
    run(fm_rest,t0,run_t0,5);
    print(" DONE\n");
  }
  // Finalize everything
  scorpio::finalize_subsystem();
}

/*=============================================================================================*/
std::shared_ptr<FieldManager>
get_test_fm(const std::shared_ptr<const AbstractGrid>& grid)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  using FR = FieldRequest;
  using SL = std::list<std::string>;

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid);

  auto scalar_1d = grid->get_vertical_layout(true);
  auto scalar_2d = grid->get_2d_scalar_layout();
  auto scalar_3d = grid->get_3d_scalar_layout(true);
  auto vector_3d = grid->get_3d_vector_layout(true,2);
  auto rad_vector_3d = grid->get_3d_vector_layout(true,3,"SWBND");

  const std::string& gn = grid->name();

  FieldIdentifier fid1("field_1",scalar_2d,    m,   gn);
  FieldIdentifier fid2("field_2",scalar_1d,    kg,  gn);
  FieldIdentifier fid3("field_3",scalar_3d,    kg/m,gn);
  FieldIdentifier fid4("field_4",vector_3d,    kg/m,gn);
  FieldIdentifier fid5("field_5",rad_vector_3d,m*m, gn);

  // Register fields with fm
  fm->registration_begins();
  fm->register_field(FR{fid1,SL{"output"}});
  fm->register_field(FR{fid2,SL{"output"}});
  fm->register_field(FR{fid3,SL{"output"}});
  fm->register_field(FR{fid4,SL{"output"}});
  fm->register_field(FR{fid5,SL{"output"}});
  fm->registration_ends();

  // Initialize fields to -1.0, and set initial time stamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);
  for (const auto& fn : {"field_1","field_2","field_3","field_4","field_5"} ) {
    fm->get_field(fn).deep_copy(-1.0);
    fm->get_field(fn).sync_to_host();
  }

  return fm;
}

std::shared_ptr<FieldManager>
clone_fm(const std::shared_ptr<const FieldManager>& src) {
  auto copy = std::make_shared<FieldManager>(src->get_grid());
  copy->registration_begins();
  copy->registration_ends();
  for (auto it : *src) {
    copy->add_field(it.second->clone());
  }

  return copy;
}

/*=================================================================================================*/
template<typename Engine>
void randomize_fields (const FieldManager& fm, Engine& engine)
{
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.01,0.99);

  // Initialize these fields
  const auto& f1 = fm.get_field("field_1");
  const auto& f2 = fm.get_field("field_2");
  const auto& f3 = fm.get_field("field_3");
  const auto& f4 = fm.get_field("field_4");
  const auto& f5 = fm.get_field("field_5");
  randomize(f1,engine,pdf);
  randomize(f2,engine,pdf);
  randomize(f3,engine,pdf);
  randomize(f4,engine,pdf);
  randomize(f5,engine,pdf);
}

/*=============================================================================================*/
std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& comm, const Int num_gcols, const Int num_levs)
{
  auto gm = create_mesh_free_grids_manager(comm,0,0,num_levs,num_gcols);
  gm->build_grids();
  return gm;
}
/*===================================================================================================*/
void time_advance (const FieldManager& fm,
                   const std::list<ekat::CaseInsensitiveString>& fnames,
                   const int dt) {
  for (const auto& fname : fnames) {
    auto f  = fm.get_field(fname);
    f.sync_to_host();
    auto fl = f.get_header().get_identifier().get_layout();
    switch (fl.rank()) {
      case 1:
        {
          auto v = f.get_view<Real*,Host>();
          for (int i=0; i<fl.dim(0); ++i) {
            v(i) += dt;
          }
        }
        break;
      case 2:
        {
          auto v = f.get_view<Real**,Host>();
          for (int i=0; i<fl.dim(0); ++i) {
            for (int j=0; j<fl.dim(1); ++j) {
              v(i,j) += dt;
            }
          }
        }
        break;
      case 3:
        {
          auto v = f.get_view<Real***,Host>();
          for (int i=0; i<fl.dim(0); ++i) {
            for (int j=0; j<fl.dim(1); ++j) {
              for (int k=0; k<fl.dim(2); ++k) {
                if (fname == "field_5") {
                  // field_5 is used to test restarts w/ filled values, so
                  // we cycle between filled and unfilled states.
                  v(i,j,k) = (v(i,j,k)==FillValue) ? dt :
                    ( (v(i,j,k)==1.0) ? 2.0*dt : FillValue );
                } else {
                              v(i,j,k) += dt;
                }
              }
            }
          }
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unexpected field of rank " + std::to_string(fl.rank()) + ".\n");
    }
    f.sync_to_dev();
    auto ft = f.get_header_ptr()->get_tracking();
    auto ts = ft.get_time_stamp();
    ft.update_time_stamp(ts+dt);
  }
}

std::shared_ptr<FieldManager>
backup_fm (const std::shared_ptr<FieldManager>& src_fm)
{
  // Now, create a copy of the field manager current status, for comparisong
  auto dst_fm = get_test_fm(src_fm->get_grid());
  for (const auto& fn : {"field_1","field_2","field_3","field_4","field_5"} ) {
          auto f_dst = dst_fm->get_field(fn);
    const auto f_src = src_fm->get_field(fn);
    f_dst.deep_copy(f_src);

    auto src_ts = f_src.get_header().get_tracking().get_time_stamp();
    f_dst.get_header().get_tracking().update_time_stamp(src_ts);
  }
  return dst_fm;
}


} // namespace scream
