#include <catch2/catch.hpp>
#include <iostream>
#include <fstream> 

#include "ekat/util/ekat_string_utils.hpp"
#include "scream_config.h"
#include "share/scream_types.hpp"

#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_parameter_list.hpp"
namespace {
using namespace scream;
using namespace ekat::units;
using input_type = AtmosphereInput;

std::shared_ptr<FieldManager<Real>>
get_test_fm(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<FieldManager<Real>>
backup_fm(const std::shared_ptr<FieldManager<Real>>& src);

std::shared_ptr<UserProvidedGridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);

ekat::ParameterList get_om_params(const ekat::Comm& comm, const bool check);

ekat::ParameterList
get_in_params(const std::string& type, const ekat::Comm& comm);

void randomize_fields (const FieldManager<Real>& fm, const int seed);

void time_advance (const FieldManager<Real>& fm,
                   const std::list<ekat::CaseInsensitiveString>& fnames,
                   const double dt);

TEST_CASE("restart","io")
{
  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");

  // Note to AaronDonahue:  You are trying to figure out why you can't change the number of cols and levs for this test.  
  // Something having to do with freeing up and then resetting the io_decompositions.
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_gcols = 2*io_comm.size();
  Int num_levs = 3;

  // First set up a field manager and grids manager to interact with the output functions
  auto grid_man = get_test_gm(io_comm,num_gcols,num_levs);
  auto grid = grid_man->get_grid("Physics");
  auto field_manager = get_test_fm(grid);
  randomize_fields(*field_manager,seed);

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Create an Output manager for testing output
  OutputManager output_manager;
  auto output_params = get_om_params(io_comm,false);
  output_manager.setup(io_comm,output_params,field_manager,false);

  // Construct a timestamp
  util::TimeStamp time (0,0,0,0);
  const auto& out_fields = field_manager->get_groups_info().at("output")->m_fields_names;

  // We advance the fields, by adding dt to each entry of the fields at each time step
  // The restart data is written every 5 time steps, while the output freq is 10.
  // We run for 17 steps, which means that after 17 steps we should have both a restart
  // and a history restart files, with solution at t=15. 

  // Time-advance all fields
  const Real dt = 1.0;
  const int nsteps = 15;
  for (int i=0; i<nsteps; ++i) {
    time_advance(*field_manager,out_fields,dt);
    time += dt;
    output_manager.run(time);
  }

  // Create a copy of the field manager current status, for checking restart
  auto fm_15 = backup_fm(field_manager);

  // Restart the simulation from t=15. The restart group should match
  // what was in the fm at t=15, while field_3, which is not in the restart
  // group, should be 
  auto fm_res = get_test_fm(grid);
  auto res_params = get_in_params("Restart",io_comm);
  input_type ins_input(io_comm,res_params,fm_res);
  ins_input.read_variables();

  // Restart group fields must all match
  for (const auto& fn : {"field_1", "field_2", "field_4"}) {
    auto res_f = field_manager->get_field(fn);
    auto tgt_f = fm_15->get_field(fn);

    // The stored values must match
    REQUIRE ( views_are_equal(res_f,tgt_f) );

    // The time stamps must match
    const auto& res_ts = res_f.get_header().get_tracking().get_time_stamp();
    const auto& tgt_ts = tgt_f.get_header().get_tracking().get_time_stamp();
    REQUIRE (res_ts==tgt_ts);
  }
  // "field_3" is not in restart group, so it should contain all -1
  // Create a dummy field of equal layout, and simply call 'views_are_equal'
  auto f3 = fm_res->get_field("field_3");
  Field<Real> f_check(f3.get_header().get_identifier());
  f_check.allocate_view();
  f_check.deep_copy(-1.0);
  REQUIRE (views_are_equal(f_check,f3));

  // THIS IS HACKY BUT VERY IMPORTANT!
  // E3SM relies on the 'rpointer.atm' file to write/read the name of the restart files.
  // As of this point, rpointer contains the restart info for the timestep 15.
  // But when we run the next 5 time steps, we will reach another restart checkpoint,
  // at which point the rpointer file will be updated, and the info about the
  // restart files at timestep 15 will be lost.
  // To overcome this, we open the rpointer fiile NOW, store its content in a string,
  // run the next 5 timesteps, and then OVERWRITE the rpointer file with the content
  // we saved from the timestep=15 one.
  std::ifstream rpointer_file_in;
  rpointer_file_in.open("rpointer.atm");
  std::string content, line;
  while (rpointer_file_in >> line) {
    content += line + "\n";
  }
  rpointer_file_in.close();

  // Continue initial simulation for 5 more steps, to get to the next output step
  for (int i=0; i<5; ++i) {
    time_advance(*field_manager,out_fields,dt);
    time += dt;
    output_manager.run(time);
  }
  output_manager.finalize();

  // Restore the rpointer file as it was after timestep=15
  std::ofstream rpointer_file_out;
  rpointer_file_out.open("rpointer.atm", std::ios_base::trunc | std::ios_base::out);
  rpointer_file_out << content;
  rpointer_file_out.close();

  // Creating a new scorpio output (via the output manager) should
  // restart the output from the saved history.

  // Create Output manager, and read the restart
  util::TimeStamp time_res (0,0,0,15);
  auto output_params_res = get_om_params(io_comm,true);
  OutputManager output_manager_res;
  output_manager_res.setup(io_comm,output_params_res,fm_res,true);

  // Run 5 more steps from the restart, to get to the next output step.
  // We should be generating the same output file as before.
  for (int i=0; i<5; ++i) {
    time_advance(*fm_res,out_fields,dt);
    time_res += dt;
    output_manager_res.run(time_res);
  }
  output_manager_res.finalize();

  // Finalize everything
  scorpio::eam_pio_finalize();
  grid_man->clean_up();
} // TEST_CASE restart

/*===================================================================================================================*/
std::shared_ptr<FieldManager<Real>> get_test_fm(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;
  using FR = FieldRequest;
  using SL = std::list<std::string>;

  // Create a fm
  auto fm = std::make_shared<FieldManager<Real>>(grid);

  const int num_lcols = grid->get_num_local_dofs();
  const int num_levs = grid->get_num_vertical_levels();

  // Create some fields for this fm
  std::vector<FieldTag> tag_h  = {COL};
  std::vector<FieldTag> tag_v  = {LEV};
  std::vector<FieldTag> tag_2d = {COL,LEV};
  std::vector<FieldTag> tag_3d = {COL,CMP,LEV};

  std::vector<Int>     dims_h  = {num_lcols};
  std::vector<Int>     dims_v  = {num_levs};
  std::vector<Int>     dims_2d = {num_lcols,num_levs};
  std::vector<Int>     dims_3d = {num_lcols,2,num_levs};

  const std::string& gn = grid->name();

  FieldIdentifier fid1("field_1",FL{tag_h,dims_h},m,gn);
  FieldIdentifier fid2("field_2",FL{tag_v,dims_v},kg,gn);
  FieldIdentifier fid3("field_3",FL{tag_2d,dims_2d},kg/m,gn);
  FieldIdentifier fid4("field_4",FL{tag_3d,dims_3d},kg/m,gn);

  // Register fields with fm
  fm->registration_begins();
  fm->register_field(FR{fid1,SL{"output","restart"}});
  fm->register_field(FR{fid2,SL{"output","restart"}});
  fm->register_field(FR{fid3,"output"});
  fm->register_field(FR{fid4,SL{"output","restart"}});
  fm->registration_ends();

  // Initialize fields to -1.0, and set initial time stamp
  util::TimeStamp time (0,0,0,0);
  fm->init_fields_time_stamp(time);
  for (const auto& fn : {"field_1","field_2","field_3","field_4"} ) {
    fm->get_field(fn).deep_copy(-1.0);
  }

  return fm;
}

std::shared_ptr<FieldManager<Real>>
backup_fm (const std::shared_ptr<FieldManager<Real>>& src_fm)
{
  // Now, create a copy of the field manager current status, for comparisong
  auto dst_fm = get_test_fm(src_fm->get_grid());
  for (const auto& fn : {"field_1","field_2","field_3","field_4"} ) {
          auto f_dst = dst_fm->get_field(fn);
    const auto f_src = src_fm->get_field(fn);
    f_dst.deep_copy(f_src);

    auto src_ts = f_src.get_header().get_tracking().get_time_stamp();
    f_dst.get_header().get_tracking().update_time_stamp(src_ts);
  }
  return dst_fm;
}

/*===================================================================================================================*/
void randomize_fields (const FieldManager<Real>& fm, const int seed)
{
  using rngAlg = std::mt19937_64;
  rngAlg engine(seed);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.01,0.99);

  // Initialize these fields
  const auto& f1 = fm.get_field("field_1");
  const auto& f2 = fm.get_field("field_2");
  const auto& f3 = fm.get_field("field_3");
  const auto& f4 = fm.get_field("field_4");
  randomize(f1,engine,pdf);
  randomize(f2,engine,pdf);
  randomize(f3,engine,pdf);
  randomize(f4,engine,pdf);
}

/*===================================================================================================================*/
std::shared_ptr<UserProvidedGridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs)
{

  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("User Provided",create_user_provided_grids_manager);
  auto dummy_grid = create_point_grid("Physics",num_gcols,num_levs,io_comm);
  auto upgm = std::make_shared<UserProvidedGridsManager>();
  upgm->set_grid(dummy_grid);

  return upgm;
}
/*===================================================================================================================*/
ekat::ParameterList get_om_params(const ekat::Comm& comm, const bool check)
{
  ekat::ParameterList om_params("Output Manager");
  om_params.set<Int>("PIO Stride",1);

  std::vector<std::string> fileNames (1);
  fileNames[0] = std::string("io_test_restart") + (check ? "_check" : "") + "_np" + std::to_string(comm.size()) + ".yaml";
  om_params.set("Output YAML Files",fileNames);
  auto& res_sub = om_params.sublist("Restart Control");
  auto& freq_sub = res_sub.sublist("FREQUENCY");
  freq_sub.set<Int>("OUT_N",5);
  freq_sub.set<std::string>("OUT_OPTION","Steps");

  return om_params;  
}
/*===================================================================================================================*/
ekat::ParameterList get_in_params(const std::string& type, const ekat::Comm& comm)
{
  ekat::ParameterList in_params("Input Parameters");
  std::string filename;
  if (type=="Restart")
  {
    std::ifstream rpointer_file;
    rpointer_file.open("rpointer.atm");
    bool found = false;
    while (rpointer_file >> filename)
    {
      if (filename.find(".r.") != std::string::npos)
      {
        found = true;
        break;
      }
    }
    EKAT_REQUIRE_MSG(found,"ERROR! rpointer.atm file does not contain a restart file."); 
  
    std::vector<std::string> fnames = {"field_1", "field_2", "field_4"};
    in_params.set("FIELDS",fnames);
  } else if (type=="Final") {
    filename = "io_output_restart_np" + std::to_string(comm.size()) + ".AVERAGE.Steps_x10.0000-01-01.000020.nc";
    std::vector<std::string> fnames = {"field_1", "field_2", "field_3", "field_4"};
    in_params.set("FIELDS",fnames);
  }
  in_params.set<std::string>("FILENAME",filename);
  in_params.set<std::string>("GRID","Physics");
  return in_params;
}
/*===================================================================================================================*/

void time_advance (const FieldManager<Real>& fm,
                   const std::list<ekat::CaseInsensitiveString>& fnames,
                   const double dt) {
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
                v(i,j,k) += dt;
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

} // undefined namespace
