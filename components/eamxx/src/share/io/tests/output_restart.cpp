#include <catch2/catch.hpp>

#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util//scream_setup_random_test.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <iostream>
#include <fstream>

namespace scream {

std::shared_ptr<FieldManager>
get_test_fm(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);

template<typename Engine>
void randomize_fields (const FieldManager& fm, Engine& engine);

void time_advance (const FieldManager& fm,
                   const std::list<ekat::CaseInsensitiveString>& fnames,
                   const int dt);

std::shared_ptr<FieldManager>
backup_fm (const std::shared_ptr<FieldManager>& src_fm);

TEST_CASE("output_restart","io")
{
  // Note to AaronDonahue:  You are trying to figure out why you can't change the number of cols and levs for this test.  
  // Something having to do with freeing up and then resetting the io_decompositions.
  ekat::Comm io_comm(MPI_COMM_WORLD);
  Int num_gcols = 2*io_comm.size();
  Int num_levs = 3;

  auto engine = setup_random_test(&io_comm);

  // First set up a field manager and grids manager to interact with the output functions
  auto gm = get_test_gm(io_comm,num_gcols,num_levs);
  auto grid = gm->get_grid("Point Grid");
  auto field_manager = get_test_fm(grid);
  randomize_fields(*field_manager,engine);
  const auto& out_fields = field_manager->get_groups_info().at("output")->m_fields_names;

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());
  scorpio::eam_init_pio_subsystem(fcomm);

  // Timestamp of the simulation initial time
  util::TimeStamp t0 ({2000,1,1},{0,0,0});

  // Create an Output manager for testing output
  std::string param_filename = "io_test_restart.yaml";
  ekat::ParameterList output_params;
  ekat::parse_yaml_file(param_filename,output_params);
  output_params.set<std::string>("Floating Point Precision","real");
  OutputManager output_manager;
  output_manager.setup(io_comm,output_params,field_manager,gm,t0,t0,false);

  // We advance the fields, by adding dt to each entry of the fields at each time step
  // The output restart data is written every 5 time steps, while the output freq is 10.
  // We run for 15 steps, which means that after 15 steps we should have a history restart
  // file, with output history right in the middle between two output steps.

  // Time-advance all fields
  const int dt = 1;
  const int nsteps = 15;
  auto time = t0;
  for (int i=0; i<nsteps; ++i) {
    time_advance(*field_manager,out_fields,dt);
    time += dt;
    output_manager.run(time);
  }

  // THIS IS HACKY BUT VERY IMPORTANT!
  // E3SM relies on the 'rpointer.atm' file to write/read the name of the model/output
  // restart files. As of this point, rpointer contains the restart info for the timestep 15.
  // But when we run the next 5 time steps, we will reach another checkpoint step,
  // at which point the rpointer file will be updated, and the info about the
  // restart files at timestep 15 will be lost.
  // To overcome this, we open the rpointer fiile NOW, store its content in a string,
  // run the next 5 timesteps, and then OVERWRITE the rpointer file with the content
  // we saved from the timestep=15 one.
  std::string rpointer_content;
  if (io_comm.am_i_root()) {
    std::ifstream rpointer_file_in;
    rpointer_file_in.open("rpointer.atm");
    std::string line;
    while (rpointer_file_in >> line) {
      rpointer_content += line + "\n";
    }
    rpointer_file_in.close();
  }

  // Create a copy of the FM at the current state (used later to emulate a "restarted state")
  auto fm_res = backup_fm(field_manager);

  // Continue initial simulation for 5 more steps, to get to the next output step
  for (int i=0; i<5; ++i) {
    time_advance(*field_manager,out_fields,dt);
    time += dt;
    output_manager.run(time);
  }
  output_manager.finalize();

  // Restore the rpointer file as it was after timestep=15
  if (io_comm.am_i_root()) {
    std::ofstream rpointer_file_out;
    rpointer_file_out.open("rpointer.atm", std::ios_base::trunc | std::ios_base::out);
    rpointer_file_out << rpointer_content;
    rpointer_file_out.close();
  }

  // Now we redo the timesteps 16-20 with a fresh new output manager,
  // but we specify to the output manager that this is a restarted simulation.
  // NOTE: we use fm_res (the copy of field_manager at t=15), since we don't want
  //       to have to do a state restart, which would/could mix issues related to
  //       model restart with the testing of output history restart.
  util::TimeStamp time_res ({2000,1,1},{0,0,15},5);
  std::string param_filename_res = "io_test_restart_check.yaml";

  ekat::ParameterList output_params_res;
  ekat::parse_yaml_file(param_filename_res,output_params_res);
  output_params_res.set<std::string>("Floating Point Precision","real");

  OutputManager output_manager_res;
  output_manager_res.setup(io_comm,output_params_res,fm_res,gm,time_res,t0,false);

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
} 

/*=============================================================================================*/
std::shared_ptr<FieldManager> get_test_fm(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  using FL = FieldLayout;
  using FR = FieldRequest;
  using SL = std::list<std::string>;

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid);

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
  fm->register_field(FR{fid1,SL{"output"}});
  fm->register_field(FR{fid2,SL{"output"}});
  fm->register_field(FR{fid3,SL{"output"}});
  fm->register_field(FR{fid4,SL{"output"}});
  fm->registration_ends();

  // Initialize fields to -1.0, and set initial time stamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);
  for (const auto& fn : {"field_1","field_2","field_3","field_4"} ) {
    fm->get_field(fn).deep_copy(-1.0);
    fm->get_field(fn).sync_to_host();
  }

  return fm;
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
  randomize(f1,engine,pdf);
  randomize(f2,engine,pdf);
  randomize(f3,engine,pdf);
  randomize(f4,engine,pdf);
}

/*=============================================================================================*/
std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs)
{
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",num_gcols);
  gm_params.set("number_of_vertical_levels",num_levs);
  auto gm = create_mesh_free_grids_manager(io_comm,gm_params);
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

std::shared_ptr<FieldManager>
backup_fm (const std::shared_ptr<FieldManager>& src_fm)
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


} // namespace scream
