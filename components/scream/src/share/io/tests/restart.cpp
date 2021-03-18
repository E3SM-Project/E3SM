#include <catch2/catch.hpp>
#include <iostream>
#include <fstream> 

#include "scream_config.h"
#include "share/scream_types.hpp"

#include "share/io/output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"

#include "ekat/ekat_parameter_list.hpp"
namespace {
using namespace scream;
using namespace ekat::units;
using input_type = AtmosphereInput;

std::shared_ptr<FieldRepository<Real>>    get_test_repo(const Int num_lcols, const Int num_levs);
std::shared_ptr<UserProvidedGridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);
ekat::ParameterList                       get_om_params(const Int casenum, const ekat::Comm& comm);
ekat::ParameterList                       get_in_params(const std::string& type, const ekat::Comm& comm);
void                                      Initialize_field_repo(const FieldRepository<Real>& repo, const Int num_lcols, const Int num_levs);

TEST_CASE("restart","io")
{
  // Note to AaronDonahue:  You are trying to figure out why you can't change the number of cols and levs for this test.  
  // Something having to do with freeing up and then resetting the io_decompositions.
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_gcols = 2*io_comm.size();
  Int num_levs = 2;

  // First set up a field manager and grids manager to interact with the output functions
  std::shared_ptr<UserProvidedGridsManager> grid_man   = get_test_gm(io_comm,num_gcols,num_levs);
  int num_lcols = grid_man->get_grid("Physics")->get_num_local_dofs();
  std::shared_ptr<FieldRepository<Real>>    field_repo = get_test_repo(num_lcols,num_levs);

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Create an Output manager for testing output
  OutputManager m_output_manager;
  auto output_params = get_om_params(2,io_comm);
  m_output_manager.set_params(output_params);
  m_output_manager.set_comm(io_comm);
  m_output_manager.set_grids(grid_man);
  m_output_manager.set_repo(field_repo);
  m_output_manager.init();

  // Construct a timestamp
  util::TimeStamp time (0,0,0,0);

  //  Cycle through data and write output
  const auto& out_fields = field_repo->get_groups_info().at("output");
  Int max_steps = 17;  // Go a few steps past the last restart write to make sure that the last written file is on the 15th step.
  Real dt = 1.0;
  for (Int ii=0;ii<max_steps;++ii) {
    for (const auto& fname : out_fields->m_fields_names) {
      auto f = field_repo->get_field(fname,"Physics");
      f.sync_to_host();
      auto f_host = f.get_view<Host>();
      for (size_t jj=0;jj<f_host.size();++jj)
      {
        f_host(jj) += dt;
      }
      f.sync_to_dev();
    }
    time += dt;
    m_output_manager.run(time);
  }
  m_output_manager.finalize();

  // At this point we should have produced 2 files, a restart and a restart history file.
  util::TimeStamp time_res (0,0,0,15);
  OutputManager m_output_manager_res;
  m_output_manager_res.set_params(output_params);
  m_output_manager_res.set_comm(io_comm);
  m_output_manager_res.set_grids(grid_man);
  m_output_manager_res.set_repo(field_repo);
  m_output_manager_res.set_runtype_restart(true);
  m_output_manager_res.init();
  auto res_params = get_in_params("Restart",io_comm);
  // reinit fields
  Initialize_field_repo(*field_repo,num_lcols,num_levs);
  // grab restart data
  input_type ins_input(io_comm,res_params,field_repo,grid_man);
  ins_input.pull_input();
  // Note, that only field_1 and field_2 were marked for restart.  Check to make sure values
  // in the field manager reflect those fields as restarted from 15 and field_3 as being
  // freshly restarted:
  auto field1 = field_repo->get_field("field_1","Physics");
  auto field2 = field_repo->get_field("field_2","Physics");
  auto field3 = field_repo->get_field("field_3","Physics");
  auto field1_hst = field1.get_view<Host>();
  auto field2_hst = field2.get_view<Host>();
  auto field3_hst = field3.get_reshaped_view<Real**,Host>();
  field1.sync_to_host();
  field2.sync_to_host();
  field3.sync_to_host();
  Real tol = pow(10,-6);
  for (Int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(field1_hst(ii)-(15+ii))<tol);
    for (Int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(field2_hst(jj)   -((jj+1)/10.+15))<tol);
      REQUIRE(std::abs(field3_hst(ii,jj)-((jj+1)/10.+ii))<tol);  //Note, field 3 is not restarted, so doesn't have the +15
    }
  }
  // Finish the last 5 steps
  for (Int ii=0;ii<5;++ii) {
    for (const auto& fname : out_fields->m_fields_names) {
      auto f = field_repo->get_field(fname,"Physics");
      f.sync_to_host();
      auto f_host = f.get_view<Host>();
      for (size_t jj=0;jj<f_host.size();++jj)
      {
        f_host(jj) += dt;
      }
      f.sync_to_dev();
    }
    time_res += dt;
    m_output_manager_res.run(time_res);
  }
  m_output_manager_res.finalize();
  
  // We have now finished running a restart that should have loaded a restart history mid-way through averaging.  The
  // final output should be stored in a file: io_output_restart.Average.Steps_x10.0000-01-01.000020.nc
  // and reflect averaged values over the last 10 steps for fields 1 and 2 which were restarted and field 3 which started
  // with the initial value.  Note that we only took 5 steps, so field 3 should reflect a a value that stored steps 11-15
  // the restart history file, but was re-initialized and run for the last 5 steps.
  auto avg_params = get_in_params("Final",io_comm);
  input_type avg_input(io_comm,avg_params,field_repo,grid_man);
  avg_input.pull_input();
  field1.sync_to_host();
  field2.sync_to_host();
  field3.sync_to_host();
  Real avg_val;
  for (Int ii=0;ii<num_lcols;++ii) {
    avg_val = (20+11)/2.+ii;
    REQUIRE(std::abs(field1_hst(ii)-avg_val)<tol);
    for (Int jj=0;jj<num_levs;++jj) {
      avg_val = (20+11)/2.+(jj+1)/10.0;
      REQUIRE(std::abs(field2_hst(jj)-avg_val)<tol);
      avg_val = (15+11)/4. + (5+1)/4. + ii + (jj+1)/10.;
      REQUIRE(std::abs(field3_hst(ii,jj)-avg_val)<tol);
    }
  }
  
  // Finalize everything
  scorpio::eam_pio_finalize();
  grid_man->clean_up();
} // TEST_CASE restart
/*===================================================================================================================*/
std::shared_ptr<FieldRepository<Real>> get_test_repo(const Int num_lcols, const Int num_levs)
{
  using namespace ShortFieldTagsNames;

  // Create a repo
  std::shared_ptr<FieldRepository<Real>>  repo = std::make_shared<FieldRepository<Real>>();
  // Create some fields for this repo
  std::vector<FieldTag> tag_h  = {COL};
  std::vector<FieldTag> tag_v  = {LEV};
  std::vector<FieldTag> tag_2d = {COL,LEV};

  std::vector<Int>     dims_h  = {num_lcols};
  std::vector<Int>     dims_v  = {num_levs};
  std::vector<Int>     dims_2d = {num_lcols,num_levs};

  using FL = FieldLayout;
  const std::string gn = "Physics";
  FieldIdentifier fid1("field_1",FL{tag_h,dims_h},m,gn);
  FieldIdentifier fid2("field_2",FL{tag_v,dims_v},kg,gn);
  FieldIdentifier fid3("field_3",FL{tag_2d,dims_2d},kg/m,gn);

  // Register fields with repo
  repo->registration_begins();
  repo->register_field(fid1,{"output","restart"});
  repo->register_field(fid2,{"output","restart"});
  repo->register_field(fid3,{"output"});
  repo->registration_ends();

  // Initialize these fields
  Initialize_field_repo(*repo,num_lcols,num_levs);
  // Update timestamp
  auto f1 = repo->get_field(fid1);
  auto f2 = repo->get_field(fid2);
  auto f3 = repo->get_field(fid3);
  util::TimeStamp time (0,0,0,0);
  repo->init_fields_time_stamp(time);

  return repo;
}
/*===================================================================================================================*/
void Initialize_field_repo(const FieldRepository<Real>& repo, const Int num_lcols, const Int num_levs)
{

  // Initialize these fields
  const auto& f1 = repo.get_field("field_1","Physics");
  const auto& f2 = repo.get_field("field_2","Physics");
  const auto& f3 = repo.get_field("field_3","Physics");
  auto f1_hst = f1.get_view<Host>();
  auto f2_hst = f2.get_view<Host>();
  auto f3_hst = f3.get_reshaped_view<Real**,Host>();
  for (int ii=0;ii<num_lcols;++ii) {
    f1_hst(ii) = ii;
    for (int jj=0;jj<num_levs;++jj) {
      f2_hst(jj) = (jj+1)/10.0;
      f3_hst(ii,jj) = (ii) + (jj+1)/10.0;
    }
  }
  f1.sync_to_dev();
  f2.sync_to_dev();
  f3.sync_to_dev();
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
ekat::ParameterList get_om_params(const Int casenum, const ekat::Comm& comm)
{
  ekat::ParameterList om_params("Output Manager");
  om_params.set<Int>("PIO Stride",1);
  if (casenum == 1)
  {
    std::vector<std::string> fileNames = { "io_test_instant","io_test_average",
                                           "io_test_max",    "io_test_min" };
    for (auto& name : fileNames) {
      name += "_np" + std::to_string(comm.size()) + ".yaml";
    }

    om_params.set("Output YAML Files",fileNames);
  } 
  else if (casenum == 2)
  {
    std::vector<std::string> fileNames = { "io_test_restart_np" + std::to_string(comm.size()) + ".yaml" };
    om_params.set("Output YAML Files",fileNames);
    auto& res_sub = om_params.sublist("Restart Control");
    auto& freq_sub = res_sub.sublist("FREQUENCY");
    freq_sub.set<Int>("OUT_N",5);
    freq_sub.set<std::string>("OUT_OPTION","Steps");
  }
  else
  {
    EKAT_REQUIRE_MSG(false,"Error, incorrect case number for get_om_params");
  }

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
  
    auto& f_list = in_params.sublist("FIELDS");
    f_list.set<Int>("Number of Fields",2);
    for (int ii=1;ii<=2+1;++ii)
    {
      f_list.set<std::string>("field "+std::to_string(ii),"field_"+std::to_string(ii));
    }
  } else if (type=="Final")
  {
    filename = "io_output_restart_np" + std::to_string(comm.size()) + ".Average.Steps_x10.0000-01-01.000020.nc";
    auto& f_list = in_params.sublist("FIELDS");
    f_list.set<Int>("Number of Fields",3);
    for (int ii=1;ii<=3+1;++ii)
    {
      f_list.set<std::string>("field "+std::to_string(ii),"field_"+std::to_string(ii));
    }
  }
  in_params.set<std::string>("FILENAME",filename);
  in_params.set<std::string>("GRID","Physics");
  return in_params;
}
/*===================================================================================================================*/
} // undefined namespace
