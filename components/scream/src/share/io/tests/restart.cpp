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
#include "share/field/field_manager.hpp"

#include "ekat/ekat_parameter_list.hpp"
namespace {
using namespace scream;
using namespace ekat::units;
using input_type = AtmosphereInput;

std::shared_ptr<FieldManager<Real>>    get_test_fm(std::shared_ptr<const AbstractGrid> grid);
std::shared_ptr<UserProvidedGridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);
ekat::ParameterList                       get_om_params(const Int casenum, const ekat::Comm& comm);
ekat::ParameterList                       get_in_params(const std::string& type, const ekat::Comm& comm);
void                                      Initialize_field_manager(const FieldManager<Real>& fm, const Int num_lcols, const Int num_levs);

TEST_CASE("restart","io")
{
  // Note to AaronDonahue:  You are trying to figure out why you can't change the number of cols and levs for this test.  
  // Something having to do with freeing up and then resetting the io_decompositions.
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_gcols = 2*io_comm.size();
  Int num_levs = 2;

  // First set up a field manager and grids manager to interact with the output functions
  auto grid_man = get_test_gm(io_comm,num_gcols,num_levs);
  auto grid = grid_man->get_grid("Physics");
  int num_lcols = grid->get_num_local_dofs();
  auto field_manager = get_test_fm(grid);

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Create an Output manager for testing output
  OutputManager m_output_manager;
  auto output_params = get_om_params(2,io_comm);
  m_output_manager.set_params(output_params);
  m_output_manager.set_comm(io_comm);
  m_output_manager.set_grids(grid_man);
  m_output_manager.set_field_mgr(field_manager);
  m_output_manager.init();

  // Construct a timestamp
  util::TimeStamp time (0,0,0,0);

  //  Cycle through data and write output
  const auto& out_fields = field_manager->get_groups_info().at("output");
  Int max_steps = 17;  // Go a few steps past the last restart write to make sure that the last written file is on the 15th step.
  Real dt = 1.0;
  for (Int ii=0;ii<max_steps;++ii) {
    for (const auto& fname : out_fields->m_fields_names) {
      auto f  = field_manager->get_field(fname);
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
  m_output_manager_res.set_field_mgr(field_manager);
  m_output_manager_res.set_runtype_restart(true);
  m_output_manager_res.init();
  auto res_params = get_in_params("Restart",io_comm);
  // reinit fields
  Initialize_field_manager(*field_manager,num_lcols,num_levs);
  // grab restart data
  input_type ins_input(io_comm,res_params,field_manager,grid_man);
  ins_input.pull_input();
  // Note, that only field_1, field_2, and field_4 were marked for restart.  Check to make sure values
  // in the field manager reflect those fields as restarted from 15 and field_3 as being
  // freshly restarted:
  auto field1 = field_manager->get_field("field_1");
  auto field2 = field_manager->get_field("field_2");
  auto field3 = field_manager->get_field("field_3");
  auto field4 = field_manager->get_field("field_4");
  auto field1_hst = field1.get_view<Real*,Host>();
  auto field2_hst = field2.get_view<Real*,Host>();
  auto field3_hst = field3.get_view<Real**,Host>();
  auto field4_hst = field4.get_view<Real***,Host>();
  field1.sync_to_host();
  field2.sync_to_host();
  field3.sync_to_host();
  field4.sync_to_host();
  Real tol = pow(10,-6);
  for (Int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(field1_hst(ii)-(15+ii))<tol);
    for (Int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(field2_hst(jj)   -((jj+1)/10.+15))<tol);
      REQUIRE(std::abs(field3_hst(ii,jj)-((jj+1)/10.+ii))<tol);  //Note, field 3 is not restarted, so doesn't have the +15
      std::cout << "f  : " << field4_hst(ii,0,jj) << "\n";
      std::cout << "tgt: " << ((jj+1)/10.+ii+15) << "\n";
      REQUIRE(std::abs(field4_hst(ii,0,jj)-((jj+1)/10.+ii+15))<tol);
      REQUIRE(std::abs(field4_hst(ii,1,jj)-(-((jj+1)/10.+ii)+15))<tol);
    }
  }
  // Finish the last 5 steps
  for (Int ii=0;ii<5;++ii) {
    for (const auto& fname : out_fields->m_fields_names) {
      auto f = field_manager->get_field(fname);
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
  input_type avg_input(io_comm,avg_params,field_manager,grid_man);
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

  // Initialize these fields
  Initialize_field_manager(*fm,num_lcols,num_levs);
  // Update timestamp
  util::TimeStamp time (0,0,0,0);
  fm->init_fields_time_stamp(time);

  return fm;
}
/*===================================================================================================================*/
void Initialize_field_manager(const FieldManager<Real>& fm, const Int num_lcols, const Int num_levs)
{

  // Initialize these fields
  const auto& f1 = fm.get_field("field_1");
  const auto& f2 = fm.get_field("field_2");
  const auto& f3 = fm.get_field("field_3");
  const auto& f4 = fm.get_field("field_4");
  auto f1_hst = f1.get_view<Real*, Host>();
  auto f2_hst = f2.get_view<Real*, Host>();
  auto f3_hst = f3.get_view<Real**,Host>();
  auto f4_hst = f4.get_view<Real***,Host>();
  for (int ii=0;ii<num_lcols;++ii) {
    f1_hst(ii) = ii;
    for (int jj=0;jj<num_levs;++jj) {
      f2_hst(jj) = (jj+1)/10.0;
      f3_hst(ii,jj) = (ii) + (jj+1)/10.0;
      f4_hst(ii,0,jj) = ii + (jj+1)/10.0;
      f4_hst(ii,1,jj) = -( ii + (jj+1)/10.0 );
    }
  }
  f1.sync_to_dev();
  f2.sync_to_dev();
  f3.sync_to_dev();
  f4.sync_to_dev();
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
  
    std::vector<std::string> fnames = {"field_1", "field_2", "field_4"};
    in_params.set("FIELDS",fnames);
  } else if (type=="Final") {
    filename = "io_output_restart_np" + std::to_string(comm.size()) + ".Average.Steps_x10.0000-01-01.000020.nc";
    std::vector<std::string> fnames = {"field_1", "field_2", "field_3", "field_4"};
    in_params.set("FIELDS",fnames);
  }
  in_params.set<std::string>("FILENAME",filename);
  in_params.set<std::string>("GRID","Physics");
  return in_params;
}
/*===================================================================================================================*/
} // undefined namespace
