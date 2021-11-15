#include <catch2/catch.hpp>
#include <memory>

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_string_utils.hpp"
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

#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"
#include "scream_config.h"

#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/ekat_assert.hpp"

namespace {

using namespace scream;
using namespace ekat::units;
using input_type = AtmosphereInput;
// Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
const int packsize = 2;
using Pack         = ekat::Pack<Real,packsize>;

std::shared_ptr<FieldManager<Real>>
get_test_fm(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);

ekat::ParameterList get_in_params(const std::string type,
                                  const ekat::Comm& comm,
                                  const util::TimeStamp& t_first_write);

TEST_CASE("input_output_basic","io")
{

  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_gcols = 2*io_comm.size();
  Int num_levs = 3;

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // First set up a field manager and grids manager to interact with the output functions
  auto gm = get_test_gm(io_comm,num_gcols,num_levs);
  auto grid = gm->get_grid("Point Grid");
  int num_lcols = grid->get_num_local_dofs();
  auto field_manager = get_test_fm(grid);

  // Construct a timestamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});
  util::TimeStamp time = t0;

  std::vector<std::string> fileNames = { "io_test_instant","io_test_average",
                                         "io_test_max",    "io_test_min",
                                         "io_test_multisnap" };

  // Create an Output manager for testing output
  std::vector<OutputManager> output_managers;
  for (const auto& fname : fileNames) {
    ekat::ParameterList params;
    ekat::parse_yaml_file(fname+"_np" + std::to_string(io_comm.size()) + ".yaml",params);
    output_managers.emplace_back();
    auto& om = output_managers.back();
    om.setup(io_comm,params,field_manager,gm,t0,false,false);
    io_comm.barrier();
  }

  //  Cycle through data and write output
  const auto& out_fields = field_manager->get_groups_info().at("output");
  Int max_steps = 10;
  Int dt = 1;
  for (Int ii=0;ii<max_steps;++ii) {
    time += dt;
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
        default:
          EKAT_ERROR_MSG ("Error! Unexpected field rank.\n");
      }
      f.sync_to_dev();
    }
    for (auto& om : output_managers) {
      om.run(time);
    }
  }
  for (auto& om : output_managers) {
    om.finalize();
  }

  // Because the next test for input in this sequence relies on the same fields.  We set all field values
  // to nan to ensure no false-positive tests if a field is simply not read in as input.
  for (const auto& fname : out_fields->m_fields_names) {
    auto f = field_manager->get_field(fname);
    f.deep_copy(ekat::ScalarTraits<Real>::invalid());
  }

  // At this point we should have 5 files output:
  // 1 file each for averaged, instantaneous, min and max data.
  // And 1 file with multiple time snaps of instantaneous data.
  // Cycle through each output and make sure it is correct.
  // We can use the produced output files to simultaneously check output quality and the
  // ability to read input.
  // NOTE: single-snap file start outputing at the final time only,
  //       while multi-snap starts to output at the first time step
  auto ins_params = get_in_params("Instant",io_comm,time);
  auto avg_params = get_in_params("Average",io_comm,time);
  auto min_params = get_in_params("Min",io_comm,time);
  auto max_params = get_in_params("Max",io_comm,time);
  auto multi_params = get_in_params("Multisnap",io_comm,t0+dt);
  Real tol = 100*std::numeric_limits<Real>::epsilon();
  // TODO: Create a small nc dummy file and a separate unit test which tests all input functions.
  // Test that pio_inq_dimlen is correct, using a file from one of the above parameter lists.
  {
    auto test_filename = ins_params.get<std::string>("Filename");
    scorpio::register_file(test_filename,scorpio::Read);
    Int test_gcols_len = scorpio::get_dimlen_c2f(test_filename.c_str(),"ncol");
    REQUIRE(test_gcols_len==num_gcols);
    scorpio::eam_pio_closefile(test_filename);
  }

  auto f1 = field_manager->get_field("field_1");
  auto f2 = field_manager->get_field("field_2");
  auto f3 = field_manager->get_field("field_3");
  auto f4 = field_manager->get_field("field_packed");
  auto f1_host = f1.get_view<Real*,Host>();
  auto f2_host = f2.get_view<Real*,Host>();
  auto f3_host = f3.get_view<Real**,Host>();
  auto f4_host = f4.get_view<Real**,Host>();

  // Check instant output
  input_type ins_input(io_comm,ins_params,field_manager);
  ins_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  f4.sync_to_host();

  for (int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(f1_host(ii)-(max_steps*dt+ii))<tol);
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f3_host(ii,jj)-(ii+max_steps*dt + (jj+1)/10.))<tol);
      REQUIRE(std::abs(f4_host(ii,jj)-(ii+max_steps*dt + (jj+1)/10.))<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    REQUIRE(std::abs(f2_host(jj)-(max_steps*dt + (jj+1)/10.))<tol);
  }
  ins_input.finalize();

  // Check average output
  input_type avg_input(io_comm,avg_params,field_manager);
  avg_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  Real avg_val;
  for (int ii=0;ii<num_lcols;++ii) {
    avg_val = (max_steps+1)/2.0*dt + ii; // Sum(x0+i*dt,i=1...N) = N*x0 + dt*N*(N+1)/2, AVG = Sum/N, note x0=ii in this case
    REQUIRE(std::abs(f1_host(ii)-avg_val)<tol);
    for (int jj=0;jj<num_levs;++jj) {
      avg_val = (max_steps+1)/2.0*dt + (jj+1)/10.+ii;  //note x0=(jj+1)/10+ii in this case.
      REQUIRE(std::abs(f3_host(ii,jj)-avg_val)<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    avg_val = (max_steps+1)/2.0*dt + (jj+1)/10.;  //note x0=(jj+1)/10 in this case.
    REQUIRE(std::abs(f2_host(jj)-avg_val)<tol);
  }
  avg_input.finalize();

  // Check max output
  // The max should be equivalent to the instantaneous because this function is monotonically increasing.
  input_type max_input(io_comm,max_params,field_manager);
  max_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  for (int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(f1_host(ii)-(max_steps*dt+ii))<tol);
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f3_host(ii,jj)-(ii+max_steps*dt + (jj+1)/10.))<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    REQUIRE(std::abs(f2_host(jj)-(max_steps*dt + (jj+1)/10.))<tol);
  }
  max_input.finalize();
  // Check min output
  // The min should be equivalent to the first step because this function is monotonically increasing.
  input_type min_input(io_comm,min_params,field_manager);
  min_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  for (int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(f1_host(ii)-(dt+ii))<tol);
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f3_host(ii,jj)-(dt+ii + (jj+1)/10.))<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    REQUIRE(std::abs(f2_host(jj)-(dt+(jj+1)/10.))<tol);
  }
  min_input.finalize();

  // Check multisnap output; note, tt starts at 1 instead of 0 to follow netcdf time dimension indexing.
  input_type multi_input(io_comm,multi_params,field_manager);
  for (int tt = 1; tt<=std::min(max_steps,10); tt++) {
    multi_input.read_variables(tt);
    f1.sync_to_host();
    f2.sync_to_host();
    f3.sync_to_host();
    f4.sync_to_host();

    for (int ii=0;ii<num_lcols;++ii) {
      if(std::abs(f1_host(ii)-(tt*dt+ii))>=tol) {
        printf("f1_host(%d): %f\n",ii,f1_host(ii));
        printf("expected: %f\n",Real(tt*dt+ii));
      }

      REQUIRE(std::abs(f1_host(ii)-(tt*dt+ii))<tol);
      for (int jj=0;jj<num_levs;++jj) {
        REQUIRE(std::abs(f3_host(ii,jj)-(ii+tt*dt + (jj+1)/10.))<tol);
        REQUIRE(std::abs(f4_host(ii,jj)-(ii+tt*dt + (jj+1)/10.))<tol);
      }
    }
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f2_host(jj)-(tt*dt + (jj+1)/10.))<tol);
    }
  }
  multi_input.finalize();
  
  // All Done 
  scorpio::eam_pio_finalize();
} // TEST_CASE output_instance
/* ----------------------------------*/

/*===================================================================================================================*/
std::shared_ptr<FieldManager<Real>> get_test_fm(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;
  using FR = FieldRequest;

  // Create a fm
  auto fm = std::make_shared<FieldManager<Real>>(grid);

  const int num_lcols = grid->get_num_local_dofs();
  const int num_levs = grid->get_num_vertical_levels();

  // Create some fields for this fm
  std::vector<FieldTag> tag_h  = {COL};
  std::vector<FieldTag> tag_v  = {LEV};
  std::vector<FieldTag> tag_2d = {COL,LEV};

  std::vector<Int>     dims_h  = {num_lcols};
  std::vector<Int>     dims_v  = {num_levs};
  std::vector<Int>     dims_2d = {num_lcols,num_levs};

  const std::string& gn = grid->name();

  FieldIdentifier fid1("field_1",FL{tag_h,dims_h},m,gn);
  FieldIdentifier fid2("field_2",FL{tag_v,dims_v},kg,gn);
  FieldIdentifier fid3("field_3",FL{tag_2d,dims_2d},kg/m,gn);
  FieldIdentifier fid4("field_packed",FL{tag_2d,dims_2d},kg/m,gn);

  // Register fields with fm
  // Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
  fm->registration_begins();
  fm->register_field(FR{fid1,"output"});
  fm->register_field(FR{fid2,"output"});
  fm->register_field(FR{fid3,"output"});
  fm->register_field(FR{fid4,"output",Pack::n}); // Register field as packed
  fm->registration_ends();

  // Make sure that field 4 is in fact a packed field
  auto field4 = fm->get_field(fid4);
  auto fid4_padding = field4.get_header().get_alloc_properties().get_padding();
  REQUIRE(fid4_padding > 0);

  // Initialize these fields
  auto f1 = fm->get_field(fid1);
  auto f2 = fm->get_field(fid2);
  auto f3 = fm->get_field(fid3);
  auto f4 = fm->get_field(fid4);
  auto f1_host = f1.get_view<Real*,Host>(); 
  auto f2_host = f2.get_view<Real*,Host>(); 
  auto f3_host = f3.get_view<Real**,Host>();
  auto f4_host = f4.get_view<Pack**,Host>(); 
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  f4.sync_to_host();
  for (int ii=0;ii<num_lcols;++ii) {
    f1_host(ii) = ii;
    for (int jj=0;jj<num_levs;++jj) {
      f2_host(jj) = (jj+1)/10.0;
      f3_host(ii,jj) = (ii) + (jj+1)/10.0;
      int ipack = jj / packsize;
      int ivec  = jj % packsize;
      f4_host(ii,ipack)[ivec] = (ii) + (jj+1)/10.0;
    }
  }
  // Update timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);
  // Sync back to device
  f1.sync_to_dev();
  f2.sync_to_dev();
  f3.sync_to_dev();
  f4.sync_to_dev();

  return fm;
}
/*===================================================================================================================*/
std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs)
{
  ekat::ParameterList gm_params;
  gm_params.sublist("Mesh Free").set("Number of Global Columns",num_gcols);
  gm_params.sublist("Mesh Free").set("Number of Vertical Levels",num_levs);
  auto gm = create_mesh_free_grids_manager(io_comm,gm_params);
  gm->build_grids(std::set<std::string>{"Point Grid"});
  return gm;
}
/*===================================================================================================================*/
ekat::ParameterList get_in_params(const std::string type,
                                  const ekat::Comm& comm,
                                  const util::TimeStamp& t_first_write)
{
  using vos_type = std::vector<std::string>;
  ekat::ParameterList in_params("Input Parameters");
  bool multisnap = type=="Multisnap";

  std::string filename =
        "io_" + std::string(multisnap ? "multisnap" : "output")
      + "_test_np" + std::to_string(comm.size())
      + "." + (multisnap ? ekat::upper_case("Instant") : ekat::upper_case(type))
      + ".Steps_x1" + std::string(multisnap ? "" : "0")
      + "." + t_first_write.to_string() + ".nc";

  in_params.set<std::string>("Filename",filename);
  in_params.set<vos_type>("Fields",{"field_1", "field_2", "field_3", "field_packed"});
  return in_params;
}
/*===================================================================================================================*/
} // undefined namespace
