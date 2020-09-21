#include <catch2/catch.hpp>

#include "scream_config.h"
#include "share/scream_types.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scorpio.hpp"

#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/simple_grid.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"

#include "ekat/ekat_pack.hpp"

namespace {
/* ======= Internal Functions needed for test ======*/
using namespace scream;
using gid_type       = long;          // TODO: use int64_t? int? template class on gid_type?
using device_type    = DefaultDevice; // TODO: template class on device type
using kokkos_types   = KokkosTypes<device_type>;
using dofs_list_type = kokkos_types::view_1d<gid_type>;
void set_spatial_vectors(Real* x, Real* y, Real* z, std::vector<int> dims3d, const dofs_list_type gids, const Int total_cols);
void update_index_output(Kokkos::View<Real*> index_1d, Kokkos::View<Real**> index_2d, Kokkos::View<Real***> index_3d, std::vector<int> dims3d, const Int tt, const dofs_list_type gids, bool& init);
void update_data_output(Kokkos::View<Real*> data_1d, Kokkos::View<Real**> data_2d, Kokkos::View<Real***> data_3d, Real* x, Real* y, Real* z, const std::vector<int> dims3d, const Real t, const dofs_list_type gids);
 
/* ================================================================================================================ */
TEST_CASE("scorpio_yaml", "") {

  using namespace scream;
  using namespace scream::scorpio;
  using namespace ekat::units;
  using Device = DefaultDevice;

  Int num_cols = 10;
  Int num_levs = 5;
  Int num_comp = 2;

  int compid=0;  // For CIME based builds this will be the integer ID assigned to the atm by the component coupler.  For testing we simply set to 0
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("User Provided",create_user_provided_grids_manager);
  // Set the dummy grid in the UserProvidedGridManager
  // Recall that this class stores *static* members, so whatever
  // we set here, will be reflected in the GM built by the factory.
  UserProvidedGridsManager upgm;
  auto dummy_grid= std::make_shared<SimpleGrid>("Physics",num_cols,num_levs,io_comm);
  upgm.set_grid(dummy_grid);
  auto gids       = upgm.get_grid("Physics")->get_dofs_gids();
  auto phys_lt    = upgm.get_grid("Physics")->get_native_dof_layout();

/* The first step is to establish a Field Manager Repo to work with.  This example is fashioned from the 'field_repo'
 * test from /share/tests/field_tests.hpp                                                                          */ 
  std::vector<FieldTag> tag1d  = {FieldTag::Column};
  std::vector<FieldTag> tag1db = {FieldTag::VerticalLevel};
  std::vector<FieldTag> tag1dc = {FieldTag::Component};
  std::vector<FieldTag> tag2d  = {FieldTag::Column, FieldTag::VerticalLevel};
  std::vector<FieldTag> tag3d  = {FieldTag::Column, FieldTag::VerticalLevel, FieldTag::Component};

  FieldIdentifier fid_x("x", tag1d,  m);
  FieldIdentifier fid_y("y", tag1db, m);
  FieldIdentifier fid_z("z", tag1dc, m);
  FieldIdentifier fid_index_1d("index_1d", tag1d, kg/s);
  FieldIdentifier fid_index_2d("index_2d", tag2d, kg/s);
  FieldIdentifier fid_index_3d("index_3d", tag3d, kg/s);
  FieldIdentifier fid_data_1d("data_1d", tag1d, m/s);
  FieldIdentifier fid_data_2d("data_2d", tag2d, m/s);
  FieldIdentifier fid_data_3d("data_3d", tag3d, m/s);

  std::vector<int> dimsx = {phys_lt.dim(0)};
  std::vector<int> dimsy = {num_levs};
  std::vector<int> dimsz = {num_comp};
  std::vector<int> dims1d = {dimsx[0]};
  std::vector<int> dims2d = {dimsx[0], dimsy[0]};
  std::vector<int> dims3d = {dimsx[0], dimsy[0], dimsz[0]};

  fid_x.set_dimensions(dimsx);
  fid_y.set_dimensions(dimsy);
  fid_z.set_dimensions(dimsz);
  fid_index_1d.set_dimensions(dims1d);
  fid_index_2d.set_dimensions(dims2d);
  fid_index_3d.set_dimensions(dims3d);
  fid_data_1d.set_dimensions(dims1d);
  fid_data_2d.set_dimensions(dims2d);
  fid_data_3d.set_dimensions(dims3d);

  fid_x.set_grid_name("Physics");
  fid_y.set_grid_name("Physics");
  fid_z.set_grid_name("Physics");
  fid_index_1d.set_grid_name("Physics");
  fid_index_2d.set_grid_name("Physics");
  fid_index_3d.set_grid_name("Physics");
  fid_data_1d.set_grid_name("Physics");
  fid_data_2d.set_grid_name("Physics");
  fid_data_3d.set_grid_name("Physics");

  // Field Repo    /* TODO: Expand test to support packing */
  FieldRepository<Real,DefaultDevice>  repo;
  repo.registration_begins();
  repo.register_field(fid_x,"group_1");
  repo.register_field(fid_y,"group_1");
  repo.register_field(fid_z,"group_1");
  // All fields in the field manager are consider REAL, so index still used but converted to REAL
  repo.register_field(fid_index_1d,"group_2");
  repo.register_field(fid_index_2d,"group_2");
  repo.register_field(fid_index_3d,"group_2");
  repo.register_field(fid_data_1d,"group_3");
  repo.register_field(fid_data_2d,"group_3");
  repo.register_field(fid_data_3d,"group_3");
  repo.registration_ends();

  // Initialize spatial vectors
  auto xd = repo.get_field(fid_x).get_view(); 
  auto yd = repo.get_field(fid_y).get_view();
  auto zd = repo.get_field(fid_z).get_view();
  auto xh = Kokkos::create_mirror_view( xd );
  auto yh = Kokkos::create_mirror_view( yd );
  auto zh = Kokkos::create_mirror_view( zd );
  set_spatial_vectors(xh.data(), yh.data(), zh.data(), dims3d, gids, num_cols);
  Kokkos::deep_copy(xd,xh);
  Kokkos::deep_copy(yd,yh);
  Kokkos::deep_copy(zd,zh);

/* The next step is to initialize the set of IO class objects to handle output,
 * and sync it with the field manager we just created. */
  // Initialize the PIO subsystem:
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  // Gather control data for the test from the scorpio_control.yaml file.  Note this is similar to the scream_control.yaml file used in regular runs.
  char *input_yaml_file = "scorpio_control.yaml";
  printf("[scream] reading parameters from yaml file: %s\n",input_yaml_file);
  ekat::ParameterList scorpio_params("Scorpio Parameters");
  parse_yaml_file (input_yaml_file, scorpio_params);
  // Retrieve the list of scorpio output control files from scream_params.
  auto& output_control_files = scorpio_params.get<std::vector<std::string>>("Output YAML Files");    // First grab the list of Output files from the control YAML

  // For each output control file create a separate instance of the io class.
  std::vector<AtmosphereOutput> AtmOutput_Instances;
  for (int ii=0;ii<output_control_files.size();ii++) {
    ekat::ParameterList out_params(output_control_files[ii]);
    parse_yaml_file(output_control_files[ii],out_params);
    AtmosphereOutput output_instance(io_comm,out_params);
    AtmOutput_Instances.push_back(output_instance);
  }
  // Initialize the output files:
  for (auto& out_ins : AtmOutput_Instances)
  {
    out_ins.init(repo,upgm);
  }

  // Run and record data:
  bool index_init = true;
  Real dt = 1.0;
  Int max_tstep = 10;
  for (int tstep = 0;tstep<max_tstep;++tstep)
  {
    Real l_time = tstep*dt;
    // Update index variables
    auto index_1d_dev = repo.get_field(fid_index_1d).get_view(); 
    auto index_2d_dev = repo.get_field(fid_index_2d).get_reshaped_view<Real**>();
    auto index_3d_dev = repo.get_field(fid_index_3d).get_reshaped_view<Real***>();
    auto index_1d_h = Kokkos::create_mirror_view( index_1d_dev );
    auto index_2d_h = Kokkos::create_mirror_view( index_2d_dev );
    auto index_3d_h = Kokkos::create_mirror_view( index_3d_dev );
    update_index_output(index_1d_dev, index_2d_dev, index_3d_dev, dims3d, (int) dt, gids, index_init);
    Kokkos::deep_copy(index_1d_dev,index_1d_h);
    Kokkos::deep_copy(index_2d_dev,index_2d_h);
    Kokkos::deep_copy(index_3d_dev,index_3d_h);
    // Update data variables
    auto data_1d_dev = repo.get_field(fid_data_1d).get_view(); 
    auto data_2d_dev = repo.get_field(fid_data_2d).get_reshaped_view<Real**>();
    auto data_3d_dev = repo.get_field(fid_data_3d).get_reshaped_view<Real***>();
    auto data_1d_h = Kokkos::create_mirror_view( data_1d_dev );
    auto data_2d_h = Kokkos::create_mirror_view( data_2d_dev );
    auto data_3d_h = Kokkos::create_mirror_view( data_3d_dev );
    update_data_output(data_1d_dev, data_2d_dev, data_3d_dev, xh.data(), yh.data(), zh.data(), dims3d, l_time, gids);
    Kokkos::deep_copy(data_1d_dev,data_1d_h);
    Kokkos::deep_copy(data_2d_dev,data_2d_h);
    Kokkos::deep_copy(data_3d_dev,data_3d_h);

    for (auto& out_ins : AtmOutput_Instances)
    {
      out_ins.run(repo,upgm,l_time);
    }
  }

  //Finished with PIO, finalize the system
  for (auto& out_ins : AtmOutput_Instances)
  {
    out_ins.finalize();
    out_ins.check_status();
  }
  eam_pio_finalize();
  upgm.clean_up();
} // TEST_CASE scorpio_yaml


/* ================================================================================================================ */
/* ====================                           INTERNAL FUNCTIONS                           ==================== */
/* ================================================================================================================ */


void set_spatial_vectors(Real* x, Real* y, Real* z, std::vector<int> dims3d, const dofs_list_type gids, const Int total_cols) 
{
  Real pi = 2*acos(0.0);

  for (int ii=0;ii<dims3d[0];ii++)
  {
    Int g_ii = gids(ii);
    x[ii] = 2.0*pi/total_cols*(g_ii+1);
  }
  for (int jj=0;jj<dims3d[1];jj++) 
  {
    y[jj] = 4.0*pi/dims3d[1]*(jj+1);
  }
  for (int kk=0;kk<dims3d[2];kk++)
  {
    z[kk] = 100*(kk+1);
  }
}

void update_index_output(Kokkos::View<Real*> index_1d, Kokkos::View<Real**> index_2d, Kokkos::View<Real***> index_3d, const std::vector<int> dims3d, const Int dt, const dofs_list_type gids,bool& init)
{
  if (init)
  {
    for (int ii=0;ii<dims3d[0];ii++)
    {
      Int g_ii = gids(ii);
      index_1d(ii) = g_ii*1.;
      for (int jj=0;jj<dims3d[1];jj++) 
      {
        index_2d(ii,jj) = g_ii*1. + jj*100.;
        for (int kk=0;kk<dims3d[2];kk++)
        {
          index_3d(ii,jj,kk) = g_ii*1. + jj*100. + kk*1000.;
        }
      }
    }
    init = false;
  } else {
    for (int ii=0;ii<dims3d[0];ii++)
    {
      index_1d(ii) += 10000.*dt;
      for (int jj=0;jj<dims3d[1];jj++) 
      {
        index_2d(ii,jj) += 10000.*dt;
        for (int kk=0;kk<dims3d[2];kk++)
        {
          index_3d(ii,jj,kk) += 10000.*dt;
        }
      }
    }
  }
}

void update_data_output(Kokkos::View<Real*> data_1d, Kokkos::View<Real**> data_2d, Kokkos::View<Real***> data_3d, Real* x, Real* y, Real* z, const std::vector<int> dims3d, const Real t, const dofs_list_type gids)
{
  for (int ii=0;ii<dims3d[0];ii++)
  {
    Int g_ii = gids(ii);
    data_1d(ii) = 0.1 * cos(x[g_ii]+t);
    for (int jj=0;jj<dims3d[1];jj++) 
    {
      data_2d(ii,jj) = sin(y[jj]+t) * data_1d(ii);
      for (int kk=0;kk<dims3d[2];kk++)
      {
        data_3d(ii,jj,kk) = z[kk] + data_2d(ii,jj);
      }
    }
  }
}
/* ========== Unit Test Ideas to Implement ========== */
/*
 * 1) Same fields but written on physics for one file and dynamics for another: Tests output grid.
 * 2) Same run but output in instant, averaged, min, max: Tests output type.
 * 3) Set max steps per file to 4 for a 10 step run, make sure steps 1-4, 4-8 and 9-10 are on three separate files: Tests max steps/file
 *    -xtra test, do 1 snap per file averaged over 3 steps.  Tests average works with multiple files.
 * 4) Write restart, i.e. use DAG to determine which fields are required by a process and write as output.  Should be done with a full AD test.
 * 5) Have # PIO tasks < Total number of tasks.  Thus any one PIO task is responsible for writing the info from >1 MPI tasks.  Tests the ability
 *    to use PIO_STRIDE or to distribute PIO tasks over nodes.
 */    
} // namespace
