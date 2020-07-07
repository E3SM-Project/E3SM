#include <catch2/catch.hpp>

#include "scream_config.h"
#include "ekat/scream_pack.hpp"
#include "mct_coupling/scream_scorpio_interface.hpp"
#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/ekat_md_array.hpp"


namespace {

TEST_CASE("scorpio_interface_output", "") {

  using namespace scream;
  using namespace scream::scorpio;
  using ekat::util::data;

  // Create the set of SCORPIO output files and their respective
  // dimensions and variables.
  int compid=0;
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  // Register the set of output files:
  std::string outfilename = "scorpio_output_test.nc";
  register_outfile(outfilename);
  // Register the set of dimensions per output file
  int xlen=10, ylen=5, zlen=2;
  register_dimension(outfilename,"x","horizontal distance",xlen);
  register_dimension(outfilename,"y","vertical distance",ylen);
  register_dimension(outfilename,"z","height",zlen);
  register_dimension(outfilename,"time","time",0);
  // Register the set of variables per output file
  std::string vec_time[] = {"time"};
  std::string vec_x[]    = {"x"};
  std::string vec_y[]    = {"y"};
  std::string vec_z[]    = {"z"};
  std::string vec_xt[]   = {"x","time"}; 
  std::string vec_xyt[]  = {"x","y","time"}; 
  std::string vec_xyzt[] = {"x","y","z","time"};
 
  register_variable(outfilename,"time","time",1,vec_time, PIO_REAL,"t");
  register_variable(outfilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  register_variable(outfilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  register_variable(outfilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  register_variable(outfilename,"data_1d","test value for 1d field",2,vec_xt, PIO_REAL,"xt-real");
  register_variable(outfilename,"data_2d","test value for 2d field",3,vec_xyt, PIO_REAL,"xyt-real");
  register_variable(outfilename,"data_3d","test value for 3d field",4,vec_xyzt, PIO_REAL,"xyzt-real");
  register_variable(outfilename,"index_1d","test value for 1d field",2,vec_xt, PIO_INT,"xt-int");
  register_variable(outfilename,"index_2d","test value for 2d field",3,vec_xyt, PIO_INT,"xyt-int");
  register_variable(outfilename,"index_3d","test value for 3d field",4,vec_xyzt, PIO_INT,"xyzt-int");
  // Finished with the initialization of variables in output file
  eam_pio_enddef();
  // Create data to be written
  std::array<Real,10> x_data;
  std::array<Real, 5> y_data;
  std::array<Real, 2> z_data;
  std::array<Int,1> xdim = {10};
  std::array<Int,1> ydim = {5};
  std::array<Int,1> zdim = {3};
  std::array<Int,1> dimlen_1d = {10};
  std::array<Int,2> dimlen_2d = {5,10};
  std::array<Int,3> dimlen_3d = {2,5,10};
  ekat::util::md_array<Real,10>       test_data_1d;
  ekat::util::md_array<Real, 5,10>    test_data_2d;
  ekat::util::md_array<Real, 2, 5,10> test_data_3d;
  ekat::util::md_array<Int,10>        test_index_1d;
  ekat::util::md_array<Int, 5,10>     test_index_2d;
  ekat::util::md_array<Int, 2, 5,10>  test_index_3d;
  Real pi = 2*acos(0.0);

  for (int ii=0;ii<x_data.size();ii++) {
    x_data[ii] = 2.0*pi/x_data.size()*(ii+1);
  }
  for (int jj=0;jj<5;jj++) {
    y_data[jj] = 4.0*pi/y_data.size()*(jj+1);
  }
  for (int kk=0;kk<2;kk++) {
    z_data[kk] = 100*(kk+1);
  }
  // Write dimension data and initial fields
  grid_write_data_array(outfilename,"x",xdim,ekat::util::data(x_data));
  grid_write_data_array(outfilename,"y",ydim,ekat::util::data(y_data));
  grid_write_data_array(outfilename,"z",zdim,ekat::util::data(z_data));
  sync_outfile(outfilename); 
  // write multiple timesteps of data for comparison:
  Real dt = 1.0;
  for (int tt=0;tt<3;tt++) {
    for (int ii=0;ii<x_data.size();ii++) {
      test_data_1d[ii]  = 0.1 * cos(x_data[ii] + tt*dt); // phase shift by dt
      test_index_1d[ii] = 10000*tt + ii;
      for (int jj=0;jj<5;jj++) {
        test_data_2d[jj][ii]  = test_data_1d[ii] * sin(y_data[jj] + tt*dt); //phase shift by dt
        test_index_2d[jj][ii] = test_index_1d[ii] + 100*jj;
        for (int kk=0;kk<2;kk++) {
          test_data_3d[kk][jj][ii]  = test_data_2d[jj][ii] + z_data[kk];
          test_index_3d[kk][jj][ii] = test_index_2d[jj][ii] + 1000*kk;
        } //kk
      } //jj
    } //ii
    pio_update_time(outfilename,tt*dt);
    grid_write_data_array(outfilename,"index_1d",dimlen_1d,ekat::util::data(test_index_1d));
    grid_write_data_array(outfilename,"index_2d",dimlen_2d,ekat::util::data(test_index_2d));
    grid_write_data_array(outfilename,"index_3d",dimlen_3d,ekat::util::data(test_index_3d));
    grid_write_data_array(outfilename,"data_1d",dimlen_1d,ekat::util::data(test_data_1d));
    grid_write_data_array(outfilename,"data_2d",dimlen_2d,ekat::util::data(test_data_2d));
    grid_write_data_array(outfilename,"data_3d",dimlen_3d,ekat::util::data(test_data_3d));
    sync_outfile(outfilename); 
  } //tt

  eam_pio_finalize();
} // TEST scorpio_interface_output

} //namespace
