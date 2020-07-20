#include <catch2/catch.hpp>

#include "scream_config.h"
#include "ekat/scream_pack.hpp"
#include "mct_coupling/scream_scorpio_interface.hpp"
#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/ekat_md_array.hpp"



namespace {

Real f_x(const Real x, const Real t); 
Real f_y(const Real y, const Real t); 
Real f_z(const Real z, const Real t); 
Int ind_x(const Int ii);
Int ind_y(const Int jj);
Int ind_z(const Int kk);
Int ind_t(const Int tt);

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
  const char* vec_time[] = {"time"};
  const char* vec_x[]    = {"x"};
  const char* vec_y[]    = {"y"};
  const char* vec_z[]    = {"z"};
  const char* vec_xt[]   = {"x","time"};
  const char* vec_xyt[]  = {"x","y","time"};
  const char* vec_xyzt[] = {"x","y","z","time"};
 
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
  eam_pio_enddef(outfilename);
  // Create data to be written
  std::array<Real,10> x_data;
  std::array<Real, 5> y_data;
  std::array<Real, 2> z_data;
  std::array<Int,1> xdim = {10};
  std::array<Int,1> ydim = {5};
  std::array<Int,1> zdim = {2};
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

  for (decltype(x_data)::size_type ii=0;ii<x_data.size();ii++) {
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
    for (decltype(x_data)::size_type ii=0;ii<x_data.size();ii++) {
      test_data_1d[ii]  = f_x(x_data[ii],tt*dt);
      test_index_1d[ii] = ind_x(ii) + ind_t(tt);
      for (int jj=0;jj<5;jj++) {
        test_data_2d[jj][ii]  = f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt);
        test_index_2d[jj][ii] = ind_y(jj) + ind_x(ii) + ind_t(tt);
        for (int kk=0;kk<2;kk++) {
          test_data_3d[kk][jj][ii]  = f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt) + f_z(z_data[kk],tt*dt);
          test_index_3d[kk][jj][ii] = ind_z(kk) + ind_y(jj) + ind_x(ii) + ind_t(tt);
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
/* ================================================================================================================ */
TEST_CASE("scorpio_interface_input", "") {

  using namespace scream;
  using namespace scream::scorpio;
  using ekat::util::data;

  int compid=0;
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  // Register the set of output files:
  std::string infilename = "scorpio_output_test.nc";
  register_infile(infilename);

  const char* vec_time[] = {"time"};
  const char* vec_x[]    = {"x"};
  const char* vec_y[]    = {"y"};
  const char* vec_z[]    = {"z"};
  const char* vec_xy[]   = {"x","y"}; 
  const char* vec_xyz[]  = {"x","y","z"};
 
  register_variable(infilename,"time","time",1,vec_time, PIO_REAL,"t");
  register_variable(infilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  register_variable(infilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  register_variable(infilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  register_variable(infilename,"data_1d","test value for 1d field",1,vec_x, PIO_REAL,"x-real");
  register_variable(infilename,"data_2d","test value for 2d field",2,vec_xy, PIO_REAL,"xy-real");
  register_variable(infilename,"data_3d","test value for 3d field",3,vec_xyz, PIO_REAL,"xyz-real");
  register_variable(infilename,"index_1d","test value for 1d field",1,vec_x, PIO_INT,"x-int");
  register_variable(infilename,"index_2d","test value for 2d field",2,vec_xy, PIO_INT,"xy-int");
  register_variable(infilename,"index_3d","test value for 3d field",3,vec_xyz, PIO_INT,"xyz-int");

  // Create data to be written
  std::array<Real,10> x_data;
  std::array<Real, 5> y_data;
  std::array<Real, 2> z_data;
  std::array<Int,1> xdim = {10};
  std::array<Int,1> ydim = {5};
  std::array<Int,1> zdim = {2};
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

  grid_read_data_array(infilename,"x",xdim,ekat::util::data(x_data));
  grid_read_data_array(infilename,"y",ydim,ekat::util::data(y_data));
  grid_read_data_array(infilename,"z",zdim,ekat::util::data(z_data));

  for (decltype(x_data)::size_type ii=0;ii<x_data.size();ii++) {
    REQUIRE( x_data[ii] == 2.0*pi/x_data.size()*(ii+1) );
  }
  for (decltype(y_data)::size_type jj=0;jj<5;jj++) {
    REQUIRE( y_data[jj] == 4.0*pi/y_data.size()*(jj+1) );
  }
  for (decltype(z_data)::size_type kk=0;kk<2;kk++) {
    REQUIRE( z_data[kk] == 100*(kk+1) );
  }

  Real dt = 1.0;
  for (int tt=0;tt<3;tt++) {
    pio_update_time(infilename,-999.0);
    grid_read_data_array(infilename,"data_1d",dimlen_1d,ekat::util::data(test_data_1d));
    grid_read_data_array(infilename,"data_2d",dimlen_2d,ekat::util::data(test_data_2d));
    grid_read_data_array(infilename,"data_3d",dimlen_3d,ekat::util::data(test_data_3d));
    grid_read_data_array(infilename,"index_1d",dimlen_1d,ekat::util::data(test_index_1d));
    grid_read_data_array(infilename,"index_2d",dimlen_2d,ekat::util::data(test_index_2d));
    grid_read_data_array(infilename,"index_3d",dimlen_3d,ekat::util::data(test_index_3d));
    for (int ii=0;ii<x_data.size();ii++) {
      REQUIRE(test_data_1d[ii] == f_x(x_data[ii],tt*dt));
      REQUIRE(test_index_1d[ii]== ind_x(ii) + ind_t(tt));
      for (int jj=0;jj<5;jj++) {
        REQUIRE(test_index_2d[jj][ii] == ind_y(jj) + ind_x(ii) + ind_t(tt));
        REQUIRE(test_data_2d[jj][ii] == (f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt)));
        for (int kk=0;kk<2;kk++) {
          REQUIRE(test_index_3d[kk][jj][ii] == ind_z(kk) + ind_y(jj) + ind_x(ii) + ind_t(tt));
          REQUIRE(test_data_3d[kk][jj][ii] ==((f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt))+f_z(z_data[kk],tt*dt)));
        } //kk
      } //jj
    } //ii
  } //tt
  eam_pio_finalize();
} // TEST scorpio_interface_input
/* ================================================================================================================ */
/*                                   Local functions to be used for tests:                                          */
Real f_x(const Real x, const Real t) {
  Real f;
  f = 0.1 * cos(x+t);
  return f;
}
Real f_y(const Real y, const Real t) {
  Real f;
  f = sin(y+t);
  return f;
}
Real f_z(const Real z, const Real t) {
  Real f;
  f = z;
  return f;
}
Int ind_x(const Int ii) {
  return ii;
}
Int ind_y(const Int jj) {
  return jj*100;
}
Int ind_z(const Int kk) {
  return 1000*kk;
}
Int ind_t(const Int tt) {
  return 10000*tt;
}
/* ================================================================================================================ */
} //namespace
