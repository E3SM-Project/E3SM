#include "mo_garand_atmos_io.h"
#include "simple_netcdf.hpp"
#include <netcdf.h>

using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

// Read in the data, then use only the first column, and copy it to all of the model columns
// In the end, all model columns will be identical
void read_atmos(std::string input_file, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs, int ncol) {
  simple_netcdf::SimpleNetCDF io;
  io.open(input_file , NC_NOWRITE);

  int nlay = io.getDimSize("lay");
  int nlev = io.getDimSize("lev");

  //p_lay = real2d("p_lay",ncol,nlay);
  //t_lay = real2d("t_lay",ncol,nlay);
  //p_lev = real2d("p_lev",ncol,nlev);
  //t_lev = real2d("t_lev",ncol,nlev);

  real2d tmp2d;
  // p_lay
  io.read(tmp2d,"p_lay");
  // for (int ilay=1 ; ilay <= nlay ; ilay++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
    p_lay(icol,ilay) = tmp2d(1,ilay);
  });
  // t_lay
  io.read(tmp2d,"t_lay");
  // for (int ilay=1 ; ilay <= nlay ; ilay++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA ( int ilay, int icol) {
    t_lay(icol,ilay) = tmp2d(1,ilay);
  });
  // p_lev
  tmp2d = real2d();  // Reset tmp2d to avoid warnings about reallocating during file read
  io.read(tmp2d,"p_lev");
  // for (int ilev=1 ; ilev <= nlev ; ilev++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( SimpleBounds<2>(nlev,ncol) , YAKL_LAMBDA ( int ilev, int icol) {
    p_lev(icol,ilev) = tmp2d(1,ilev);
  });
  // t_lev
  io.read(tmp2d,"t_lev");
  // for (int ilev=1 ; ilev <= nlev ; ilev++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( SimpleBounds<2>(nlev,ncol) , YAKL_LAMBDA( int ilev, int icol) {
    t_lev(icol,ilev) = tmp2d(1,ilev);
  });

  int ngas = 8;
  string1d gas_names("gas_names",ngas);
  gas_names(1) = std::string("h2o");
  gas_names(2) = std::string("co2");
  gas_names(3) = std::string("o3" );
  gas_names(4) = std::string("n2o");
  gas_names(5) = std::string("co" );
  gas_names(6) = std::string("ch4");
  gas_names(7) = std::string("o2" );
  gas_names(8) = std::string("n2" );

  // Initialize GasConcs object with an "ncol" given from the calling program
  gas_concs.init(gas_names,ncol,nlay);

  tmp2d = real2d();     // Reset the tmp2d variable
  for (int igas=1 ; igas <= ngas ; igas++) {
    std::string vmr_name = "vmr_"+gas_names(igas);
    if ( ! io.varExists(vmr_name) ) { stoprun("ERROR: gas does not exist in input file"); }
    // Read in 2-D varaible
    io.read(tmp2d,vmr_name);
    // Create 1-D variable with just the first column
    real1d tmp1d("tmp1d",nlay);
    // for (int i=1 ; i <= nlay ; i++) {
    parallel_for( SimpleBounds<1>(nlay) , YAKL_LAMBDA (int i) {
      tmp1d(i) = tmp2d(1,i);
    });
    // Call set_vmr with only the first column from the data file copied among all of the model columns
    gas_concs.set_vmr( gas_names(igas) , tmp1d );
  }

  io.close();
}



void write_sw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, real2d const &flux_dir, int ncol) {
  simple_netcdf::SimpleNetCDF io;
  io.open(fileName , NC_WRITE);
  io.write(flux_up  , "sw_flux_up_result"  , {"col_new","lev"});
  io.write(flux_dn  , "sw_flux_dn_result"  , {"col_new","lev"});
  io.write(flux_dir , "sw_flux_dir_result" , {"col_new","lev"});
  io.close();
}



void write_lw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, int ncol) {
  simple_netcdf::SimpleNetCDF io;
  io.open(fileName , NC_WRITE);
  io.write(flux_up , "lw_flux_up_result" , {"col_new","lev"});
  io.write(flux_dn , "lw_flux_dn_result" , {"col_new","lev"});
  io.close();
}

