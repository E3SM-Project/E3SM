#include "ekat/scream_session.hpp"
#include "scream_config.h"

#include "ekat/scream_parse_yaml_file.hpp"
#include "ekat/util/ekat_md_array.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/scream_pack.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/se_grid.hpp"
#include "control/atmosphere_driver.hpp"

#include "dynamics/register_dynamics.hpp"
#include "physics/register_physics.hpp"

#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/shoc/scream_shoc_interface.hpp"

#include "control/tests/dummy_grid.hpp"

#include "mct_coupling/ScreamContext.hpp"
#include "mct_coupling/scream_scorpio_interface.hpp"
#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_types.hpp"

using scream::Real;

extern "C"
{

/*===============================================================================================*/
// WARNING: make sure input_yaml_file is a null-terminated string!
void scream_init (const MPI_Fint& f_comm,
                  const int& start_ymd,
                  const int& start_tod,
                  const char* input_yaml_file,
                  const int& compid) {
  using namespace scream;
  using namespace scream::control;
  using namespace scream::scorpio;

  // First of all, disable all fpes we may have enabled.
  // Store the mask, so we can restore before returning.
  int fpe_mask = get_enabled_fpes();
  disable_all_fpes();

  // First of all, initialize the scream session
  scream::initialize_scream_session();

  // Create the context
  auto& c = ScreamContext::singleton();

  // Create the C MPI_Comm from the Fortran one
  MPI_Comm mpi_comm_c = MPI_Comm_f2c(f_comm);
  auto& atm_comm = c.create<scream::Comm>(mpi_comm_c);


  // Create a parameter list for inputs
  printf("[scream] reading parameterr from yaml file: %s\n",input_yaml_file);
  ParameterList scream_params("Scream Parameters");
  parse_yaml_file (input_yaml_file, scream_params);
  scream_params.print();

  error::runtime_check(scream_params.isSublist("Atmosphere Driver"),
       "Error! Sublist 'Atmosphere Driver' not found inside '" +
       std::string(input_yaml_file) + "'.\n");

  auto& ad_params = scream_params.sublist("Atmosphere Driver");

  // Need to register products in the factories *before* we attempt to create any.
  // In particular, register all atm processes, and all grids managers.
  register_dynamics();
  register_physics();

  // Create the bare ad, then init it
  auto& ad = c.create<AtmosphereDriver>();

  // Recall that e3sm uses the int YYYYMMDD to store a date
  std::cout << "start_ymd: " << start_ymd << "\n";
  const int dd = start_ymd % 100;
  const int mm = (start_ymd / 100) % 100;
  const int yy = start_ymd / 10000;
  util::TimeStamp time (yy,mm,dd,start_tod);

  // Init and run (to finalize, wait till checks are completed,
  // or you'll clear the field repo!)
  ad.initialize(atm_comm,ad_params,time);

  // Restore the FPE flag as it was when control was handed to us.
  disable_all_fpes();
  enable_fpes(fpe_mask);

  //                 SCORPIO                        //
  // Create the set of SCORPIO output files and their respective
  // dimensions and variables.
  eam_init_pio_subsystem(f_comm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  // Register the set of output files:
  register_outfile("example_pio_structured_v2.nc");
  // Register the set of dimensions per output file
  register_dimension("example_pio_structured_v2.nc","x","horizontal distance",10);
  register_dimension("example_pio_structured_v2.nc","y","vertical distance",3);
  register_dimension("example_pio_structured_v2.nc","z","height",2);
  register_dimension("example_pio_structured_v2.nc","time","time",0);
  // Register the set of variables per output file
  std::string vec_time[] = {"time"};
  std::string vec_x[]    = {"x"};
  std::string vec_y[]    = {"y"};
  std::string vec_z[]    = {"z"};
  std::string vec_xt[]   = {"x","time"}; 
  std::string vec_yt[]   = {"y","time"}; 
  std::string vec_xyt[]  = {"x","y","time"}; 
  std::string vec_xyzt[]  = {"x","y","z","time"}; 
  register_variable("example_pio_structured_v2.nc","time","time",1,vec_time, PIO_REAL,"t");
  register_variable("example_pio_structured_v2.nc","x","answer to space and time",1,vec_x, PIO_REAL,"x-real");
  register_variable("example_pio_structured_v2.nc","y","answer to space and time",1,vec_y, PIO_REAL,"y-real");
  register_variable("example_pio_structured_v2.nc","z","answer to space and time",1,vec_z, PIO_REAL,"z-real");
  register_variable("example_pio_structured_v2.nc","bar","answer to space and time",3,vec_xyt, PIO_REAL,"xyt-real");
  register_variable("example_pio_structured_v2.nc","foo","answer to space and time",3,vec_xyt, PIO_REAL,"xyt-real");
  register_variable("example_pio_structured_v2.nc","bar_flip","answer to space and time",3,vec_xyt, PIO_REAL,"yxt-real");
  register_variable("example_pio_structured_v2.nc","foo_flip","answer to space and time",3,vec_xyt, PIO_REAL,"yxt-real");
  register_variable("example_pio_structured_v2.nc","foo_big","answer to space and time",4,vec_xyzt, PIO_REAL,"xyzt-real");

  eam_pio_enddef();

  std::array<Real,10> x_data;
  std::array<Real, 3> y_data;
  std::array<Real, 2> z_data = {100,200};
  Int m_rank; 
  MPI_Comm_rank(mpi_comm_c,&m_rank);
  for (int ii=0;ii<10;ii++) {
    x_data[ii] = ii + 1.0;
  }
  for (int jj=0;jj<3;jj++) {
    y_data[jj] = jj + 1.0;
  }
  grid_write_data_array("example_pio_structured_v2.nc","x",10,x_data.data());
  grid_write_data_array("example_pio_structured_v2.nc","y",3,y_data.data());
  grid_write_data_array("example_pio_structured_v2.nc","z",2,z_data.data());
  
  //                 SCORPIO DONE                    //
}
/*===============================================================================================*/
void scream_run (const Real& dt) {
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  using namespace scream;
  using namespace scream::control;
  using namespace scream::scorpio;
  using ekat::util::data;

  // First of all, enable only scream fpes.
  // Store the mask, so we can restore before returning.
  int fpe_mask = get_enabled_fpes();
  disable_all_fpes();
  enable_default_fpes();

  // Get the context
  auto& c = ScreamContext::singleton();

  // Get the AD, and run it
  auto& ad = c.getNonConst<AtmosphereDriver>();
  ad.run(dt);

  // TEST
  ekat::util::md_array<Real,10>     foo_data_1d;
  ekat::util::md_array<Real,10,3>   foo_data_2d;
  ekat::util::md_array<Real,2,3,10> foo_data_3d;
  ekat::util::md_array<Real,10,3>   foo_data_2d_inv;
  for (int ii=0;ii<10;ii++) {
      foo_data_1d[ii] = ii + 1.0;
    for (int jj=0;jj<3;jj++) {
      foo_data_2d[jj][ii] = 10.0*(jj+1.0) + ii + 1.0;
      foo_data_2d_inv[ii][jj] = 10.0*(jj+1.0) + ii + 1.0;
      for (int kk=0;kk<2;kk++) {
        foo_data_3d[kk][jj][ii] = 100.0*(kk+1.0) + 10.0*(jj+1.0) + (ii+1.0);
      }
    }
  }
  std::array<Int,2> dimlen = {3,10};
  std::array<Int,3> dimlen_3d = {3,10,2};
  pio_update_time("example_pio_structured_v2.nc",dt);
  grid_write_data_array("example_pio_structured_v2.nc","foo",dimlen,ekat::util::data(foo_data_2d));
  grid_write_data_array("example_pio_structured_v2.nc","bar",dimlen,ekat::util::data(foo_data_2d_inv));
  grid_write_data_array("example_pio_structured_v2.nc","bar_flip",dimlen,ekat::util::data(foo_data_2d_inv));
  grid_write_data_array("example_pio_structured_v2.nc","foo_flip",dimlen,ekat::util::data(foo_data_2d));
  grid_write_data_array("example_pio_structured_v2.nc","foo_big",dimlen_3d,ekat::util::data(foo_data_3d));
  sync_outfile("example_pio_structured_v2.nc"); 
  // TEST END

  (void) dt;

  // Restore the FPE flag as it was when control was handed to us.
  disable_all_fpes();
  enable_fpes(fpe_mask);
}
/*===============================================================================================*/
void scream_finalize (/* args ? */) {
  using namespace scream;
  using namespace scream::control;

  // First of all, enable only scream fpes.
  // Store the mask, so we can restore before returning.
  int fpe_mask = get_enabled_fpes();
  disable_all_fpes();
  enable_default_fpes();

  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  // Get the context
  auto& c = ScreamContext::singleton();

  // Get the AD, and finalize it
  auto& ad = c.getNonConst<AtmosphereDriver>();
  ad.finalize();

  // Clean up also P3 stuff
  p3::P3GlobalForFortran::deinit();

  // Restore the FPE flag as it was when control was handed to us.
  disable_all_fpes();
  enable_fpes(fpe_mask);
}

} // extern "C"
