#include "mct_coupling/ScreamContext.hpp"

#include "dynamics/register_dynamics.hpp"

#include "share/io/register_io.hpp"

#include "physics/register_physics.hpp"
#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/shoc/scream_shoc_interface.hpp"

// #include "control/tests/dummy_grid.hpp"
#include "control/atmosphere_driver.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/se_grid.hpp"
#include "share/scream_types.hpp"
#include "share/scream_session.hpp"
#include "scream_config.h"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

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

  // First of all, disable all fpes we may have enabled.
  // Store the mask, so we can restore before returning.
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();

  // First of all, initialize the scream session
  scream::initialize_scream_session();

  // Create the context
  auto& c = ScreamContext::singleton();

  // Create the C MPI_Comm from the Fortran one
  MPI_Comm mpi_comm_c = MPI_Comm_f2c(f_comm);
  auto& atm_comm = c.create<ekat::Comm>(mpi_comm_c);


  // Create a parameter list for inputs
  printf("[scream] reading parameters from yaml file: %s\n",input_yaml_file);
  ekat::ParameterList scream_params("Scream Parameters");
  parse_yaml_file (input_yaml_file, scream_params);
  scream_params.print();

  ekat::error::runtime_check(scream_params.isSublist("Atmosphere Driver"),
       "Error! Sublist 'Atmosphere Driver' not found inside '" +
       std::string(input_yaml_file) + "'.\n");

  auto& ad_params = scream_params.sublist("Atmosphere Driver");

  // Need to register products in the factories *before* we attempt to create any.
  // In particular, register all atm processes, and all grids managers.
  if ( scream_params.isSublist("Output YAML Files") ) 
  {
    auto& io_file_params = scream_params.sublist("Output YAML Files");
    register_io( io_file_params );
  }
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
  ekat::disable_all_fpes();
  ekat::enable_fpes(fpe_mask);

}
/*===============================================================================================*/
void scream_run (const Real& dt) {
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  using namespace scream;
  using namespace scream::control;

  // First of all, enable only scream fpes.
  // Store the mask, so we can restore before returning.
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  ekat::enable_default_fpes();

  // Get the context
  auto& c = ScreamContext::singleton();

  // Get the AD, and run it
  auto& ad = c.getNonConst<AtmosphereDriver>();
  ad.run(dt);

  (void) dt;

  // Restore the FPE flag as it was when control was handed to us.
  ekat::disable_all_fpes();
  ekat::enable_fpes(fpe_mask);
}
/*===============================================================================================*/
void scream_finalize (/* args ? */) {
  using namespace scream;
  using namespace scream::control;

  // First of all, enable only scream fpes.
  // Store the mask, so we can restore before returning.
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  ekat::enable_default_fpes();

  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  // Get the context
  auto& c = ScreamContext::singleton();

  // Get the AD, and finalize it
  auto& ad = c.getNonConst<AtmosphereDriver>();
  ad.finalize();

  // Clean up also P3 stuff
  p3::P3GlobalForFortran::deinit();

  // Restore the FPE flag as it was when control was handed to us.
  ekat::disable_all_fpes();
  ekat::enable_fpes(fpe_mask);
}

} // extern "C"
