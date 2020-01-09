#include "interface/ScreamContext.hpp"
#include "control/atmosphere_driver.hpp"
#include "share/mpi/scream_comm.hpp"
#include "share/parameter_list.hpp"
#include "share/scream_session.hpp"

extern "C"
{

void scream_init (const MPI_Fint& f_comm, const int& start_ymd, const int& start_tod) {

  // First of all, initialize the scream session
  scream::initialize_scream_session();

  // // Get the context
  // auto& c = scream::ScreamContext::singleton();

  // // Create the C MPI_Comm from the Fortran one
  // MPI_Comm mpi_comm_c = MPI_Comm_f2c(f_comm);
  // auto& comm = c.create<scream::Comm>(mpi_comm_c);

  // // Fill the params, somehow
  // scream::ParameterList ad_params;

  // // Init the initial time, somehow
  // // You probably need more than yymmdd to init the time stamp
  // // E.g: you have to formalize switches for leap/nonleap years
  // scream::util::TimeStamp t0;

  // // Create the bare ad, then init it
  // // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  // auto& ad = c.create<scream::control::AtmosphereDriver>();
  // ad.initialize(comm, ad_params, t0);

  (void) start_ymd;
  (void) start_tod;
  (void) f_comm;
}

void scream_run (const double& dt) {
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.

  // // Get the context
  // auto& c = scream::ScreamContext::singleton();

  // // Get the AD, and run it
  // auto& ad = c.getNonConst<scream::control::AtmosphereDriver>();
  // ad.run(dt);
  (void) dt;
}

void scream_finalize (/* args ? */) {
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  // // Get the context
  // auto& c = scream::ScreamContext::singleton();

  // // Get the AD, and finalize it
  // auto& ad = c.getNonConst<scream::control::AtmosphereDriver>();
  // ad.finalize();
}

} // extern "C"
