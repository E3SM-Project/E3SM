#include "atmosphere_driver.hpp"

#include <share/atmosphere_process.hpp>
#include <share/error_defs.hpp>

namespace scream {
namespace control {

void AtmosphereDriver::initialize ( /* inputs? */ ) {
  // Initialize the FieldRepository (FR) for host fields
  // (i.e., fields that are I/O w.r.t the coupler)

  // Create processes, store them in the requested order

  // Check: make sure there is a dynamics process.
  bool dynamics_found = false;
  for (const auto& atm_subcomp : m_atm_processes) {
    if (atm_subcomp->type()==AtmosphereProcessType::Dynamics) {
      // We already found a Dynamics process. Issue an error.
      errors::runtime_abort("More than one Dynamics process found. This makes no sense. Aborting.\n", -1); 
      dynamics_found = true;
    }
  }
  if (!dynamics_found) {
    // We did not find a Dynamics process. Issue an error.
    errors::runtime_abort("No Dynamics process found. This makes no sense. Aborting.\n", -2); 
  }

  // Initialize the processes
  for (auto& atm_subcomp : m_atm_processes) {
    atm_subcomp->initialize( /* inputs ? */ );
  }

  // Initialize the FR for device fields
  // (i.e., fields that are used by processes during computations)

  // Other? Perhaps setup a registry of fields that need copying
  // host->device at the beginning of run and device->host at the end of run?
}

void AtmosphereDriver::run ( /* inputs ? */ ) {
  // We have some flexibility over the order in which some things happen here.
  // For instance, do we want to copy (host->device) input fields from coupler
  // right away or should we do it 'just in time'? This may be irrelevant now
  // (and therefore do it right away), but it *may*  be relevant if we enable
  // parallel processes execution (in which case, we may delegate each of these copies
  // to (one of) the processes that use it).

  // Copy (host->device) data from coupler views

  // Run processes.
  // Question: are we going to enforce same data structure for all processes?
  //           If yes, then this step is easy, otherwise we need
  //           to sync data between each each run call.
  // Note: if we want to allow parallel processes, we need to
  //       consider it here. Perhaps a big if, like
  //
  // if (parallel_parametrizations) {
  //   dispatch processes depending on mpi rank
  // } else {
  //   plow through processes list in the order in which they are stored.
  // }

  // Copy (device->host) data to coupler views
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  // Finalize all processes
  for (auto& atm_subcomp : m_atm_processes) {
    atm_subcomp->finalize( /* inputs? */ );
  }
}

int driver_stub() {
  // Test that we can at least create an AD...
  AtmosphereDriver ad;
  (void) ad;

  // Return the answer to life, universe, everything...
  return 42;
}

}  // namespace control
}  // namespace scream
