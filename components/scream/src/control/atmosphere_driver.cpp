#include "atmosphere_driver.hpp"

#include <share/atmosphere_process.hpp>
#include <share/error_defs.hpp>
#include <share/util/string_utils.hpp>

#include "atm_processes_scheduling.hpp"
#include "surface_coupling.hpp"
#include "dynamics/atmosphere_dynamics.hpp"

namespace scream {
namespace control {

void AtmosphereDriver::initialize ( /* inputs? */ ) {

  // Create processes, store them in the requested order
  create_atm_processes();

  // Initialize the processes
  for (auto& atm_process : m_atm_processes) {
    atm_process->initialize( /* inputs ? */ );
  }

  // Initialize the FR for device fields
  // (i.e., fields that are used by processes during computations)
}

void AtmosphereDriver::run ( /* inputs ? */ ) {
  // Here we have to call the 'run' method on each process that is stored on this rank.
  for (auto& atm_process : m_atm_processes) {
    // Run the current process
    atm_process->run();
  }
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  // Finalize all processes
  for (auto& atm_process : m_atm_processes) {
    atm_process->finalize( /* inputs? */ );
  }
}


void AtmosphereDriver::create_atm_processes ()  {
  // Parse some sort of data, to create the schedule of all the atmosphere
  // processes, including order and parallel/sequential distribution

  ProcessesSchedule schedule;
  // Parse some sort of data to fill schedule

  // Check: make sure that in the global schedule there is one (and only one) dynamics process,
  //        as well as one (and only one) surface coupling process. Also, check for no repetitions.
  error::runtime_check(schedule.check_process_present("Dynamics"),"No Dynamics process found. This makes no sense. Aborting",-1);
  error::runtime_check(schedule.check_process_present("Surface_Coupling"),"No SurfaceCoupling process found. This makes no sense. Aborting",-1);
  schedule.check_unique_processes();

  // Assign to each rank a pipeline of processes.
  std::vector<std::pair<std::string,Comm>> my_pipeline = schedule.get_processes_pipeline(m_atm_comm);

  for (auto& name : my_pipeline) {
    auto name_upper = util::upper_case(name.first);
    if (name_upper=="DYNAMICS") {
      m_atm_processes.push_back(std::make_shared<AtmosphereDynamics>());
    } else if (name_upper=="SURFACE_COUPLING") {
      m_atm_processes.push_back(std::make_shared<SurfaceCoupling>());
    } else {
      std::string err_msg = "Error! Uknown/unsupported atmosphere process '" + name.first + "'.\n";
      error::runtime_abort(err_msg,-1);
    }
  }
}

// ==================== STUB for testing ================= //

int driver_stub() {
  // Test that we can at least create an AD...
  AtmosphereDriver ad;
  (void) ad;

  // Return the answer to life, universe, everything...
  return 42;
}

}  // namespace control
}  // namespace scream
