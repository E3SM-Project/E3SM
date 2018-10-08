#include "atmosphere_driver.hpp"

#include <share/atmosphere_process.hpp>
#include <share/error_defs.hpp>
#include <share/util/string_utils.hpp>

#include "control/atm_processes_scheduling.hpp"
#include "control/surface_coupling.hpp"
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
  for (auto& atm_process : m_atm_processes) {
    const auto& inputs  = atm_process->get_required_fields();
    const auto& outputs = atm_process->get_required_fields();

    for (auto& header : inputs) {
      m_device_field_repo.register_field(header);
    }
    for (auto& header : outputs) {
      m_device_field_repo.register_field(header);
    }
  }

  // Make sure all fields are allocated (with the correct dimensions)
  /*  How/where to we allocate the fields? I see two scenarios:
        1) fields are not (all) yet dimensioned. Two sub-options:
             a) the field repo dimensions all the feilds during registration_complete
             b) the driver calls get_field(header)->set_dimensions(...) for all fields
           Option a requires the repo to know all the extents (e.g., number of elements,
           columns, tracers,...). This information is probably in the driver already, so
           we *may* just do b...
        2) fields are all dimensioned at register time. This allows the registration to
           not only create the Field's, but also to create the underlying Kokkos::View's.
           This can be an advantage for atm_process-specific dimensions, that the
           driver does not need to know about, but can potentially create issues if
           for some reason (read: bug) different processes think there are a different
           number of element/columns/tracers/...
      Probably the best solution is a compromise: let each process initialize the dimensions
      that are specific to it (e.g., number of chemical species for microphysics), but leave
      it to the driver (or field_repo) to initialize dimensions for tags that are of 'general
      purpose', such as Element, GaussPoint, etc.
      NOTE: the compromise approach requires to change a bit the FieldHeader implementation.
            In particular, as of now, FieldHeader does not allow multiple calls to set_dimensions,
            to avoid changing dimensions when the views are already allocated. But we could
            probably easily relax this restriction, while keepin the dimensions in check; we
            probably would need to add a bool flag, stating whether all dimensions have been
            set, using -1 (or similar) as placeholder for dimensions yet to be set.
   */
  m_device_field_repo.registration_complete();

  // Set all the fields in the processes needing them (before, they only had headers)
  for (auto& atm_process : m_atm_processes) {
    const auto& inputs  = atm_process->get_required_fields();
    const auto& outputs = atm_process->get_required_fields();

    for (auto& header : inputs) {
      auto g = m_device_field_repo.get_field(header);
      auto g2 = g.get_const();

      atm_process->set_required_field(m_device_field_repo.get_field(header).get_const());
    }
    for (auto& header : outputs) {
      atm_process->set_computed_field(m_device_field_repo.get_field(header));
    }
  }

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
