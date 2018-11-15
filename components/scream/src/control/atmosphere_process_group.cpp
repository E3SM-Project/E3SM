#include "atmosphere_process_group.hpp"

#include <share/util/string_utils.hpp>

namespace scream {

AtmosphereProcessGroup::AtmosphereProcessGroup (const ParameterList& params) {
  // Get number of processes in the group and the scheduling type (Sequential vs Parallel)
  m_group_size = params.get<int>("Number of Entries");

  // Create the individual atmosphere processes
  for (int i=0; i<m_group_size; ++i) {
    const auto& params_i = params.sublist(util::strint("Process",i));
    const std::string& process_name = params_i.get<std::string>("Process Name");
    m_atm_processes.emplace_back(AtmosphereProcessFactory::instance().create(process_name,params_i));
  }

  m_group_schedule_type = params.get<GroupScheduleType>("Schedule Type");

  error::runtime_check(m_group_schedule_type==GroupScheduleType::Sequential, "Error! Parallel schedule not yet implemented.\n");
}

void AtmosphereProcessGroup::initialize (const Comm& comm) {
  m_comm = comm;

  // Now that we have the comm for the processes in the group, we can initialize them
  for (int i=0; i<m_group_size; ++i) {
    // The comm to be passed to the processes construction is
    //  - the same as the input comm if num_entries=1 or sched_type=Sequential
    //  - a sub-comm of the input comm otherwise
    Comm proc_comm = comm;
    if (m_group_size>1 && m_group_schedule_type==GroupScheduleType::Parallel) {
      // This is what's going to happen:
      //  - the processes in the group are going to be run in parallel
      //  - each rank is assigned ONE atm process
      //  - all the atm processes not assigned to this rank are filled with
      //    an instance of RemoteProcessStub (to be implemented), which is a do-nothing class,
      //    only responsible to keep track of dependencies
      //  - the input parameter list should specify for each atm process the number
      //    of mpi ranks dedicated to it. Obviously, these numbers should add up
      //    to the size of the input communicator.
      // We therefore construct entries_comm following the specs of the processs
      // we have to handle on this rank,  and pass the same comm to all the processes.

      // In order to figure out what atm process this rank should handle,
      // we need to know how many ranks are assigned to each atm process
    }
  
    m_atm_processes[i]->initialize(proc_comm);

    // Add inputs/outputs to the list of inputs/outputs of this group
    for (const auto& id : m_atm_processes[i]->get_required_fields()) {
      m_required_fields.insert(id);
    }
    for (const auto& id : m_atm_processes[i]->get_computed_fields()) {
      m_computed_fields.insert(id);
    }
  }
}

void AtmosphereProcessGroup::run        (/* what inputs? */) {
  for (auto atm_proc : m_atm_processes) {
    atm_proc->run(/* what inputs? */);
  }
}

void AtmosphereProcessGroup::finalize   (/* what inputs? */) {
  for (auto atm_proc : m_atm_processes) {
    atm_proc->finalize(/* what inputs? */);
  }
}

void AtmosphereProcessGroup::set_required_field_impl (const Field<const Real*, ExecMemSpace, MemoryManaged>& f) {
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->requires(f.get_header().get_identifier())) {
      atm_proc->set_required_field(f);
    }
  }
}

void AtmosphereProcessGroup::set_computed_field_impl (const Field<      Real*, ExecMemSpace, MemoryManaged>& f) {
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->computes(f.get_header().get_identifier())) {
      atm_proc->set_computed_field(f);
    }
  }
}

// The AtmosphereProcessGroup creator function (for the factory)
AtmosphereProcess* create_process_group (const ParameterList& params) {
  return new AtmosphereProcessGroup(params);
}

} // namespace scream
