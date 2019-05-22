#include "atmosphere_process_group.hpp"

#include "share/util/string_utils.hpp"

namespace scream {

AtmosphereProcessGroup::
AtmosphereProcessGroup (const Comm& comm, const ParameterList& params)
 : m_comm(comm)
{
  // Get number of processes in the group and the scheduling type (Sequential vs Parallel)
  m_group_size = params.get<int>("Number of Entries");

  // Create the individual atmosphere processes
  for (int i=0; i<m_group_size; ++i) {
    // The comm to be passed to the processes construction is
    //  - the same as the input comm if num_entries=1 or sched_type=Sequential
    //  - a sub-comm of the input comm otherwise
    Comm proc_comm = m_comm;
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
      error::runtime_abort("Error! Parallel schedule type not yet implemented.\n");
    }

    const auto& params_i = params.sublist(util::strint("Process",i));
    const std::string& process_name = params_i.get<std::string>("Process Name");
    m_atm_processes.emplace_back(AtmosphereProcessFactory::instance().create(process_name,proc_comm,params_i));

    // NOTE: the shared_ptr of the new atmosphere process *MUST* have been created correctly. Namely, the creation
    //       process must have set up enable_shared_from_this's status correctly. This is done by the library-provided
    //       templated function 'create_atm_process<T>. However, if the user has decided to roll his/her own creator
    //       function to be registered in the AtmosphereProcessFactory, he/she may have forgot to set the self pointer
    //       in the process. To make sure this is not the case, we check that the weak_ptr in the newly created
    //       atmosphere process (which comes through inheritance from enable_shared_from_this) is valid.
    scream_require_msg(!m_atm_processes.back()->weak_from_this().expired(),
                       "Error! The newly created std::shared_ptr<AtmosphereProcess> does not correctly setup the 'enable_shared_from_this' interface.\n"
                       "       Did you by change register your own creator function in the AtmosphereProccessFactory class?\n"
                       "       If so, don't. Instead, use the instantiation of create_atmosphere_process<T>, with T = YourAtmProcessClassName.\n");

    // Update the grid types of the group, given the needs of the newly created process
    for (const auto& name : m_atm_processes.back()->get_required_grids()) {
      m_required_grids.insert(name);
    }
  }

  if (params.get<std::string>("Schedule Type") == "Sequential") {
    m_group_schedule_type = GroupScheduleType::Sequential;
  } else if (params.get<std::string>("Schedule Type") == "Parallel") {
    m_group_schedule_type = GroupScheduleType::Parallel;
    error::runtime_abort("Error! Parallel schedule not yet implemented.\n");
  } else {
    error::runtime_abort("Error! Invalid 'Schedule Type'. Available choices are 'Parallel' and 'Sequential'.\n");
  }

}

void AtmosphereProcessGroup::set_grid (const std::shared_ptr<const GridsManager> grids_manager) {

  for (int i=0; i<m_group_size; ++i) {
    m_atm_processes[i]->set_grid(grids_manager);

    // Add inputs/outputs to the list of inputs/outputs of this group
    for (const auto& id : m_atm_processes[i]->get_required_fields()) {
      m_required_fields.insert(id);
    }
    for (const auto& id : m_atm_processes[i]->get_computed_fields()) {
      m_computed_fields.insert(id);
    }
  }
}

void AtmosphereProcessGroup::initialize (const util::TimeStamp& t0) {
  // Now that we have the comm for the processes in the group, we can initialize them
  for (int i=0; i<m_group_size; ++i) {
    m_atm_processes[i]->initialize(t0);
  }
}

void AtmosphereProcessGroup::run (const double dt) {
  for (auto atm_proc : m_atm_processes) {
    atm_proc->run(dt);
  }
}

void AtmosphereProcessGroup::finalize (/* what inputs? */) {
  for (auto atm_proc : m_atm_processes) {
    atm_proc->finalize(/* what inputs? */);
  }
}

void AtmosphereProcessGroup::register_fields (FieldRepository<Real, device_type>& field_repo) const {
  for (const auto& atm_proc : m_atm_processes) {
    atm_proc->register_fields(field_repo);

    // Make sure processes are not calling methods they shouldn't on the repo
    scream_require_msg(field_repo.repository_state()==RepoState::Open,
                         "Error! Atmosphere processes are *not* allowed to modify the state of the repository.\n");

    // Check that the required fields are indeed in the repo now
    for (const auto& id : atm_proc->get_required_fields()) {
      scream_require_msg(field_repo.has_field(id), "Error! Process '" + atm_proc->name() + "' failed to register required field '" + id.get_identifier() + "'.\n");
    }
    for (const auto& id : atm_proc->get_computed_fields()) {
      scream_require_msg(field_repo.has_field(id), "Error! Process '" + atm_proc->name() + "' failed to register computed field '" + id.get_identifier() + "'.\n");
    }
  }
}

void AtmosphereProcessGroup::set_required_field_impl (const Field<const Real, device_type>& f) {
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->requires(f.get_header().get_identifier())) {
      atm_proc->set_required_field(f);
    }
  }
}

void AtmosphereProcessGroup::set_computed_field_impl (const Field<Real, device_type>& f) {
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->computes(f.get_header().get_identifier())) {
      atm_proc->set_computed_field(f);
    }
  }
}

} // namespace scream
