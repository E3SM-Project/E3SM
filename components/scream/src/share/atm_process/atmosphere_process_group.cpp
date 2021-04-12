#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_string_utils.hpp"

namespace scream {

AtmosphereProcessGroup::
AtmosphereProcessGroup (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_comm(comm)
{
  // Get number of processes in the group and the scheduling type (Sequential vs Parallel)
  m_group_size = params.get<int>("Number of Entries");
  EKAT_REQUIRE_MSG (m_group_size>0, "Error! Invalid group size.\n");

  if (m_group_size>1) {
    if (params.get<std::string>("Schedule Type") == "Sequential") {
      m_group_schedule_type = ScheduleType::Sequential;
    } else if (params.get<std::string>("Schedule Type") == "Parallel") {
      m_group_schedule_type = ScheduleType::Parallel;
      ekat::error::runtime_abort("Error! Parallel schedule not yet implemented.\n");
    } else {
      ekat::error::runtime_abort("Error! Invalid 'Schedule Type'. Available choices are 'Parallel' and 'Sequential'.\n");
    }
  } else {
    // Pointless to handle this group as parallel, if only one process is in it
    m_group_schedule_type = ScheduleType::Sequential;
  }

  // Create the individual atmosphere processes
  m_group_name = "Group [";
  m_group_name += m_group_schedule_type==ScheduleType::Sequential
                ? "Sequential]:" : "Parallel]:";
  for (int i=0; i<m_group_size; ++i) {
    // The comm to be passed to the processes construction is
    //  - the same as the input comm if num_entries=1 or sched_type=Sequential
    //  - a sub-comm of the input comm otherwise
    ekat::Comm proc_comm = m_comm;
    if (m_group_schedule_type==ScheduleType::Parallel) {
      // This is what's going to happen:
      //  - the processes in the group are going to be run in parallel
      //  - each rank is assigned ONE atm process
      //  - all the atm processes not assigned to this rank are filled with
      //    an instance of RemoteProcessStub (to be implemented), which is a do-nothing class,
      //    only responsible to keep track of dependencies
      //  - the input parameter list should specify for each atm process the number
      //    of mpi ranks dedicated to it. Obviously, these numbers should add up
      //    to the size of the input communicator.
      //  - this class is then responsible of 'combining' the results togehter,
      //    including remapping input/output fields to/from the sub-comm
      //    distribution.
      ekat::error::runtime_abort("Error! Parallel schedule type not yet implemented.\n");
    }

    const auto& params_i = params.sublist(ekat::strint("Process",i));
    const std::string& process_name = params_i.get<std::string>("Process Name");
    m_atm_processes.emplace_back(AtmosphereProcessFactory::instance().create(process_name,proc_comm,params_i));

    // NOTE: the shared_ptr of the new atmosphere process *MUST* have been created correctly.
    //       Namely, the creation process must have set up enable_shared_from_this's status correctly.
    //       This is done by the library-provided templated function 'create_atm_process<T>.
    //       However, if the user has decided to roll his/her own creator function to be registered
    //       in the AtmosphereProcessFactory, he/she may have forgot to set the self pointer in the process.
    //       To make sure this is not the case, we check that the weak_ptr in the newly created
    //       atmosphere process (which comes through inheritance from enable_shared_from_this) is valid.
    EKAT_REQUIRE_MSG(!m_atm_processes.back()->weak_from_this().expired(),
        "Error! The newly created std::shared_ptr<AtmosphereProcess> did not correctly setup\n"
        "       the 'enable_shared_from_this' interface.\n"
        "       Did you by chance register your own creator function in the AtmosphereProccessFactory class?\n"
        "       If so, don't. Instead, use the instantiation of create_atmosphere_process<T>,\n"
        "       with T = YourAtmProcessClassName.\n");

    // Update the grid types of the group, given the needs of the newly created process
    for (const auto& name : m_atm_processes.back()->get_required_grids()) {
      m_required_grids.insert(name);
    }

    m_group_name += " ";
    m_group_name += m_atm_processes.back()->name();
  }
}

void AtmosphereProcessGroup::set_grids (const std::shared_ptr<const GridsManager> grids_manager) {

  // The atm process group (APG) simply 'concatenates' required/computed
  // fields of the stored process. There is a single exception to this
  // rule, in case of sequential splitting: if an atm proc requires a
  // field that is computed by a previous atm proc in the group, that
  // field is not exposed as a required field of the group.

  for (auto& atm_proc : m_atm_processes) {
    atm_proc->set_grids(grids_manager);

    // Add inputs/outputs to the list of inputs of the group
    for (const auto& req : atm_proc->get_required_fields()) {
      process_required_field(req);
    }
    for (const auto& req : atm_proc->get_computed_fields()) {
      add_computed_field(req);
    }
    for (const auto& req : atm_proc->get_required_groups()) {
      add_required_group(req);
    }
    for (const auto& req : atm_proc->get_updated_groups()) {
      add_updated_group(req);
    }
  }
}

void AtmosphereProcessGroup::initialize_impl (const TimeStamp& t0) {
  for (auto& atm_proc : m_atm_processes) {
    atm_proc->initialize(t0);
  }
}

void AtmosphereProcessGroup::run_impl (const Real dt) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    run_sequential(dt);
  } else {
    run_parallel(dt);
  }
}

void AtmosphereProcessGroup::run_sequential (const Real dt) {
  // Get the timestamp at the beginning of the step and advance it.
  auto ts = timestamp();
  ts += dt;

  for (auto atm_proc : m_atm_processes) {
    // Run the process
    atm_proc->run(dt);
  }
}

void AtmosphereProcessGroup::run_parallel (const Real /* dt */) {
  EKAT_REQUIRE_MSG (false,"Error! Parallel splitting not yet implemented.\n");
}

void AtmosphereProcessGroup::finalize_impl (/* what inputs? */) {
  for (auto atm_proc : m_atm_processes) {
    atm_proc->finalize(/* what inputs? */);
  }
}

void AtmosphereProcessGroup::
set_required_group (const FieldGroup<const Real>& group)
{
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    auto atm_proc = m_atm_processes[iproc];

    if (atm_proc->requires_group(group.m_info->m_group_name,group.grid_name())) {
      atm_proc->set_required_group(group);
    }
  }
}

void AtmosphereProcessGroup::
set_updated_group (const FieldGroup<Real>& group)
{
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    auto atm_proc = m_atm_processes[iproc];

    if (atm_proc->updates_group(group.m_info->m_group_name,group.grid_name())) {
      atm_proc->set_updated_group(group);
    }
  }
}

void AtmosphereProcessGroup::
register_fields (const std::map<std::string,std::shared_ptr<FieldManager<Real>>>& field_mgrs) const {
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    const auto& atm_proc = m_atm_processes[iproc];
    atm_proc->register_fields(field_mgrs);

#ifdef SCREAM_DEBUG
    // Make sure processes are not calling methods they shouldn't on the repo
    for (const auto& it : field_mgrs) {
      EKAT_REQUIRE_MSG(it.second->repository_state()==RepoState::Open,
          "Error! Atmosphere processes are *not* allowed to modify the state of the repository.\n");

      // Check that the required fields are indeed in the repo now
      for (const auto& id : atm_proc->get_required_fields()) {
        EKAT_REQUIRE_MSG(it.second->has_field(id),
            "Error! Process '" + atm_proc->name() + "' failed to register required field '" + id.get_id_string() + "'.\n");
      }
      for (const auto& id : atm_proc->get_computed_fields()) {
        EKAT_REQUIRE_MSG(it.second->has_field(id),
            "Error! Process '" + atm_proc->name() + "' failed to register computed field '" + id.get_id_string() + "'.\n");
      }
    }
#endif
  }
}
void AtmosphereProcessGroup::set_required_field_impl (const Field<const Real>& f) {
  const auto& fid = f.get_header().get_identifier();
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->requires_field(fid)) {
      atm_proc->set_required_field(f);
    }
  }
}

void AtmosphereProcessGroup::set_computed_field_impl (const Field<Real>& f) {
  const auto& fid = f.get_header().get_identifier();
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->computes_field(fid)) {
      atm_proc->set_computed_field(f);
    }
    // In sequential scheduling, some fields may be computed by
    // a process and used by the next one. In this case, the field
    // does not figure as 'input' for the group, but we still
    // need to set it in the processes that need it.
    if (atm_proc->requires_field(fid)) {
      atm_proc->set_required_field(f);
    }
  }
}

void AtmosphereProcessGroup::
process_required_group (const GroupRequest& req) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    if (updates_group(req.name,req.grid)) {
      // Some previous atm proc updated this group, so it's not an 'input'
      // of the atm group as a whole. However, we might need a different
      // pack size. So, instead of adding to the required groups,
      // we add to the computed fields. This way we don't modify the inputs
      // of the group, and still manage to communicate to the AD the pack size
      // that we need.
      // NOTE; we don't have a way to check if all the fields in the group
      //       are computed by previous processes, since we don't have
      //       the list of all fields in this group.
      add_updated_group(req);
    } else {
      add_required_group(req);
    }
  } else {
    // In parallel schedule, the inputs of all processes are inputs of the group
    add_required_group(req);
  }
}

void AtmosphereProcessGroup::
process_required_field (const FieldIdentifier& fid) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    if (computes_field(fid)) {
      // Some previous atm proc computes this field, so it's not an 'input'
      // of the group as a whole. However, we might need a different pack size,
      // or want to add it to a different group. So, instead of adding to
      // the required fields, we add to the computed fields. This way we
      // don't modify the inputs of the group, and still manage to communicate
      // to the AD the pack size and group affiliations that we need
      add_computed_field(fid);
    } else {
      add_required_field(fid);
    }
  } else {
    // In parallel schedule, the inputs of all processes are inputs of the group
    add_required_field(fid);
  }
}

} // namespace scream
