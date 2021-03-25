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

    // NOTE: the shared_ptr of the new atmosphere process *MUST* have been created correctly. Namely, the creation
    //       process must have set up enable_shared_from_this's status correctly. This is done by the library-provided
    //       templated function 'create_atm_process<T>. However, if the user has decided to roll his/her own creator
    //       function to be registered in the AtmosphereProcessFactory, he/she may have forgot to set the self pointer
    //       in the process. To make sure this is not the case, we check that the weak_ptr in the newly created
    //       atmosphere process (which comes through inheritance from enable_shared_from_this) is valid.
    EKAT_REQUIRE_MSG(!m_atm_processes.back()->weak_from_this().expired(),
                       "Error! The newly created std::shared_ptr<AtmosphereProcess> does not correctly setup the 'enable_shared_from_this' interface.\n"
                       "       Did you by chance register your own creator function in the AtmosphereProccessFactory class?\n"
                       "       If so, don't. Instead, use the instantiation of create_atmosphere_process<T>, with T = YourAtmProcessClassName.\n");

    // Update the grid types of the group, given the needs of the newly created process
    for (const auto& name : m_atm_processes.back()->get_required_grids()) {
      m_required_grids.insert(name);
    }

    m_group_name += " ";
    m_group_name += m_atm_processes.back()->name();
  }

#ifdef SCREAM_DEBUG
  m_field_repo     = nullptr;
  m_bkp_field_repo = nullptr;
#endif
}

void AtmosphereProcessGroup::set_grids (const std::shared_ptr<const GridsManager> grids_manager) {

  // The atm process group (APG) acts a bit different than other atm processes.
  // There is a catch when it comes to exposing a list of required/computed
  // fields: if (and only if) the splitting is sequential, and contains
  // 2+ processes, it is possible that some inputs of atmosphere processes
  // in the list are computed by processes that precede them.
  // In this case, such fields are not really 'requirements' of this
  // group, but simply "outputs" of the group (which happen to also
  // be reused by other atm processes in the group).

  // The (bare) names of the computed fields. We use this to figure out
  // if an 'input' to a following process should be marked as input to
  // the APG as a whole, depending on whether one of the previous atm
  // processes in the group is computing the field.
  // NOTE: this is only in the case of sequential scheduling
  std::set<std::string> computed;

  std::pair<decltype(m_required_fields)::iterator,bool>  it_bool;
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    auto atm_proc = m_atm_processes[iproc];
    atm_proc->set_grids(grids_manager);

    // Add inputs to the list of inputs of the group
    for (const auto& fid : atm_proc->get_required_fields()) {
      process_required_field(fid);
    }

    // Add outputs to the list of outputs of the group
    for (const auto& fid : atm_proc->get_computed_fields()) {
      m_computed_fields.insert(fid);
    }
  }
}

void AtmosphereProcessGroup::initialize_impl (const TimeStamp& t0) {
  // Now that we have the comm for the processes in the group, we can initialize them
  for (int i=0; i<m_group_size; ++i) {
    m_atm_processes[i]->initialize(t0);
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

  for (int iproc=0; iproc<m_group_size; ++iproc) {
    // Run the process
    auto atm_proc = m_atm_processes[iproc];
    atm_proc->run(dt);

#ifdef SCREAM_DEBUG
    if (m_field_repo!=nullptr && m_bkp_field_repo!=nullptr) {
      // Check that this process did not update any field it should not have updated
      const auto& computed = atm_proc->get_computed_fields();
      for (const auto& it : (*m_field_repo)) {
        for (const auto& f : it.second) {
          const auto& fid = f.first;
          const auto& gn = fid.get_grid_name();
          auto& f_old = m_bkp_field_repo->get_field(fid);
          bool field_is_unchanged = views_are_equal(f_old,f.second);
          if (ekat::contains(computed,fid)) {
            // For fields that changed, make sure the time stamp has been updated
            const auto& ts = f.second.get_header().get_tracking().get_time_stamp();
            EKAT_REQUIRE_MSG(field_is_unchanged || ts==timestamp(),
                               "Error! Process '" + atm_proc->name() + "' updated field '" +
                                fid.get_id_string() + "', but it did not update its time stamp.\n");
          } else {
            // For non computed fields, make sure the field in m_repo matches
            // the copy in m_bkp_repo, and update the copy in m_bkp_repo.
            // There are three ok scenarios: 1) field is unchanged, 2) field is computed,
            // 3) field is a remapped input.
            EKAT_REQUIRE_MSG(field_is_unchanged,
                               "Error! Process '" + atm_proc->name() + "' updated field '" +
                                fid.get_id_string() + "', which it wasn't allowed to update.\n");

          }
        }
      }
      // Update the computed fields in the bkp field repo
      for (auto& fid : computed) {
        auto src = m_field_repo->get_field(fid).get_view();
        auto dst = m_bkp_field_repo->get_field(fid).get_view();

        Kokkos::deep_copy(dst,src);
      }
    }
#endif
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
  const auto& name = group.m_info->m_group_name;
  const auto& grid = group.m_grid_name;

  for (int iproc=0; iproc<m_group_size; ++iproc) {
    auto atm_proc = m_atm_processes[iproc];

    if (atm_proc->requires_group(name,grid)) {
      atm_proc->set_required_group(group);
      // Some groups might be optional, so don't error out if they're empty
      if (not group.m_info->empty()) {
        // If the group is 'bundled', we remap the bundled group,
        // otherwise remap each individual field.
        if (group.m_info->m_bundled) {
          const auto& f = *group.m_bundle;
          const auto& fid = f.get_header().get_identifier();
          process_required_field(fid);
        } else {
          // We also might need to add each field of the group to the remap in/out
          for (const auto& it : group.m_fields) {
            const auto& f = *it.second;
            const auto& fid = f.get_header().get_identifier();
            process_required_field(fid);
          }
        }
      }
    }
  }
}

void AtmosphereProcessGroup::
set_updated_group (const FieldGroup<Real>& group)
{
  const auto& name = group.m_info->m_group_name;
  const auto& grid = group.m_grid_name;

  for (int iproc=0; iproc<m_group_size; ++iproc) {
    auto atm_proc = m_atm_processes[iproc];

    if (atm_proc->updates_group(name,grid)) {
      atm_proc->set_updated_group(group);
      // Some groups might be optional, so don't error out if they're empty
      if (not group.m_info->empty()) {
        // If the group is 'bundled', we remap the bundled group,
        // otherwise remap each individual field.
        if (group.m_info->m_bundled) {
          const auto& f = *group.m_bundle;
          const auto& fid = f.get_header().get_identifier();
          m_computed_fields.insert(fid);
        } else {
          // We also might need to add each field of the group to the remap in/out
          for (const auto& it : group.m_fields) {
            const auto& f = *it.second;
            const auto& fid = f.get_header().get_identifier();
            m_computed_fields.insert(fid);
          }
        }
      }
    }
  }
}

void AtmosphereProcessGroup::register_fields (FieldRepository<Real>& field_repo) const {
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    const auto& atm_proc = m_atm_processes[iproc];
    atm_proc->register_fields(field_repo);

#ifdef SCREAM_DEBUG
    // Make sure processes are not calling methods they shouldn't on the repo
    EKAT_REQUIRE_MSG(field_repo.repository_state()==RepoState::Open,
                         "Error! Atmosphere processes are *not* allowed to modify the state of the repository.\n");

    // Check that the required fields are indeed in the repo now
    for (const auto& id : atm_proc->get_required_fields()) {
      EKAT_REQUIRE_MSG(field_repo.has_field(id), "Error! Process '" + atm_proc->name() + "' failed to register required field '" + id.get_id_string() + "'.\n");
    }
    for (const auto& id : atm_proc->get_computed_fields()) {
      EKAT_REQUIRE_MSG(field_repo.has_field(id), "Error! Process '" + atm_proc->name() + "' failed to register computed field '" + id.get_id_string() + "'.\n");
    }
#endif
  }
}

std::set<AtmosphereProcessGroup::GroupRequest>
AtmosphereProcessGroup::get_required_groups () const {
  std::set<AtmosphereProcessGroup::GroupRequest> groups;
  for (auto atm_proc : m_atm_processes) {
    const auto& groups_i = atm_proc->get_required_groups();
    groups.insert(groups_i.begin(),groups_i.end());
  }
  return groups;
}

std::set<AtmosphereProcessGroup::GroupRequest>
AtmosphereProcessGroup::get_updated_groups () const {
  std::set<AtmosphereProcessGroup::GroupRequest> groups;
  for (auto atm_proc : m_atm_processes) {
    const auto& groups_i = atm_proc->get_updated_groups();
    groups.insert(groups_i.begin(),groups_i.end());
  }
  return groups;
}

void AtmosphereProcessGroup::set_required_field_impl (const Field<const Real>& f) {
  const auto& fid = f.get_header().get_identifier();
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->requires(fid)) {
      atm_proc->set_required_field(f);
    }
  }
}

void AtmosphereProcessGroup::set_computed_field_impl (const Field<Real>& f) {
  const auto& fid = f.get_header().get_identifier();
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->computes(fid)) {
      atm_proc->set_computed_field(f);
    }
    // In sequential scheduling, some fields may be computed by
    // a process and used by the next one. In this case, the field
    // does not figure as 'input' for the group, but we still
    // need to set it in the processes that need it.
    if (atm_proc->requires(fid)) {
      atm_proc->set_required_field(f);
    }
  }
}

#ifdef SCREAM_DEBUG
void AtmosphereProcessGroup::
set_field_repos (const FieldRepository<Real>& repo,
                 const FieldRepository<Real>& bkp_repo) {
  m_field_repo = &repo;
  m_bkp_field_repo = &bkp_repo;
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->type()==AtmosphereProcessType::Group) {
      auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_proc);
      EKAT_REQUIRE_MSG(static_cast<bool>(group),
                         "Error! Something went is wrong with the atmosphere process\n" +
                         ("          " + atm_proc->name() + "\n") +
                         "       Its type returns 'Group', but casting to AtmosphereProcessGroup\n" +
                         "       returns a null shared_ptr.\n");
      group->set_field_repos (repo,bkp_repo);
    }
  }
}
#endif

void AtmosphereProcessGroup::
process_required_field (const FieldIdentifier& fid) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    // If the schedule is sequential, we do not add inputs if they are computed
    // by a previous process (they are not an 'input' to the group).
    if (!ekat::contains(m_computed_fields,fid)) {
      m_required_fields.insert(fid);
    }
  } else {
    // In parallel schedule, the inputs of all processes are inputs of the group
    m_required_fields.insert(fid);
  }
}

} // namespace scream
