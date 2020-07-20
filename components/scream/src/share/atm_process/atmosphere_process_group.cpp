#include "share/atm_process/atmosphere_process_group.hpp"

#include "ekat/util/scream_std_utils.hpp"
#include "ekat/util/string_utils.hpp"
#include "ekat/util/scream_std_utils.hpp"

namespace scream {

AtmosphereProcessGroup::
AtmosphereProcessGroup (const Comm& comm, const ParameterList& params)
 : m_comm(comm)
{
  // Get number of processes in the group and the scheduling type (Sequential vs Parallel)
  m_group_size = params.get<int>("Number of Entries");
  scream_require_msg (m_group_size>0, "Error! Invalid group size.\n");

  if (m_group_size>1) {
    if (params.get<std::string>("Schedule Type") == "Sequential") {
      m_group_schedule_type = ScheduleType::Sequential;
    } else if (params.get<std::string>("Schedule Type") == "Parallel") {
      m_group_schedule_type = ScheduleType::Parallel;
      error::runtime_abort("Error! Parallel schedule not yet implemented.\n");
    } else {
      error::runtime_abort("Error! Invalid 'Schedule Type'. Available choices are 'Parallel' and 'Sequential'.\n");
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
    Comm proc_comm = m_comm;
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
  // There are two catches when it comes to exposing a list of required/computed
  // fields.
  //  - catch 1: if (and only if) the splitting is sequential, and contains
  //    2+ processes, it is possible that some inputs of the non-first atm
  //    processes are computed by one (or more) of the previous ones.
  //    In this case, such fields are not really 'requirements' of this
  //    group, but simply "outputs" of the group (which happen to also
  //    be reused by other atm processes in the group).
  //  - catch 2: the APG is responsible of mapping input/output fields
  //    from/to the reference grid. For this reason, if one of the atm
  //    processes claims to need a field F1 which is defined on a grid other
  //    than the reference one, the APG will expose as a requirement F2,
  //    the same field defined on the reference grid, and will expose F1
  //    as a 'computed' field. Similarly, if a process computes F1 on a
  //    a grid other than the reference one, the APG will also expose
  //    a computed field F2, the same as F1 but defined on the reference grid.
  //    The extra computed fields will be "computed" by calling the
  //    proper AbstractRemapper implementation for the reference-actual
  //    grids pair. The APG stores internally a bunch of remappers,
  //    which calls during the `run` call, to perform the remappings.

  // Resizing the remappers vectors
  m_inputs_remappers.resize(m_group_size);
  m_outputs_remappers.resize(m_group_size);

  // The (bare) names of the computed fields. We use this to figure out
  // if an 'input' to a following process should be marked as input to
  // the APG as a whole, depending on whether one of the previous atm
  // processes in the group is computing the field.
  // NOTE: this is only in the case of sequential scheduling
  std::set<std::string> computed;

  // The reference grid
  const auto ref_grid = grids_manager->get_reference_grid();
  m_ref_grid_name = ref_grid->name();

  std::pair<decltype(m_required_fields)::iterator,bool>  it_bool;
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    auto atm_proc = m_atm_processes[iproc];
    atm_proc->set_grids(grids_manager);
    const auto proc_name = atm_proc->name();

    auto& remap_in = m_inputs_remappers[iproc];
    auto& remap_out = m_outputs_remappers[iproc];

    // Any atm proc sub-group will take care of remapping input/outputs.
    // We only have to deal with 'individual' atm processes
    const auto type = atm_proc->type();
    const bool is_group = type==AtmosphereProcessType::Group;
    if (!is_group) {
      const auto& atm_proc_grids = atm_proc->get_required_grids();
      // Note: in general, as of today, no atm process other than APG stores more than
      //       one grid, so this for loop is somewhat overkill. Still, it doesn't hurt.
      for (const auto& gname : atm_proc_grids) {
        const auto grid = grids_manager->get_grid(gname);
        remap_in[grid->name()] = grids_manager->create_remapper(ref_grid,grid);
        remap_out[grid->name()] = grids_manager->create_remapper(grid,ref_grid);
        remap_in[grid->name()]->registration_begins();
        remap_out[grid->name()]->registration_begins();
      }
    }

    // Add inputs to the list of inputs of the group
    for (const auto& fid : atm_proc->get_required_fields()) {
      // Whether this fid is on the reference grid or not, we expose
      // as input a copy of fid on the reference grid.
      // If fid is not on the reference grid, we add fid to our 'computed' fields.
      const bool is_ref_grid = (fid.get_grid_name()==ref_grid->name());
      const FieldIdentifier ref_fid = is_ref_grid ? fid : create_ref_fid(fid,remap_in[fid.get_grid_name()]);

      if (m_group_schedule_type==ScheduleType::Sequential) {
        // If the schedule is sequential, we do not add inputs if they are computed
        // by a previous process (they are not an 'input' to the group).
        if (computed.find(fid.name())==computed.end()) {
          it_bool = m_required_fields.insert(ref_fid);
        } else {
          m_internal_fields.insert(fid);
        }
      } else {
        // In parallel schedule, the inputs of all processes are inputs of the group
        it_bool = m_required_fields.insert(ref_fid);
      }

      // If this field id is not on the reference grid, we need to remap,
      // which means that we "compute" the input on this particular grid
      if (!is_ref_grid) {
        // We "compute" this field (on the non-reference grid) by means of remapping
        computed.insert(fid.name());
        m_computed_fields.insert(fid);

        remap_in[fid.get_grid_name()]->register_field(ref_fid,fid);
      }
    }

    // Add outputs to the list of outputs of the group
    for (const auto& fid : atm_proc->get_computed_fields()) {
      m_computed_fields.insert(fid);

      // Update the list of fields computed by the processes processed so far
      computed.insert(fid.name());

      // If the grid of fid is not the reference one, we also remap it
      // to the reference grid, hence "computing" ref_fid
      // Exception: if atm proc is a group, then it should already take care of this.
      if (!is_group && fid.get_grid_name()!=ref_grid->name()) {
        auto& remapper = remap_out[fid.get_grid_name()];

        const FieldIdentifier ref_fid = create_ref_fid(fid,remapper);
        m_computed_fields.insert(ref_fid);

        remapper->register_field(fid,ref_fid);
      }
    }

    // Close registration in the remappers
    for (auto it : m_inputs_remappers[iproc]) {
      it.second->registration_ends();
    }
    for (auto it : m_outputs_remappers[iproc]) {
      it.second->registration_ends();
    }
  }
}

void AtmosphereProcessGroup::initialize (const util::TimeStamp& t0) {
  // Now that we have the comm for the processes in the group, we can initialize them
  for (int i=0; i<m_group_size; ++i) {
    m_atm_processes[i]->initialize(t0);
  }
  m_current_ts = t0;
}

void AtmosphereProcessGroup::
setup_remappers (const FieldRepository<Real, device_type>& field_repo) {
  // Now that all fields have been set, we can set the fields in the remappers
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    for (auto& it : m_inputs_remappers[iproc]) {
      const auto& remapper = it.second;
      const int num_fields = remapper->get_num_fields();
      for (int ifield=0; ifield<num_fields; ++ifield) {
        const auto& src_id = remapper->get_src_field_id(ifield);
        const auto& tgt_id = remapper->get_tgt_field_id(ifield);

        auto& src = field_repo.get_field(src_id);
        auto& tgt = field_repo.get_field(tgt_id);

        remapper->bind_field(src,tgt);
      }
    }

    for (auto& it : m_outputs_remappers[iproc]) {
      const auto& remapper = it.second;
      const int num_fields = remapper->get_num_fields();
      for (int ifield=0; ifield<num_fields; ++ifield) {
        const auto& src_id = remapper->get_src_field_id(ifield);
        const auto& tgt_id = remapper->get_tgt_field_id(ifield);

        auto& src = field_repo.get_field(src_id);
        auto& tgt = field_repo.get_field(tgt_id);

        remapper->bind_field(src,tgt);
      }
    }
  }

  // If any of the stored atm procs is an AtmosphereProcessGroup,
  // propagate the call
  for (const auto& atm_proc : m_atm_processes) {
    if (atm_proc->type()==AtmosphereProcessType::Group) {
      auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_proc);
      scream_require_msg(static_cast<bool>(group),
                         "Error! Something went is wrong with the atmosphere process\n" +
                         ("          " + atm_proc->name() + "\n") +
                         "       Its type returns 'Group', but casting to AtmosphereProcessGroup\n" +
                         "       returns a null shared_ptr.\n");
      group->setup_remappers(field_repo);
    }
  }
}

void AtmosphereProcessGroup::run (const Real dt) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    run_sequential(dt);
  } else {
    run_parallel(dt);
  }
}

void AtmosphereProcessGroup::run_sequential (const Real dt) {
  m_current_ts += dt;
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    // Reamp the inputs of this process
    const auto& inputs_remappers = m_inputs_remappers[iproc];
    for (const auto& it : inputs_remappers) {
      const auto& remapper = it.second;
      remapper->remap(true);
      for (int ifield=0; ifield<remapper->get_num_fields(); ++ifield) {  
        const auto& f = remapper->get_tgt_field(ifield);
        f.get_header_ptr()->get_tracking().update_time_stamp(m_current_ts);
#ifdef SCREAM_DEBUG
        // Update the remapped field in the bkp repo
        auto bkp_f = m_bkp_field_repo->get_field(f.get_header().get_identifier());
        bkp_f.get_header_ptr()->get_tracking().update_time_stamp(m_current_ts);
        auto src = f.get_view();
        auto dst = bkp_f.get_view();
        Kokkos::deep_copy(dst,src);
#endif
      }
    }

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
          if (util::contains(computed,fid) ||
              (inputs_remappers.find(gn)!=inputs_remappers.end() &&
               inputs_remappers.at(gn)->has_tgt_field(fid))) {
            // For fields that changed, make sure the time stamp has been updated
            const auto& ts = f.second.get_header().get_tracking().get_time_stamp();
            scream_require_msg(field_is_unchanged || ts==m_current_ts,
                               "Error! Process '" + atm_proc->name() + "' updated field '" +
                                fid.get_id_string() + "', but it did not update its time stamp.\n");
          } else {
            // For non computed fields, make sure the field in m_repo matches
            // the copy in m_bkp_repo, and update the copy in m_bkp_repo.
            // There are three ok scenarios: 1) field is unchanged, 2) field is computed,
            // 3) field is a remapped input.
            scream_require_msg(field_is_unchanged,
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

    // Remap the outputs of this process
    for (const auto& it : m_outputs_remappers[iproc]) {
      const auto& remapper = it.second;
      remapper->remap(true);
      for (int ifield=0; ifield<remapper->get_num_fields(); ++ifield) {
        const auto& f = remapper->get_tgt_field(ifield);
        f.get_header_ptr()->get_tracking().update_time_stamp(m_current_ts);
#ifdef SCREAM_DEBUG
        // Update the remapped field in the bkp repo
        auto bkp_f = m_bkp_field_repo->get_field(f.get_header().get_identifier());
        bkp_f.get_header_ptr()->get_tracking().update_time_stamp(m_current_ts);
        auto src = f.get_view();
        auto dst = bkp_f.get_view();
        Kokkos::deep_copy(dst,src);
#endif
      }
    }
  }
}

void AtmosphereProcessGroup::run_parallel (const Real /* dt */) {
  scream_require_msg (false,"Error! Parallel splitting not yet implemented.\n");
}

void AtmosphereProcessGroup::finalize (/* what inputs? */) {
  for (auto atm_proc : m_atm_processes) {
    atm_proc->finalize(/* what inputs? */);
  }
}

void AtmosphereProcessGroup::register_fields (FieldRepository<Real, device_type>& field_repo) const {
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    const auto& atm_proc = m_atm_processes[iproc];
    atm_proc->register_fields(field_repo);

#ifdef SCREAM_DEBUG
    // Make sure processes are not calling methods they shouldn't on the repo
    scream_require_msg(field_repo.repository_state()==RepoState::Open,
                         "Error! Atmosphere processes are *not* allowed to modify the state of the repository.\n");

    // Check that the required fields are indeed in the repo now
    for (const auto& id : atm_proc->get_required_fields()) {
      scream_require_msg(field_repo.has_field(id), "Error! Process '" + atm_proc->name() + "' failed to register required field '" + id.get_id_string() + "'.\n");
    }
    for (const auto& id : atm_proc->get_computed_fields()) {
      scream_require_msg(field_repo.has_field(id), "Error! Process '" + atm_proc->name() + "' failed to register computed field '" + id.get_id_string() + "'.\n");
    }
#endif

    // Additionally, register all the fields that are needed just for remap reasons.
    // At least half of these are already registered, since they are inputs/outputs
    // of the stored atm processes. But there's no problem registering the same field
    // twice in the repo.
    for (auto it : m_inputs_remappers[iproc]) {
      const int num_fields = it.second->get_num_fields();
      for (int ifield=0; ifield<num_fields; ++ifield) {
        field_repo.register_field(it.second->get_src_field_id(ifield));
        field_repo.register_field(it.second->get_tgt_field_id(ifield));
      }
    }
    for (auto it : m_outputs_remappers[iproc]) {
      const int num_fields = it.second->get_num_fields();
      for (int ifield=0; ifield<num_fields; ++ifield) {
        field_repo.register_field(it.second->get_src_field_id(ifield));
        field_repo.register_field(it.second->get_tgt_field_id(ifield));
      }
    }
  }
}

void AtmosphereProcessGroup::set_required_field_impl (const Field<const Real, device_type>& f) {
  const auto& fid = f.get_header().get_identifier();
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->requires(fid)) {
      atm_proc->set_required_field(f);
    }
  }
}

void AtmosphereProcessGroup::set_computed_field_impl (const Field<Real, device_type>& f) {
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
set_field_repos (const FieldRepository<Real, device_type>& repo,
                 const FieldRepository<Real, device_type>& bkp_repo) {
  m_field_repo = &repo;
  m_bkp_field_repo = &bkp_repo;
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->type()==AtmosphereProcessType::Group) {
      auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_proc);
      scream_require_msg(static_cast<bool>(group),
                         "Error! Something went is wrong with the atmosphere process\n" +
                         ("          " + atm_proc->name() + "\n") +
                         "       Its type returns 'Group', but casting to AtmosphereProcessGroup\n" +
                         "       returns a null shared_ptr.\n");
      group->set_field_repos (repo,bkp_repo);
    }
  }
}

bool AtmosphereProcessGroup::
views_are_equal(const field_type& f1, const field_type& f2) {
  const auto& layout = f1.get_header().get_identifier().get_layout();
  const int size = layout.size();
  const int last_dim_alloc_size = f1.get_header().get_alloc_properties().get_last_dim_alloc_size();
  const int last_alloc_dim = last_dim_alloc_size / sizeof(Real);
  const int last_dim = layout.dim(layout.rank()-1);

  const auto v1 = f1.get_view();
  const auto v2 = f2.get_view();

  using exec_space = field_type::view_type::execution_space;
  Kokkos::RangePolicy<exec_space> policy(0,size);
  int num_diffs = 0;
  Kokkos::parallel_reduce(
    policy,
    [=](const int i, int& sum){
      const int outer_dim = i / last_dim;
      const int inner_dim = i % last_dim;
      const int idx = outer_dim*last_alloc_dim + inner_dim;

      if (v1(idx)!=v2(idx)) {
        ++sum;
      }
    },num_diffs
  );
  return (num_diffs == 0);
}
#endif

void AtmosphereProcessGroup::set_internal_field (const Field<Real, device_type>& f)
{
  const auto& fid = f.get_header().get_identifier();
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    auto proc  = m_atm_processes[iproc];
    auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(proc);
    if (proc->requires(fid)) {
      proc->set_required_field(f.get_const());
    }
    if (proc->computes(fid)) {
      proc->set_computed_field(f);
    }
    if (static_cast<bool>(group)) {
      if (util::contains(group->get_internal_fields(),fid)) {
        group->set_internal_field(f);
      }
    }
  }
}

FieldIdentifier
AtmosphereProcessGroup::create_ref_fid (const FieldIdentifier& fid,
                                        const remapper_ptr_type& remapper)
{
  const auto& fid_layout = fid.get_layout();
  const auto& grid_name = fid.get_grid_name();
  if (remapper->get_src_grid()->name()==grid_name) {
    const FieldLayout ref_layout = remapper->create_tgt_layout(fid_layout);
    return FieldIdentifier (fid.name(),ref_layout,fid.get_units(),remapper->get_tgt_grid()->name());
  } else if (remapper->get_tgt_grid()->name()==grid_name) {
    const FieldLayout ref_layout = remapper->create_src_layout(fid_layout);
    return FieldIdentifier (fid.name(),ref_layout,fid.get_units(),remapper->get_src_grid()->name());
  } else {
    // Something went wrong
    scream_require_msg(false,"Error! Input FieldIdentifier's grid name is neither the source nor the target grid of the input remapper.\n");
  }
}

} // namespace scream
