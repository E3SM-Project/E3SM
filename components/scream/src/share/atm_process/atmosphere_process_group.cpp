#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <memory>

namespace scream {

AtmosphereProcessGroup::
AtmosphereProcessGroup (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Get number of processes in the group and the scheduling type (Sequential vs Parallel)
  m_group_size = params.get<int>("Number of Entries");
  EKAT_REQUIRE_MSG (m_group_size>0, "Error! Invalid group size.\n");

  if (m_group_size>1) {
    if (m_params.get<std::string>("Schedule Type") == "Sequential") {
      m_group_schedule_type = ScheduleType::Sequential;
    } else if (m_params.get<std::string>("Schedule Type") == "Parallel") {
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

  auto& apf = AtmosphereProcessFactory::instance();
  apf.register_product("group",&create_atmosphere_process<AtmosphereProcessGroup>);
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

    // Set the logger in this AP's params
    auto& params_i = m_params.sublist(ekat::strint("Process",i));
    params_i.set("Logger",this->m_atm_logger);
    const std::string& process_name = params_i.get<std::string>("Process Name");

    m_atm_processes.emplace_back(apf.create(process_name,proc_comm,params_i));

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
    for (const auto& req : atm_proc->get_required_field_requests()) {
      process_required_field(req);
    }
    for (const auto& req : atm_proc->get_computed_field_requests()) {
      add_field<Computed>(req);
    }
    for (const auto& req : atm_proc->get_required_group_requests()) {
      process_required_group(req);
    }
    for (const auto& req : atm_proc->get_computed_group_requests()) {
      add_group<Computed>(req);
    }
  }
}

void AtmosphereProcessGroup::
gather_internal_fields  () {
  // For debug purposes
  std::map<std::string,std::string> f2proc;
  for (auto& atm_proc : m_atm_processes) {
    auto apg = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_proc);
    if (apg) {
      // Have this apg build its list of internal fields first
      apg->gather_internal_fields();
    }
    const auto& ifs = atm_proc->get_internal_fields();
    for (const auto& f : ifs) {
      const auto& name = f.get_header().get_identifier().name();
      EKAT_REQUIRE_MSG (f2proc.find(name)==f2proc.end(),
          "Error! Two atm procs created the same internal field.\n"
          "  - field name: " + name + "\n"
          "  - first atm proc: " + f2proc.at(name) + "\n"
          "  - second atm proc: " + atm_proc->name() + "\n");
      add_internal_field(f);
      f2proc[name] = atm_proc->name();
    }
  }
}

void AtmosphereProcessGroup::initialize_impl (const RunType run_type) {
  for (auto& atm_proc : m_atm_processes) {
    atm_proc->initialize(timestamp(),run_type);
  }
}

void AtmosphereProcessGroup::run_impl (const int dt) {
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

  // The stored atm procs should update the timestamp if both
  //  - this is the last subcycle iteration
  //  - nobody from outside told this APG to not update timestamps
  const bool do_update = do_update_time_stamp() &&
                      (get_subcycle_iter()==get_num_subcycles()-1);
  for (auto atm_proc : m_atm_processes) {
    atm_proc->set_update_time_stamps(do_update);
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
set_required_field (const Field& f) {
  if (m_group_schedule_type==ScheduleType::Parallel) {
    // In parallel splitting, all required fields are *actual* inputs,
    // and the base class impl is fine.
    AtmosphereProcess::set_required_field(f);
  }

  // Find the first process that requires this group
  const auto& fid = f.get_header().get_identifier();
  int first_proc_that_needs_f = -1;
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    if (m_atm_processes[iproc]->has_required_field(fid)) {
      first_proc_that_needs_f = iproc;
      break;
    }
  }

  EKAT_REQUIRE_MSG (first_proc_that_needs_f>=0,
    "Error! This atmosphere process does not require the input field.\n"
    "    field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  // Loop over all the fields in the FieldGroup (FG), and see if they are all computed
  // by a previous process. If so, then this FG is not a "required" FG
  // for this AtmosphereProcessGroup.
  // NOTE: the case where the FG itself is computed by a previous atm proc
  //       is already handled during `set_grids` (see process_required_group).
  // NOTE: we still check also the groups computed by the previous procs,
  //       in case they contain some of the fields in this group.
  bool computed = false;
  for (int iproc=0; iproc<first_proc_that_needs_f; ++iproc) {
    if (m_atm_processes[iproc]->has_computed_field(fid)) {
      computed = true;
      goto endloop;
    }

    // Check if this proc computes fn as part of another group
    const auto& out_g = m_atm_processes[iproc]->get_groups_out();
    for (const auto& g : out_g) {
      if (g.grid_name()!=fid.get_grid_name()) {
        continue;
      }
      for (const auto& fn : g.m_info->m_fields_names) {
        if (fid.name()==fn) {
          computed = true;
          goto endloop;
        }
      }
    }

  }
  // Yes, goto statement are not the prettiest C++ feature, but they allow to
  // break from multiple nested loops without the need of several checks.
endloop:
  if (computed) {
    // We compute the field *before* any atm proc that requires it. Simply set it
    // in the atm procs that require it, without storing it as an input to this group.
    set_required_field_impl(f);
  } else {
    // The field is *not* computed before the 1st atm proc that requires it. We can
    // resort to the base class' version of 'set_required_field', which will store
    // the field as a required field of the whole AtmProgGroup.
    AtmosphereProcess::set_required_field(f);
  }
}

void AtmosphereProcessGroup::
set_required_group (const FieldGroup& group) {
  if (m_group_schedule_type==ScheduleType::Parallel) {
    // In parallel splitting, all required group are *actual* inputs,
    // and the base class impl is fine.
    AtmosphereProcess::set_required_group(group);
  }

  // Find the first process that requires this group
  int first_proc_that_needs_group = -1;
  for (int iproc=0; iproc<m_group_size; ++iproc) {
    if (m_atm_processes[iproc]->has_required_group(group.m_info->m_group_name,group.grid_name())) {
      first_proc_that_needs_group = iproc;
      break;
    }
  }

  EKAT_REQUIRE_MSG (first_proc_that_needs_group>=0,
    "Error! This atmosphere process does not require the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  // Loop over all the fields in the group, and see if they are all computed
  // by a previous process. If so, then this group is not a "required" group
  // for this group.
  // NOTE: the case where the group itself is computed by a previous atm proc
  // is already handled during `set_grids` (see process_required_group).
  // NOTE: we still check also the groups computed by the previous procs,
  //       in case they contain some of the fields in this group.
  std::set<std::string> computed;
  for (const auto& it : group.m_fields) {
    const auto& fn = it.first;
    const auto& fid = it.second->get_header().get_identifier();
    for (int iproc=0; iproc<first_proc_that_needs_group; ++iproc) {
      if (m_atm_processes[iproc]->has_computed_field(fid)) {
        computed.insert(fid.name());
        goto endloop;
      }

      // Check if this proc computes fn as part of another group
      const auto& out_g = m_atm_processes[iproc]->get_groups_out();
      for (const auto& g : out_g) {
        if (g.grid_name()!=group.grid_name()) {
          continue;
        }
        for (const auto& f : g.m_info->m_fields_names) {
          if (f==fn) {
            computed.insert(fn);
            goto endloop;
          }
        }
      }

    }
  // Yes, goto statement are not the prettiest C++ feature, but they allow to
  // break from multiple nested loops without the need of several checks.
endloop:
    if (computed.size()==group.m_info->m_fields_names.size()) {
      break;
    }
  }

  if (computed.size()==group.m_info->m_fields_names.size()) {
    // We compute all the fields in the group *before* any atm proc
    // that requires this group. Simply set the group in those atm
    // proc that require it, without storing it as an input to this group.
    set_required_group_impl(group);
  } else {
    // There's at least one field of the group that is not computed
    // before the 1st atm proc that requires this group. We can resort
    // to the base class' version of 'set_required_group', which will store
    // the group as a required field group of the whole AtmProgGroup.
    AtmosphereProcess::set_required_group(group);
  }
}

void AtmosphereProcessGroup::
set_required_group_impl (const FieldGroup& group)
{
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->has_required_group(group.m_info->m_group_name,group.grid_name())) {
      atm_proc->set_required_group(group);
    }
  }
}

void AtmosphereProcessGroup::
set_computed_group_impl (const FieldGroup& group)
{
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->has_computed_group(group.m_info->m_group_name,group.grid_name())) {
      atm_proc->set_computed_group(group);
    }
    // In sequential scheduling, some groups may be computed by
    // a process and used by the next one. In this case, the group
    // may not figure as 'input' for the group, but we still
    // need to set it in the processes that need it.
    if (atm_proc->has_required_group(group.m_info->m_group_name,group.grid_name())) {
      atm_proc->set_required_group(group.get_const());
    }
  }
}

void AtmosphereProcessGroup::set_required_field_impl (const Field& f) {
  const auto& fid = f.get_header().get_identifier();
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->has_required_field(fid)) {
      atm_proc->set_required_field(f);
    }
  }
}

void AtmosphereProcessGroup::set_computed_field_impl (const Field& f) {
  const auto& fid = f.get_header().get_identifier();
  for (auto atm_proc : m_atm_processes) {
    if (atm_proc->has_computed_field(fid)) {
      atm_proc->set_computed_field(f);
    }
    // In sequential scheduling, some fields may be computed by
    // a process and used by the next one. In this case, the field
    // does not figure as 'input' for the group, but we still
    // need to set it in the processes that need it.
    if (atm_proc->has_required_field(fid)) {
      atm_proc->set_required_field(f.get_const());
    }
  }
}

void AtmosphereProcessGroup::
process_required_group (const GroupRequest& req) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    if (has_computed_group(req.name,req.grid)) {
      // Some previous atm proc computes this group, so it's not an 'input'
      // of the atm group as a whole. However, we might need a different
      // pack size. So, instead of adding to the required groups,
      // we add to the computed ones. This way we don't modify the inputs
      // of the group, and still manage to communicate to the AD the pack size
      // that we need.
      // NOTE; we don't have a way to check if all the fields in the group
      //       are computed by previous processes, since we don't have
      //       the list of all fields in this group.
      add_group<Computed>(req);
    } else {
      add_group<Required>(req);
    }
  } else {
    // In parallel schedule, the inputs of all processes are inputs of the group
    add_group<Required>(req);
  }
}

void AtmosphereProcessGroup::
process_required_field (const FieldRequest& req) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    if (has_computed_field(req.fid)) {
      // Some previous atm proc computes this field, so it's not an 'input'
      // of the group as a whole. However, we might need a different pack size,
      // or want to add it to a different group. So, instead of adding to
      // the required fields, we add to the computed fields. This way we
      // don't modify the inputs of the group, and still manage to communicate
      // to the AD the pack size and group affiliations that we need
      add_field<Computed>(req);
    } else {
      add_field<Required>(req);
    }
  } else {
    // In parallel schedule, the inputs of all processes are inputs of the group
    add_field<Required>(req);
  }
}

size_t AtmosphereProcessGroup::requested_buffer_size_in_bytes () const
{
  size_t buf_size = 0;
  for (const auto& proc : m_atm_processes) {
    buf_size = std::max(buf_size,proc->requested_buffer_size_in_bytes());
  }

  return buf_size;
}

void AtmosphereProcessGroup::
init_buffers(const ATMBufferManager& buffer_manager) {
  for (auto& atm_proc : m_atm_processes) {
    atm_proc->init_buffers(buffer_manager);
  }
}

} // namespace scream
