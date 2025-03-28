#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/field/field_utils.hpp"

#include "share/property_checks/field_nan_check.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <memory>

namespace scream {

AtmosphereProcessGroup::
AtmosphereProcessGroup (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Get the list of procs in the group
  auto group_list = m_params.get<std::vector<std::string>>("atm_procs_list");
  m_group_size = group_list.size();
  EKAT_REQUIRE_MSG (m_group_size>0, "Error! Invalid group size.\n");

  if (m_group_size>1) {
    if (m_params.get<std::string>("schedule_type") == "Sequential") {
      m_group_schedule_type = ScheduleType::Sequential;
    } else if (m_params.get<std::string>("schedule_type") == "Parallel") {
      m_group_schedule_type = ScheduleType::Parallel;
      ekat::error::runtime_abort("Error! Parallel schedule not yet implemented.\n");
    } else {
      ekat::error::runtime_abort("Error! Invalid 'schedule_type'. Available choices are 'Parallel' and 'Sequential'.\n");
    }
  } else {
    // Pointless to handle this group as parallel, if only one process is in it
    m_group_schedule_type = ScheduleType::Sequential;
  }

  // Create the individual atmosphere processes
  m_group_name = params.name();

  // Temp to store which atm proc needs each restart extra data. We'll use it to ensure
  // only one atm proc request a particular extra data key
  strmap_t<std::string> ed2proc;

  auto& apf = AtmosphereProcessFactory::instance();
  // Ensure the "Group" atm proc is registered in the factory,
  // so we can recursively create groups. Groups are an impl detail,
  // so we don't expect users to register the APG in the factory.
  apf.register_product("group",&create_atmosphere_process<AtmosphereProcessGroup>);
  for (const auto& ap_name : group_list) {
    // The comm to be passed to the processes construction is
    //  - the same as the comm of this APG, if num_entries=1 or sched_type=Sequential
    //  - a sub-comm of this APG's comm otherwise
    ekat::Comm proc_comm = m_comm;
    if (m_group_schedule_type==ScheduleType::Parallel) {
      // This is what's going to happen when we implment this:
      //  - the processes in the group are going to be run in parallel
      //  - each rank is assigned ONE atm process
      //  - all the atm processes not assigned to this rank will be filled with
      //    an instance of "RemoteProcessStub" (to be implemented),
      //    which is a do-nothing class, only responsible to keep track of dependencies
      //  - the input parameter list should specify for each atm process the number
      //    of mpi ranks dedicated to it. Obviously, these numbers should add up
      //    to the size of the input communicator.
      //  - this class is then responsible of 'combining' the results togehter,
      //    including remapping input/output fields to/from the sub-comm
      //    distribution.
      EKAT_ERROR_MSG("Error! Parallel schedule type not yet implemented.\n");
    }

    // Get the params of this atm proc
    auto& params_i = m_params.sublist(ap_name);

    // Get type (defaults to name)
    const auto& ap_type = params_i.get<std::string>("Type",ap_name);

    // Set logger in this ap params
    params_i.set("Logger",this->m_atm_logger);

    // Create the atm proc
    auto ap = apf.create(ap_type,proc_comm,params_i);
    m_atm_processes.push_back(ap);

    // NOTE: the shared_ptr of the new atmosphere process *MUST* have been created correctly.
    //       Namely, the creation process must have set up enable_shared_from_this's status correctly.
    //       This is done by the library-provided templated function 'create_atm_process<T>.
    //       However, if the user has decided to roll his/her own creator function to be registered
    //       in the AtmosphereProcessFactory, he/she may have forgot to set the self pointer in the process.
    //       To make sure this is not the case, we check that the weak_ptr in the newly created
    //       atmosphere process (which comes through inheritance from enable_shared_from_this) is valid.
    EKAT_REQUIRE_MSG(!ap->weak_from_this().expired(),
        "Error! The newly created std::shared_ptr<AtmosphereProcess> did not correctly setup\n"
        "       the 'enable_shared_from_this' interface.\n"
        "       Did you by chance register your own creator function in the AtmosphereProccessFactory class?\n"
        "       If so, don't. Instead, use the instantiation of create_atmosphere_process<T>,\n"
        "       with T = YourAtmProcessClassName.\n");

    // Store a copy of all the restart extra data of the atm proc.
    // NOTE: any uses std::shared_ptr internally, so if the atm proc updates
    //       the extra data, it will be updated in this class too.
    for (const auto& it : ap->get_restart_extra_data()) {
      // We don't want to risk having two processes overwriting restart data, in case
      // of a "common name" var (e.g., "num_steps"). Each process should try its best
      // to provide names that are likely to be unique. Even if two procs *actyally need*
      // the same var, we can write it twice to file.
      EKAT_REQUIRE_MSG (m_restart_extra_data.count(it.first)==0,
          "Error! Two atmosphere processes use the same restart extra data key.\n"
          " - Extra data key : " + it.first + "\n"
          " - First atm proc : " + ed2proc[it.first] + "\n"
          " - Second atm proc: " + ap->name() + "\n");
      m_restart_extra_data[it.first] = it.second;
      ed2proc[it.first] = ap->name();
    }
  }
}

std::shared_ptr<AtmosphereProcessGroup::atm_proc_type>
AtmosphereProcessGroup::get_process_nonconst (const std::string& name) const {
  EKAT_REQUIRE_MSG(has_process(name), "Error! Process "+name+" requested from group "+this->name()+
                                      " but is not constained in this group.\n");

  std::shared_ptr<atm_proc_type> return_process;
  for (auto& process: m_atm_processes) {
    if (process->type() == AtmosphereProcessType::Group) {
      const auto* group = dynamic_cast<const AtmosphereProcessGroup*>(process.get());
      return group->get_process_nonconst(name);
    } else if (process->name() == name) {
        return_process = process;
        break;
    }
  }
  return return_process;
}

bool AtmosphereProcessGroup::has_process(const std::string& name) const {
  for (auto& process: m_atm_processes) {
    if (process->type() == AtmosphereProcessType::Group) {
      const auto* group = dynamic_cast<const AtmosphereProcessGroup*>(process.get());
      if (group->has_process(name)) {
        return true;
      }
    } else {
      if (process->name() == name) {
        return true;
      }
    }
  }
  return false;
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

  m_grids_mgr = grids_manager;
}

void AtmosphereProcessGroup::setup_tendencies_requests () {
  auto is_tend = [](const std::string& name) -> bool
  {
    return name.size()>5 &&
           name.substr(name.size()-5)=="_tend";
  };

  AtmosphereProcess::setup_tendencies_requests();
  for (const auto& atm_proc : m_atm_processes) {
    atm_proc->setup_tendencies_requests();

    // Redo the add_field<Computed> for all XYZ_tend fields, since
    // they were added *after* the call to set_grids.
    for (const auto& req : atm_proc->get_computed_field_requests()) {
      if (is_tend(req.fid.name())) {
        add_field<Computed>(req);
      }
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

bool AtmosphereProcessGroup::
are_column_conservation_checks_enabled () const
{
  // Loop through processes and return true if an instance is found.
  for (auto atm_proc : m_atm_processes) {

    // Recursively check processes of internal groups. If one is found, return true,
    // else continue to the next process.
    auto atm_proc_group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_proc);
    if (atm_proc_group) {
      if (atm_proc_group->are_column_conservation_checks_enabled()) {
        return true;
      } else {
        continue;
      }
    }

    // If this process is not a group, query enable_column_conservation_checks
    // and return true if true.
    if (atm_proc->has_column_conservation_check()) {
      return true;
    }
  }

  // If no process was found with enable_column_conservation_checks=true, return false.
  return false;
}

void AtmosphereProcessGroup::
setup_column_conservation_checks (const std::shared_ptr<MassAndEnergyColumnConservationCheck>& conservation_check,
                                  const CheckFailHandling                                      fail_handling_type) const
{
  // Loop over atm processes and add mass and energy checker where relevant
  for (auto atm_proc : m_atm_processes) {

    // We only want the checks on an individual processes. If
    // atm_proc is a group, recursively call this function.
    auto atm_proc_group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(atm_proc);
    if (atm_proc_group) {
      // Require that a group does not enable conservation checks.
      // Individual processes are the only ones that can define
      // this check as there is currenty no concept of boundary
      // fluxes over multiple processes implemented in the model.
      EKAT_REQUIRE_MSG(not atm_proc_group->has_column_conservation_check(),
                       "Error! The ATM process group \"" + atm_proc_group->name() + "\" attempted to enable "
                       "conservation checks. Should have enable_column_conservation_checks=false for all "
                       "process groups.\n");

      atm_proc_group->setup_column_conservation_checks(conservation_check, fail_handling_type);
      continue;
    }

    // For individual processes, first query if the checks are enabled.
    // If not, continue to the next process.
    if (not atm_proc->has_column_conservation_check()) {
      continue;
    }

    // Since the checker is column local, require that an atm
    // process that enables the check is a Physics process.
    EKAT_REQUIRE_MSG(atm_proc->type() == AtmosphereProcessType::Physics,
                     "Error! enable_column_conservation_checks=true "
                     "for non-physics process \"" + atm_proc->name() + "\". "
                     "This check is column local and therefore can only be run "
                     "on physics processes.\n");

    // Query the computed fields for this atm process and see if either the mass or energy computation
    // might be changed after the process has run. If no field used in the mass or energy calculate
    // is updated by this process, there is no need to run the check.
    const std::string phys_grid_name  = conservation_check->get_grid()->name();
    const bool updates_static_energy  = atm_proc->has_computed_field("T_mid", phys_grid_name);
    const bool updates_kinetic_energy = atm_proc->has_computed_field("horiz_winds", phys_grid_name);
    const bool updates_water_vapor    = atm_proc->has_computed_field("qv", phys_grid_name);
    const bool updates_water_liquid   = atm_proc->has_computed_field("qc", phys_grid_name) ||
                                        atm_proc->has_computed_field("qr", phys_grid_name);
    const bool updates_water_ice      = atm_proc->has_computed_field("qi", phys_grid_name);
    const bool mass_or_energy_is_updated = updates_static_energy || updates_kinetic_energy ||
                                           updates_water_vapor   || updates_water_liquid ||
                                           updates_water_ice;
    EKAT_REQUIRE_MSG(mass_or_energy_is_updated, "Error! enable_column_conservation_checks=true for "
                                                "process \"" + atm_proc->name() + "\" but mass or energy is "
                                                "not updated by the process. Set to false to avoid "
                                                "unnecessary computation.\n");

    // Require that, if a process adds the conservation check, it also defines all
    // the boundary fluxes needed to compute the mass and energy tendencies.
    const bool has_all_boundary_fluxes = atm_proc->has_computed_field("vapor_flux", phys_grid_name) &&
                                         atm_proc->has_computed_field("water_flux", phys_grid_name) &&
                                         atm_proc->has_computed_field("ice_flux",   phys_grid_name) &&
                                         atm_proc->has_computed_field("heat_flux",  phys_grid_name);
    EKAT_REQUIRE_MSG(has_all_boundary_fluxes,
                     "Error! Process \"" + atm_proc->name() + "\" enables the mass "
                     "and energy conservation check, but does not define all "
                     "the boundary fluxes required: vapor_flux, water_flux "
                     "ice_flux, heat_flux. If a flux does not have a natural definition "
                     "within the process, set to 0.\n");

    // If all conditions are satisfied, add as postcondition_check
    atm_proc->add_column_conservation_check(conservation_check, fail_handling_type);
  }
}

void AtmosphereProcessGroup::add_postcondition_nan_checks () const {
  for (auto proc : m_atm_processes) {
    auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(proc);
    if (group) {
      group->add_postcondition_nan_checks();
    } else {
      for (const auto& f : proc->get_fields_out()) {
        const auto& grid_name = f.get_header().get_identifier().get_grid_name();
        auto nan_check = std::make_shared<FieldNaNCheck>(f,m_grids_mgr->get_grid(grid_name));
        proc->add_postcondition_check(nan_check, CheckFailHandling::Fatal);
      }

      for (const auto& g : proc->get_groups_out()) {
        const auto& grid = m_grids_mgr->get_grid(g.grid_name());
        for (const auto& f : g.m_individual_fields) {
          auto nan_check = std::make_shared<FieldNaNCheck>(*f.second,grid);
          proc->add_postcondition_check(nan_check, CheckFailHandling::Fatal);
        }
      }
    }
  }
}

void AtmosphereProcessGroup::add_additional_data_fields_to_property_checks (const Field& data_field) {
  for (auto proc : m_atm_processes) {
    auto group = std::dynamic_pointer_cast<AtmosphereProcessGroup>(proc);
    if (group) {
      group->add_additional_data_fields_to_property_checks(data_field);
    } else {
      for (auto& prop_check : proc->get_precondition_checks()) {
        prop_check.second->set_additional_data_field(data_field);
      }
      for (auto& prop_check : proc->get_postcondition_checks()) {
        prop_check.second->set_additional_data_field(data_field);
      }
      if (proc->has_column_conservation_check()) {
        proc->get_column_conservation_check().second->set_additional_data_field(data_field);
      }
    }
  }
}

void AtmosphereProcessGroup::pre_process_tracer_requests () {
  // Create map from tracer name to a vector which contains the field requests for that tracer.
  std::map<std::string, std::list<FieldRequest>> tracer_requests;
  auto gather_tracer_requests = [&] (FieldRequest& req) {
    if (ekat::contains(req.groups, "tracers")) {
      tracer_requests[req.fid.name()].push_back(req);
    }
  };
  for (auto& req : m_required_field_requests){
    gather_tracer_requests(req);
  }
  for (auto& req : m_computed_field_requests) {
    gather_tracer_requests(req);
  }

  // Go through the map entry for each tracer and check that every one
  // has the same request for turbulence advection, or listed no preference.
  std::map<std::string, std::string> tracer_advection_type;
  for (auto& fr : tracer_requests) {
    auto& reqs = fr.second;

    bool turb_advect_req = false;
    bool non_turb_advect_req = false;
    for (auto& req : reqs) {
      if (ekat::contains(req.groups, "turbulence_advected_tracers")) turb_advect_req = true;
      if (ekat::contains(req.groups, "non_turbulence_advected_tracers")) non_turb_advect_req = true;
      if (turb_advect_req and non_turb_advect_req) {
        // All the info we need to error out, just break request loop
        break;
      }
    }

    if (turb_advect_req and non_turb_advect_req) {
      std::ostringstream ss;
      ss << "Error! Incompatible tracer request. Turbulence advection requests not consistent among processes.\n"
            "  - Tracer name: " + fr.first + "\n"
            "  - Requests (process name, grid name, turbulence advected request type):\n";
      for (auto& req : reqs) {
        const auto turb_advect = ekat::contains(req.groups, "turbulence_advected_tracers");
        const auto non_turb_advect = ekat::contains(req.groups, "non_turbulence_advected_tracers");
        std::string turb_advect_info =
          (turb_advect ?       "DynamicsAndTurbulence" :
            (non_turb_advect ? "DynamicsOnly" :
                                "NoPreference"));
        const auto grid_name = req.fid.get_grid_name();
        ss << "    - (" + req.calling_process + ", " + grid_name + ", " + turb_advect_info + ")\n";
      }
      EKAT_ERROR_MSG(ss.str());
    } else if (non_turb_advect_req) {
      tracer_advection_type[fr.first] = "non_turbulence_advected_tracers";
    } else {
      tracer_advection_type[fr.first] = "turbulence_advected_tracers";
    }
  }

  // Set correct tracer advection type (if not set)
  auto add_correct_tracer_group = [&] (FieldRequest& req, const std::string& advect_type) {
    if (not ekat::contains(req.groups, advect_type)) {
      req.groups.push_back(advect_type);
    }
  };
  for (auto& req : m_required_field_requests){
    if (ekat::contains(req.groups, "tracers")) {
      add_correct_tracer_group(req, tracer_advection_type.at(req.fid.name()));
    }
  }
  for (auto& req : m_computed_field_requests) {
    if (ekat::contains(req.groups, "tracers")) {
      add_correct_tracer_group(req, tracer_advection_type.at(req.fid.name()));
    }
  }
}

void AtmosphereProcessGroup::initialize_impl (const RunType run_type) {
  for (auto& atm_proc : m_atm_processes) {
    atm_proc->initialize(start_of_step_ts(),run_type);
#ifdef SCREAM_HAS_MEMORY_USAGE
    long long my_mem_usage = get_mem_usage(MB);
    long long max_mem_usage;
    m_comm.all_reduce(&my_mem_usage,&max_mem_usage,1,MPI_MAX);
    m_atm_logger->debug("[EAMxx::initialize::"+atm_proc->name()+"] memory usage: " + std::to_string(max_mem_usage) + "MB");
#endif
  }
}

void AtmosphereProcessGroup::run_impl (const double dt) {
  if (m_group_schedule_type==ScheduleType::Sequential) {
    run_sequential(dt);
  } else {
    run_parallel(dt);
  }
}

void AtmosphereProcessGroup::run_sequential (const double dt) {
  // The stored atm procs should update the timestamp if both
  //  - this is the last subcycle iteration
  //  - nobody from outside told this APG to not update timestamps
  const bool do_update = do_update_time_stamp() &&
                      (get_subcycle_iter()==get_num_subcycles()-1);
  for (auto atm_proc : m_atm_processes) {
    atm_proc->set_update_time_stamps(do_update);
    // Run the process
    atm_proc->run(dt);
#ifdef SCREAM_HAS_MEMORY_USAGE
    long long my_mem_usage = get_mem_usage(MB);
    long long max_mem_usage;
    m_comm.all_reduce(&my_mem_usage,&max_mem_usage,1,MPI_MAX);
    m_atm_logger->debug("[EAMxx::run_sequential::"+atm_proc->name()+"] memory usage: " + std::to_string(max_mem_usage) + "MB");
#endif
  }
}

void AtmosphereProcessGroup::run_parallel (const double /* dt */) {
  EKAT_REQUIRE_MSG (false,"Error! Parallel splitting not yet implemented.\n");
}

void AtmosphereProcessGroup::finalize_impl (/* what inputs? */) {
  for (auto atm_proc : m_atm_processes) {
    atm_proc->finalize(/* what inputs? */);
#ifdef SCREAM_HAS_MEMORY_USAGE
    long long my_mem_usage = get_mem_usage(MB);
    long long max_mem_usage;
    m_comm.all_reduce(&my_mem_usage,&max_mem_usage,1,MPI_MAX);
    m_atm_logger->debug("[EAMxx::finalize::"+atm_proc->name()+"] memory usage: " + std::to_string(max_mem_usage) + "MB");
#endif
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
  for (const auto& it : group.m_individual_fields) {
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
