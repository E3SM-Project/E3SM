#include "share/atm_process/atmosphere_process.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/property_checks/mass_and_energy_conservation_check.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_utils.hpp"

#ifdef EAMXX_HAS_PYTHON
#include "share/field/field_pyutils.hpp"
#include "share/core/eamxx_pysession.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

#include <ekat_assert.hpp>

#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>

ekat::logger::LogLevel str2LogLevel (const std::string& s) {
  using namespace ekat::logger;

  ekat::logger::LogLevel log_level;
  if (s=="off") {
    log_level = LogLevel::off;
  } else if (s=="trace") {
    log_level = LogLevel::trace;
  } else if (s=="debug") {
    log_level = LogLevel::debug;
  } else if (s=="info") {
    log_level = LogLevel::info;
  } else if (s=="warn") {
    log_level = LogLevel::warn;
  } else {
    EKAT_ERROR_MSG ("Invalid string value for log level: " + s + "\n");
  }
  return log_level;
};

namespace scream
{

AtmosphereProcess::
AtmosphereProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
  : m_comm       (comm)
  , m_params     (params)
{
  if (m_params.isParameter("logger")) {
    m_atm_logger = m_params.get<std::shared_ptr<logger_t>>("logger");
  } else {
    // Create a console-only logger, that logs all ranks
    m_atm_logger = console_logger(str2LogLevel(m_params.get<std::string>("log_level","trace")));
  }

  if (m_params.isParameter("number_of_subcycles")) {
    m_num_subcycles = m_params.get<int>("number_of_subcycles");
  }
  EKAT_REQUIRE_MSG (m_num_subcycles>0,
      "Error! Invalid number of subcycles in param list " + m_params.name() + ".\n"
      "  - Num subcycles: " + std::to_string(m_num_subcycles) + "\n");

  m_timer_prefix = m_params.get<std::string>("timer_prefix","EAMxx::");

  m_repair_log_level = str2LogLevel(m_params.get<std::string>("repair_log_level","warn"));

  // Info for mass and energy conservation checks
  m_conservation_data.has_column_conservation_check =
      m_params.get<bool>("enable_column_conservation_checks", false);

  // Energy fixer
  m_conservation_data.has_energy_fixer =
      m_params.get<bool>("enable_energy_fixer", false);
  m_conservation_data.has_air_sea_surface_water_thermo_fixer = 
      m_params.get<bool>("enable_air_sea_surface_water_thermo_fixer", false);

  // Energy fixer
  m_conservation_data.has_energy_fixer_debug_info =
      m_params.get<bool>("enable_energy_fixer_debug_info", false);

  //either energy fixer or column checks, but not both at the same time
  EKAT_REQUIRE_MSG ( 
      !(m_conservation_data.has_energy_fixer && m_conservation_data.has_column_conservation_check),
      "Error! In param list " + m_params.name() + " both enable_energy_fixer and"
           " enable_column_conservation_check are on, which is not allowed. \n");

  // air sea surface water thermo fixer can only be TRUE if energy fixer is TRUE
  EKAT_REQUIRE_MSG (
      !((!m_conservation_data.has_energy_fixer) && m_conservation_data.has_air_sea_surface_water_thermo_fixer),
      "Error! In param list " + m_params.name() + " enable_energy_fixer is false and"
           " enable_air_sea_surface_water_thermo_fixer is true, which is not allowed. \n");
  // energy fixer debug info can only be TRUE if energy fixer is TRUE
  EKAT_REQUIRE_MSG (
      !((!m_conservation_data.has_energy_fixer) && m_conservation_data.has_energy_fixer_debug_info),
      "Error! In param list " + m_params.name() + " enable_energy_fixer is false and"
           " enable_energy_fixer_debug_info is true, which is not allowed. \n");

  m_internal_diagnostics_level = m_params.get<int>("internal_diagnostics_level", 0);
#ifdef EAMXX_HAS_PYTHON
  if (m_params.get("py_module_name",std::string(""))!="") {
    auto& pysession = PySession::get();
    pysession.initialize();

    const auto& py_module_name = m_params.get<std::string>("py_module_name");
    const auto& py_module_path = m_params.get<std::string>("py_module_path","./");

    pysession.add_path(py_module_path);
    m_py_module = pysession.safe_import(py_module_name);
  }
#endif
}

void AtmosphereProcess::set_grids (const std::shared_ptr<const GridsManager> grids_manager) {
  m_grids_manager = grids_manager;
  
  // Initialize the field managers
  m_inputs = std::make_shared<FieldManager>(grids_manager);
  m_inputs->set_name(name() + "::inputs");
  
  m_outputs = std::make_shared<FieldManager>(grids_manager);
  m_outputs->set_name(name() + "::outputs");
  
  m_internals = std::make_shared<FieldManager>(grids_manager);
  m_internals->set_name(name() + "::internals");
  
  // Close the repos since we'll be adding fields directly
  m_inputs->registration_ends();
  m_outputs->registration_ends();
  m_internals->registration_ends();
  
  create_requests();
}

void AtmosphereProcess::initialize (const TimeStamp& t0, const RunType run_type) {
  if (this->type()!=AtmosphereProcessType::Group) {
    start_timer (m_timer_prefix + this->name() + "::init");
  }

  // Avoid logging and flushing if ap type is diag ...
  // ... because we could have 100+ of those in production runs
  if (this->type()!=AtmosphereProcessType::Diagnostic) {
    log (LogLevel::info,"  Initializing " + name() + "...");
    m_atm_logger->flush(); // During init, flush often (to help debug crashes)
  }

  m_start_of_step_ts = m_end_of_step_ts = t0;
  initialize_impl(run_type);

  // Avoid logging and flushing if ap type is diag ...
  // ... because we could have 100+ of those in production runs
  if (this->type()!=AtmosphereProcessType::Diagnostic) {
    log (LogLevel::info,"  Initializing " + name() + "... done!");
    m_atm_logger->flush(); // During init, flush often (to help debug crashes)
  }

  m_is_initialized = true;

  if (this->type()!=AtmosphereProcessType::Group) {
    stop_timer (m_timer_prefix + this->name() + "::init");
  }
}

void AtmosphereProcess::run (const double dt) {
  m_atm_logger->debug("[EAMxx::" + this->name() + "] run...");
  start_timer (m_timer_prefix + this->name() + "::run");
  if (m_params.get("enable_precondition_checks", true)) {
    // Run 'pre-condition' property checks stored in this AP
    run_precondition_checks();
  }

  // Let the derived class do the actual run
  auto dt_sub = dt / m_num_subcycles;

  // Init single step tendencies (if any) with current value of output field
  init_step_tendencies ();

  for (m_subcycle_iter=0; m_subcycle_iter<m_num_subcycles; ++m_subcycle_iter) {
    m_start_of_step_ts = m_end_of_step_ts;
    m_end_of_step_ts += dt_sub;

    //energy fixer needs both mass and energy fields that are computed by compute_column_conservation_checks_data(dt_sub)
    //but this will change with cp* (with cpdry heat_glob const is much easier to compute)
    //actually we do not need mass_before, so mass_before is redundant

    //however, we do not want to use column checks if the fixer is on. without a correction 
    //from the fixer, column checks would fail after fixer. 
    //we decided that using the column check to verify the fixer (after intorducing a new "flux")
    //is not that practical. the fixer will be verified by another call to global integral,
    //but disabled by default.
    if (has_column_conservation_check() || has_energy_fixer()) {
      // Column local mass and energy checks requires the total mass and energy
      // to be computed directly before the atm process is run, as well and store
      // the correct timestep for the process.
      compute_column_conservation_checks_data(dt_sub);
    }

    if (m_internal_diagnostics_level > 0)
      // Print hash of INPUTS before run
      print_global_state_hash(name() + "-pre-sc-" + std::to_string(m_subcycle_iter),
                              m_start_of_step_ts, true);

    // Run derived class implementation
    run_impl(dt_sub);

    if (m_internal_diagnostics_level > 0)
      // Print hash of OUTPUTS/INTERNALS after run
      print_global_state_hash(name() + "-pst-sc-" + std::to_string(m_subcycle_iter),
                              m_end_of_step_ts, false);

    if (has_energy_fixer()){
      const bool water_thermo_fixer = has_air_sea_surface_water_thermo_fixer();
      const bool debug_info = has_energy_fixer_debug_info();
      fix_energy(dt_sub, water_thermo_fixer, debug_info);
    }

    if (has_column_conservation_check()) {
      // Run the column local mass and energy conservation checks
      run_column_conservation_check();
    }
  }

  // Complete tendency calculations (if any)
  compute_step_tendencies();

  if (m_params.get("enable_postcondition_checks", true)) {
    // Run 'post-condition' property checks stored in this AP
    run_postcondition_checks();
  }

  if (m_update_time_stamps) {
    // Update all output fields time stamps
    update_time_stamps ();
  }
  stop_timer (m_timer_prefix + this->name() + "::run");
}

void AtmosphereProcess::finalize () {
  finalize_impl(/* what inputs? */);
#ifdef EAMXX_HAS_PYTHON
  if (m_py_module.has_value()) {
    // Empty these vars before finalizing the py session, or else
    // their destructor will NOT find an active py interpreter
    m_py_module.reset();
    m_py_fields_dev.clear();
    m_py_fields_host.clear();

    // Note: In case multiple places have called PySession::get().initialize(),
    // only the last call to finalize *actually* finalizes the interpreter
    PySession::get().finalize();
  }
#endif
}

void AtmosphereProcess::setup_step_tendencies (const std::string& default_grid) {
  using strvec_t = std::vector<std::string>;
  auto tend_vec = m_params.get<strvec_t>("compute_tendencies",{});
  if (tend_vec.size()==0) {
    return;
  }

  // This method will be called again during initialize (it's ok, it's cheap)
  // But we need to have it called now, so we can use "get_field_out" below.
  set_fields_and_groups_pointers ();

  // Allow to request tendency of a field on a particular grid
  // by using the syntax 'field_name@grid_name'
  auto field_grid = [&] (const std::string& tn) -> std::pair<std::string,std::string>{
    auto tokens = ekat::split(tn,'@');
    EKAT_REQUIRE_MSG (tokens.size()==1 || tokens.size()==2,
        "Error! Invalid format for tendency calculation request: " + tn + "\n"
        "  To request tendencies for F, use 'F' or 'F@grid_name' format.\n");
    return std::make_pair(tokens[0],tokens.size()==2 ? tokens[1] : default_grid);
  };


  for (const auto& tn : tend_vec) {
    auto tokens = field_grid(tn);
    auto fn = tokens.first;
    auto gn = tokens.second;

    auto f = get_field_out(fn,gn);

    const auto& tname = this->name() + "_" + fn + "_tend";

    const auto fn_gn = fn + "@" + f.get_header().get_identifier().get_grid_name();

    // Create tend and start-of-step fields
    auto& tend = m_proc_tendencies[fn_gn] = f.clone(tname);
    m_start_of_step_fields[fn_gn] = f.clone();
    add_internal_field(tend,{"ACCUMULATED","DIVIDE_BY_DT"});
  }
}

void AtmosphereProcess::set_required_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_field(f.get_header().get_identifier()),
    "Error! Input field is not required by this atm process.\n"
    "    field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "    atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  m_inputs->add_field(f);

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as customer if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as customer to the field
    add_me_as_customer(f);
  }

  set_required_field_impl (f);

  add_py_fields(f);
}

void AtmosphereProcess::set_computed_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_field(f.get_header().get_identifier()),
    "Error! Input field is not computed by this atm process.\n"
   "   field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  m_outputs->add_field(f);

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as provider if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as provider for the field
    add_me_as_provider(f);
  }

  set_computed_field_impl (f);

  add_py_fields(f);
}

void AtmosphereProcess::set_required_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not require the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  m_inputs->add_group(group);
  
  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as customer if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    if (group.m_monolithic_field) {
      add_me_as_customer(*group.m_monolithic_field);
    } else {
      for (auto& it : group.m_individual_fields) {
        add_me_as_customer(*it.second);
      }
    }
  }

  set_required_group_impl(group);

  add_py_fields(group);
}

void AtmosphereProcess::set_computed_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not compute the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  m_outputs->add_group(group);
  
  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as provider if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    if (group.m_monolithic_field) {
      add_me_as_provider(*group.m_monolithic_field);
    } else {
      for (auto& it : group.m_individual_fields) {
        add_me_as_provider(*it.second);
      }
    }
  }

  set_computed_group_impl(group);

  add_py_fields(group);
}

void AtmosphereProcess::run_property_check (const prop_check_ptr&       property_check,
                                            const CheckFailHandling     check_fail_handling,
                                            const PropertyCheckCategory property_check_category) const {
  m_atm_logger->trace("[" + this->name() + "] run_property_check '" + property_check->name() + "'...");
  auto res_and_msg = property_check->check();

  // string for output
  std::string pre_post_str;
  if (property_check_category == PropertyCheckCategory::Precondition)  pre_post_str = "pre-condition";
  if (property_check_category == PropertyCheckCategory::Postcondition) pre_post_str = "post-condition";

  if (res_and_msg.result==CheckResult::Pass) {
    // Do nothing
  } else if (res_and_msg.result==CheckResult::Repairable) {
    // Ok, we can fix this
    property_check->repair();
    log (m_repair_log_level,
      "WARNING: Failed and repaired " + pre_post_str + " property check.\n"
      "  - Atmosphere process name: " + name() + "\n"
      "  - Property check name: " + property_check->name() + "\n"
      "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n"
      "  - Message: " + res_and_msg.msg + "\n");
  } else {
    // Ugh, the test failed badly, with no chance to repair it.
    if (check_fail_handling==CheckFailHandling::Warning) {
      // Still ok, just warn the user
      log (LogLevel::warn,
        "Warning! Failed " + pre_post_str + " property check (cannot be repaired).\n"
        "  - Atmosphere process name: " + name() + "\n"
        "  - Property check name: " + property_check->name() + "\n"
        "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n"
        "  - Message: " + res_and_msg.msg);
    } else {
      // No hope. Crash. But before crashing, print all the input and output fields at the samme
      // location where the check failed.
      std::ostringstream ss;
      ss << "Error! Failed " + pre_post_str + " property check (cannot be repaired).\n"
            "  - Atmosphere process name: " + name() + "\n"
            "  - Property check name: " + property_check->name() + "\n"
            "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n"
            "  - Message: " + res_and_msg.msg;
      if (res_and_msg.fail_loc_indices.size()>0) {
        // If the location is 3d, we not use the level index to slice the fields,
        // but print a column worth of data for all fields.
        using namespace ShortFieldTagsNames;
        auto tags = res_and_msg.fail_loc_tags;
        auto idx  = res_and_msg.fail_loc_indices;
        auto itm = ekat::find(tags,LEV);
        auto iti = ekat::find(tags,ILEV);
        if (itm!=tags.end()) {
          auto pos = std::distance(tags.begin(),itm);
          tags.erase(itm);
          idx.erase(idx.begin()+pos);
        } else if (iti!=tags.end()) {
          auto pos = std::distance(tags.begin(),iti);
          tags.erase(iti);
          idx.erase(idx.begin()+pos);
        }
        ss << "\n *************************** INPUT FIELDS ******************************\n";
        ss << "\n  ------- INPUT FIELDS -------\n";
        for (const auto& fname : m_inputs->field_names()) {
          auto f = m_inputs->get_field(fname);
          if (f.get_header().get_identifier().get_layout().has_tags(tags)) {
            print_field_hyperslab (f,tags,idx,ss);
            ss << " -----------------------------------------------------------------------\n";
          }
        }
        for (const auto& gname : m_inputs->group_names()) {
          auto g = m_inputs->get_field_group(gname);
          for (const auto& f : g.m_individual_fields) {
            if (f.second->get_header().get_identifier().get_layout().has_tags(tags)) {
              print_field_hyperslab (*f.second,tags,idx,ss);
              ss << " -----------------------------------------------------------------------\n";
            }
          }
        }
        ss << "\n ************************** OUTPUT FIELDS ******************************\n";
        for (const auto& fname : m_outputs->field_names()) {
          auto f = m_outputs->get_field(fname);
          if (f.get_header().get_identifier().get_layout().has_tags(tags)) {
            print_field_hyperslab (f,tags,idx,ss);
            ss << " -----------------------------------------------------------------------\n";
          }
        }
        for (const auto& gname : m_outputs->group_names()) {
          auto g = m_outputs->get_field_group(gname);
          for (const auto& f : g.m_individual_fields) {
            if (f.second->get_header().get_identifier().get_layout().has_tags(tags)) {
              print_field_hyperslab (*f.second,tags,idx,ss);
              ss << " -----------------------------------------------------------------------\n";
            }
          }
        }
      }
      EKAT_ERROR_MSG(ss.str());
    }
  }
}

void AtmosphereProcess::run_precondition_checks () const {
  m_atm_logger->debug("[" + this->name() + "] run_precondition_checks...");
  start_timer(m_timer_prefix + this->name() + "::run-precondition-checks");
  // Run all pre-condition property checks
  for (const auto& it : m_precondition_checks) {
    run_property_check(it.second, it.first,
                       PropertyCheckCategory::Precondition);
  }
  stop_timer(m_timer_prefix + this->name() + "::run-precondition-checks");
  m_atm_logger->debug("[" + this->name() + "] run_precondition_checks...done!");
}

void AtmosphereProcess::run_postcondition_checks () const {
  m_atm_logger->debug("[" + this->name() + "] run_postcondition_checks...");
  start_timer(m_timer_prefix + this->name() + "::run-postcondition-checks");
  // Run all post-condition property checks
  for (const auto& it : m_postcondition_checks) {
    run_property_check(it.second, it.first,
                       PropertyCheckCategory::Postcondition);
  }
  stop_timer(m_timer_prefix + this->name() + "::run-postcondition-checks");
  m_atm_logger->debug("[" + this->name() + "] run_postcondition_checks...done!");
}

void AtmosphereProcess::run_column_conservation_check () const {
  m_atm_logger->debug("[" + this->name() + "] run_column_conservation_check...");
  start_timer(m_timer_prefix + this->name() + "::run-column-conservation-checks");
  // Conservation check is run as a postcondition check
  run_property_check(m_conservation.second,
                     m_conservation.first,
                     PropertyCheckCategory::Postcondition);
  stop_timer(m_timer_prefix + this->name() + "::run-column-conservation-checks");
  m_atm_logger->debug("[" + this->name() + "] run_column-conservation_checks...done!");
}

void AtmosphereProcess::init_step_tendencies () {
  if (m_start_of_step_fields.size()==0) {
    return;
  }

  start_timer(m_timer_prefix + this->name() + "::compute_tendencies");
  for (auto& [fn_gn,f_beg] : m_start_of_step_fields) {
    const auto fname = ekat::split(fn_gn,"@")[0];
    const auto gname = ekat::split(fn_gn,"@")[1];
    const auto& f     = get_field_out(fname,gname);
    f_beg.deep_copy(f);
  }
  stop_timer(m_timer_prefix + this->name() + "::compute_tendencies");
}

void AtmosphereProcess::compute_step_tendencies () {
  if (m_start_of_step_fields.size()==0) {
    return;
  }

  m_atm_logger->debug("[" + this->name() + "] computing tendencies...");
  start_timer(m_timer_prefix + this->name() + "::compute_tendencies");
  for (auto& [fn_gn,tend] : m_proc_tendencies) {
    const auto fname = ekat::split(fn_gn,"@")[0];
    const auto gname = ekat::split(fn_gn,"@")[1];
    // Note: f_beg is nonconst, so we can store step tendency in it
          auto& f_beg = m_start_of_step_fields.at(fn_gn);
    const auto& f_end = get_field_out(fname,gname);

    // Compute tend from this atm proc step, then sum into overall atm timestep tendency
    // Note: don't add -f_beg to tend during init_step_tendencies, b/c the field magnitude
    // may be MUCH larger than the tend. Instead, compute step tend, and THEN add into tend
    f_beg.update(f_end,1,-1);
    tend.update(f_beg,1,1);
    tend.get_header().get_tracking().update_time_stamp(m_end_of_step_ts);
  }
  stop_timer(m_timer_prefix + this->name() + "::compute_tendencies");
}

bool AtmosphereProcess::has_required_field (const FieldIdentifier& id) const {
  return has_required_field(id.name(),id.get_grid_name());
}

bool AtmosphereProcess::has_required_field (const std::string& name, const std::string& grid_name) const
{
  for (const auto& r : m_field_requests) {
    if (r.fid.name()==name and r.fid.get_grid_name()==grid_name and r.usage & Required)
      return true;
  }
  return false;
}

bool AtmosphereProcess::has_computed_field (const FieldIdentifier& id) const {
  return has_computed_field(id.name(),id.get_grid_name());
}

bool AtmosphereProcess::has_computed_field (const std::string& name, const std::string& grid_name) const
{
  for (const auto& r : m_field_requests) {
    if (r.fid.name()==name and r.fid.get_grid_name()==grid_name and r.usage & Computed)
      return true;
  }
  return false;
}

bool AtmosphereProcess::has_required_group (const std::string& name, const std::string& grid) const
{
  for (const auto& r : m_group_requests) {
    if (r.name==name and r.grid==grid and r.usage & Required)
      return true;
  }
  return false;
}

bool AtmosphereProcess::has_computed_group (const std::string& name, const std::string& grid) const
{
  for (const auto& r : m_group_requests) {
    if (r.name==name and r.grid==grid and r.usage & Computed)
      return true;
  }
  return false;
}

void AtmosphereProcess::log (const LogLevel lev, const std::string& msg) const {
  m_atm_logger->log(lev,msg);
}

void AtmosphereProcess::set_update_time_stamps (const bool do_update) {
  m_update_time_stamps = do_update;
}

void AtmosphereProcess::update_time_stamps () {
  const auto& t = end_of_step_ts();

  // Update *all* output fields/groups, regardless of whether
  // they were touched at all during this time step.
  // TODO: this might have to be changed
  for (const auto& fname : m_outputs->field_names()) {
    auto f = m_outputs->get_field(fname);
    f.get_header().get_tracking().update_time_stamp(t);
  }
  for (const auto& gname : m_outputs->group_names()) {
    auto g = m_outputs->get_field_group(gname);
    if (g.m_monolithic_field) {
      g.m_monolithic_field->get_header().get_tracking().update_time_stamp(t);
    } else {
      for (auto& f : g.m_individual_fields) {
        f.second->get_header().get_tracking().update_time_stamp(t);
      }
    }
  }
}

void AtmosphereProcess::add_me_as_provider (const Field& f) {
  f.get_header_ptr()->get_tracking().add_provider(name());
}

void AtmosphereProcess::add_me_as_customer (const Field& f) {
  f.get_header_ptr()->get_tracking().add_customer(name());
}

void AtmosphereProcess::
add_internal_field (Field& f, const std::vector<std::string>& groups) {
  m_internals->add_field(f);
  for (const auto& gn : groups) {
    f.get_header().get_tracking().add_group(gn);
  }
  add_py_fields(f);
}

Field AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) const {
  return m_inputs->get_field(field_name, grid_name);
}

Field AtmosphereProcess::
get_field_in(const std::string& field_name) const {
  return m_inputs->get_field(field_name);
}

Field AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) const {
  return m_outputs->get_field(field_name, grid_name);
}

Field AtmosphereProcess::
get_field_out(const std::string& field_name) const {
  return m_outputs->get_field(field_name);
}

FieldGroup AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) const {
  return m_inputs->get_field_group(group_name, grid_name);
}

FieldGroup AtmosphereProcess::
get_group_in(const std::string& group_name) const {
  return m_inputs->get_field_group(group_name);
}

FieldGroup AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) const {
  return m_outputs->get_field_group(group_name, grid_name);
}

FieldGroup AtmosphereProcess::
get_group_out(const std::string& group_name) const {
  return m_outputs->get_field_group(group_name);
}

Field AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) const {
  return m_internals->get_field(field_name, grid_name);
}

Field AtmosphereProcess::
get_internal_field(const std::string& field_name) const {
  return m_internals->get_field(field_name);
}

const Field& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) const {
  return get_field_out_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) {
  return get_field_out_impl(field_name,grid_name);
}

const Field& AtmosphereProcess::
get_field_out(const std::string& field_name) const {
  return get_field_out_impl (field_name);
}

Field& AtmosphereProcess::
get_field_out(const std::string& field_name) {
  return get_field_out_impl (field_name);
}

const FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) const {
  return get_group_in_impl (group_name,grid_name);
}

FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) {
  return get_group_in_impl (group_name,grid_name);
}

const FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name) const {
  return get_group_in_impl(group_name);
}

FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name) {
  return get_group_in_impl(group_name);
}

const FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) const {
  return get_group_out_impl(group_name,grid_name);
}

FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) {
  return get_group_out_impl(group_name,grid_name);
}

const FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name) const {
  return get_group_out_impl(group_name);
}

FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name) {
  return get_group_out_impl(group_name);
}

const Field& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) const {
  return get_internal_field_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) {
  return get_internal_field_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_internal_field(const std::string& field_name) {
  return get_internal_field_impl(field_name);
}

const Field& AtmosphereProcess::
get_internal_field(const std::string& field_name) const {
  return get_internal_field_impl(field_name);
}

void AtmosphereProcess::
add_invariant_check (const prop_check_ptr& pc, const CheckFailHandling cfh)
{
  add_precondition_check (pc,cfh);
  add_postcondition_check (pc,cfh);
}

void AtmosphereProcess::
add_precondition_check (const prop_check_ptr& pc, const CheckFailHandling cfh)
{
  // If a pc can repair, we need to make sure the repairable
  // fields are among the computed fields of this atm proc.
  // Otherwise, it would be possible for this AP to implicitly
  // update a field, without that appearing in the dag.
  for (const auto& ptr : pc->repairable_fields()) {
    const auto& fid = ptr->get_header().get_identifier();
    EKAT_REQUIRE_MSG (
        has_computed_field(fid) || has_computed_group(fid.name(),fid.get_grid_name()),
        "Error! Input property check can repair a non-computed field.\n"
        "  - Atmosphere process name: " + name() + "\n"
        "  - Property check name: " + pc->name() + "\n");
  }
  m_precondition_checks.push_back(std::make_pair(cfh,pc));
}

void AtmosphereProcess::
add_postcondition_check (const prop_check_ptr& pc, const CheckFailHandling cfh)
{
  auto cfh2str = [] (const CheckFailHandling cfh) -> std::string {
    std::string s = "";
    switch (cfh) {
      case CheckFailHandling::Fatal:
        s = "fatal";
        break;
      case CheckFailHandling::Warning:
        s = "warning";
        break;
      default:
        EKAT_ERROR_MSG ("Unexpected/unsupported CheckFailHandling value.\n");
    }

    return "";
  };

  // Avoid adding the *SAME* test twice
  for (const auto& it : m_postcondition_checks) {
    if (it.second->same_as(*pc)) {
      EKAT_REQUIRE_MSG (it.first==cfh,
          "Error! Duplicate property check with different CheckFailHandling.\n"
          "  - Atmosphere process name: " + name() + "\n"
          "  - Property check name: " + pc->name() + "\n"
          "  - Current CFH: " + cfh2str(it.first) + "\n"
          "  - Input CFH: " + cfh2str(cfh) + "\n");
      return;
    }
  }

  // If a pc can repair, we need to make sure the repairable
  // fields are among the computed fields of this atm proc.
  // Otherwise, it would be possible for this AP to implicitly
  // update a field, without that appearing in the dag.
  for (const auto& ptr : pc->repairable_fields()) {
    const auto& fid = ptr->get_header().get_identifier();
    EKAT_REQUIRE_MSG (
        has_computed_field(fid) || has_computed_group(fid.name(),fid.get_grid_name()),
        "Error! Input property check can repair a non-computed field.\n"
        "  - Atmosphere process name: " + name() + "\n"
        "  - Property check name: " + pc->name() + "\n");
  }
  m_postcondition_checks.push_back(std::make_pair(cfh,pc));
}

void AtmosphereProcess::
add_column_conservation_check(const prop_check_ptr &prop_check, const CheckFailHandling cfh)
{
  EKAT_REQUIRE_MSG(m_conservation.second == nullptr,
                   "Error! Conservation check for process \""+ name() +
                   "\" has already been added.");

  m_conservation = std::make_pair(cfh,prop_check);
}

void AtmosphereProcess::
alias_field_in (const std::string& field_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  auto f = m_inputs->get_field(field_name, grid_name);
  auto aliased = f.alias(alias_name);
  m_inputs->add_field(aliased);
}

void AtmosphereProcess::
alias_field_out (const std::string& field_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  auto f = m_outputs->get_field(field_name, grid_name);
  auto aliased = f.alias(alias_name);
  m_outputs->add_field(aliased);
}

void AtmosphereProcess::
alias_group_in (const std::string& group_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  auto g = m_inputs->get_field_group(group_name, grid_name);
  auto aliased = g.alias(alias_name);
  m_inputs->add_group(aliased);
}

void AtmosphereProcess::
alias_group_out (const std::string& group_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  auto g = m_outputs->get_field_group(group_name, grid_name);
  auto aliased = g.alias(alias_name);
  m_outputs->add_group(aliased);
}

void AtmosphereProcess::add_internal_field (Field& f, const std::vector<std::string>& groups) {
  m_internals->add_field(f);
  for (const auto& g : groups) {
    if (m_internals->has_group(g)) {
      auto fg = m_internals->get_field_group(g);
      fg.m_info->m_fields_names.insert(f.name());
    } else {
      auto fg = FieldGroup(g, f.get_header().get_identifier().get_grid_name());
      fg.m_info->m_fields_names.insert(f.name());
      fg.m_info->m_bundled = false;
      m_internals->add_group(fg);
    }
  }
}

void AtmosphereProcess::compute_column_conservation_checks_data (const double dt)
{
  EKAT_REQUIRE_MSG(m_conservation.second != nullptr,
                   "Error! User set enable_column_conservation_checks=true, "
                   "or has_energy_fixer=true, "
                   "but no conservation check class exists.\n");

  // Set dt and compute current mass and energy.
  const auto& conservation_check =
      std::dynamic_pointer_cast<MassAndEnergyConservationCheck>(m_conservation.second);
  conservation_check->set_dt(dt);
  conservation_check->compute_current_mass();
  conservation_check->compute_current_energy();
}

void AtmosphereProcess::fix_energy (const double dt, const bool water_thermo_fixer, const bool print_debug_info)
{
  EKAT_REQUIRE_MSG(m_conservation.second != nullptr,
                   "Error! User set has_energy_fixer=true, "
                   "but no conservation check class exists.\n");

  // Set dt and compute current mass and energy.
  const auto& conservation_check =
      std::dynamic_pointer_cast<MassAndEnergyConservationCheck>(m_conservation.second);

  //dt is needed to convert flux to change
  conservation_check->set_dt(dt);
  conservation_check->global_fixer(water_thermo_fixer, print_debug_info);

  if(print_debug_info){
    //print everything about the fixer only in debug mode
    m_atm_logger->info("EAMxx:: energy fixer: T tend added to each physics midlevel " + std::to_string( conservation_check->get_pb_fixer() ) + " K" );
    m_atm_logger->info("EAMxx:: energy fixer: total energy before fix " + std::to_string( conservation_check->get_total_energy_before() ) + " J");
    std::stringstream ss;
    ss << "EAMxx:: energy fixer: rel energy error after fix " << std::setprecision(15) << conservation_check->get_echeck() << "\n";
    m_atm_logger->info(ss.str());
  }
}

void AtmosphereProcess::add_py_fields (const Field& f)
{
#ifdef EAMXX_HAS_PYTHON
  if (m_py_module.has_value()) {
    const auto& grid_name = f.get_header().get_identifier().get_grid_name();
    m_py_fields_dev[f.name()][grid_name] = create_py_field<Device>(f);
    m_py_fields_host[f.name()][grid_name] = create_py_field<Host>(f);
  }
#else
  (void) f;
#endif
}

void AtmosphereProcess::add_py_fields (const FieldGroup& group)
{
#ifdef EAMXX_HAS_PYTHON
  if (m_py_module.has_value()) {
    if (group.m_monolithic_field) {
      const auto& f = group.m_monolithic_field;
      add_py_fields(*f);
    } else {
      for (const auto& [name,f] : group.m_individual_fields) {
        add_py_fields(*f);
      }
    }
  }
#else
  (void) group;
#endif
}

} // namespace scream
