#include "share/atm_process/atmosphere_process.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/property_checks/mass_and_energy_column_conservation_check.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_assert.hpp"

#include <set>
#include <stdexcept>
#include <string>

namespace scream
{

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
}

AtmosphereProcess::
AtmosphereProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
  : m_comm       (comm)
  , m_params     (params)
{
  if (m_params.isParameter("Logger")) {
    m_atm_logger = m_params.get<std::shared_ptr<logger_t>>("Logger");
  } else {
    // Create a console-only logger, that logs all ranks
    using namespace ekat::logger;
    using logger_impl_t = Logger<LogNoFile,LogAllRanks>;
    auto log_level = m_params.get<std::string>("log_level","trace");
    m_atm_logger = std::make_shared<logger_impl_t>("",str2LogLevel(log_level),m_comm);
    m_atm_logger->set_no_format();
  }

  if (m_params.isParameter("number_of_subcycles")) {
    m_num_subcycles = m_params.get<int>("number_of_subcycles");
  }
  EKAT_REQUIRE_MSG (m_num_subcycles>0,
      "Error! Invalid number of subcycles in param list " + m_params.name() + ".\n"
      "  - Num subcycles: " + std::to_string(m_num_subcycles) + "\n");

  m_timer_prefix = m_params.get<std::string>("Timer Prefix","EAMxx::");

  m_repair_log_level = str2LogLevel(m_params.get<std::string>("repair_log_level","warn"));

  // Info for mass and energy conservation checks
  m_column_conservation_check_data.has_check =
      m_params.get<bool>("enable_column_conservation_checks", false);

  m_internal_diagnostics_level = m_params.get<int>("internal_diagnostics_level", 0);
}

void AtmosphereProcess::initialize (const TimeStamp& t0, const RunType run_type) {
  if (this->type()!=AtmosphereProcessType::Group) {
    start_timer (m_timer_prefix + this->name() + "::init");
  }

  log (LogLevel::info,"  Initializing " + name() + "...");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

  set_fields_and_groups_pointers();
  m_start_of_step_ts = m_end_of_step_ts = t0;
  initialize_impl(run_type);

  // Create all start-of-step fields needed for tendencies calculation
  for (const auto& it : m_proc_tendencies) {
    const auto& tname = it.first;
    const auto& fname = m_tend_to_field.at(tname);
    m_start_of_step_fields[fname] = get_field_out(fname).clone();
  }

  log (LogLevel::info,"  Initializing " + name() + "... done!");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

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

    if (has_column_conservation_check()) {
      // Column local mass and energy checks requires the total mass and energy
      // to be computed directly before the atm process is run, as well and store
      // the correct timestep for the process.
      compute_column_conservation_checks_data(dt_sub);
    }

    if (m_internal_diagnostics_level > 0)
      // Print hash of INPUTS before run
      print_global_state_hash(name() + "-pre-sc-" + std::to_string(m_subcycle_iter),
                              m_start_of_step_ts,
                              true, false, false);

    // Run derived class implementation
    run_impl(dt_sub);

    if (m_internal_diagnostics_level > 0)
      // Print hash of OUTPUTS/INTERNALS after run
      print_global_state_hash(name() + "-pst-sc-" + std::to_string(m_subcycle_iter),
                              m_end_of_step_ts,
                              true, true, true);

    if (has_column_conservation_check()) {
      // Run the column local mass and energy conservation checks
      run_column_conservation_check();
    }
  }

  // Complete tendency calculations (if any)
  compute_step_tendencies(dt);

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

void AtmosphereProcess::finalize (/* what inputs? */) {
  finalize_impl(/* what inputs? */);
}

void AtmosphereProcess::setup_tendencies_requests () {
  using vos_t = std::vector<std::string>;
  auto tend_vec = m_params.get<vos_t>("compute_tendencies",{});
  if (tend_vec.size()==0) {
    return;
  }

  // Will convert to list, since lists insert/erase are safe inside loops
  std::list<std::string> tend_list;

  // Allow user to specify [all] as a shortcut for all updated fields
  if (tend_vec.size()==1 && tend_vec[0]==ci_string("all")) {
    for (const auto& r_out : get_computed_field_requests()) {
      for (const auto& r_in : get_required_field_requests()) {
        if (r_in.fid==r_out.fid) {
          tend_list.push_back(r_out.fid.name());
          break;
        }
      }
    }
  } else {
    for (const auto& n : tend_vec) {
      // We don't have grid name info here, so we'll check that n
      // is indeed an updated field a few lines below
      tend_list.push_back(n);
    }
  }

  // Allow to request tendency of a field on a particular grid
  // by using the syntax 'field_name#grid_name'
  auto field_grid = [&] (const std::string& tn) -> std::pair<std::string,std::string>{
    auto tokens = ekat::split(tn,'#');
    EKAT_REQUIRE_MSG (tokens.size()==1 || tokens.size()==2,
        "Error! Invalid format for tendency calculation request: " + tn + "\n"
        "  To request tendencies for F, use 'F' or 'F#grid_name' format.\n");
    return std::make_pair(tokens[0],tokens.size()==2 ? tokens[1] : "");
  };

  auto multiple_matches_error = [&] (const std::string& tn) {
    auto tokens = field_grid(tn);
    auto fn = tokens.first;
    auto gn = tokens.second;

    std::set<std::pair<std::string,std::string>> matches_set;
    std::vector<std::string> matches_vec;
    for (auto it : get_computed_field_requests()) {
      if (it.fid.name()==fn && (gn=="" || gn==it.fid.get_grid_name())) {
        auto it_bool = matches_set.insert(std::make_pair(fn,it.fid.get_grid_name()));
        if (it_bool.second) {
          auto gname = it_bool.first->second;
          matches_vec.push_back(fn + "(grid=" + gname + ")");
        }
      }
    }
    EKAT_ERROR_MSG(
        "Error! Could not uniquely determine field for tendency calculation.\n"
        " - atm proc name: " + this->name() + "\n"
        " - tend requested: " + tn + "\n"
        " - all matches: " + ekat::join(matches_vec,", ") + "\n"
        "If atm process computes " + fn + " on multiple grids, use " + fn + "#grid_name\n");
  };

  // Parse all the tendencies requested, looking for a corresponding
  // field with the requested name. If more than one is found (must be
  // that we have same field on multiple grids), error out.
  using namespace ekat::units;
  using strlist_t = std::list<std::string>;
  for (const auto& tn : tend_list) {
    std::string grid_found = "";
    auto tokens = field_grid(tn);
    auto fn = tokens.first;
    auto gn = tokens.second;
    const auto& tname = this->name() + "_" + tn + "_tend";
    for (auto it : get_computed_field_requests()) {
      const auto& fid = it.fid;
      if (fid.name()==fn && (gn=="" || gn==fid.get_grid_name())) {
        // Get fid specs
        const auto& layout = fid.get_layout();
        const auto  units  = fid.get_units() / s;
        const auto& gname  = fid.get_grid_name();
        const auto& dtype  = fid.data_type();

        // Check we did not process another field with same name
        if (grid_found!="" && gname!=grid_found) {
          multiple_matches_error(tn);
        }

        // Check this computed field is also required (hence, updated)
        EKAT_REQUIRE_MSG (has_required_field(fn,gname),
            "Error! Tendency calculation available only for updated fields.\n"
            " - atm proc name: " + this->name() + "\n"
            " - field name   : " + fn + "\n"
            " - grid name    : " + gn + "\n");

        // Create tend FID and request field
        FieldIdentifier t_fid(tname,layout,units,gname,dtype);
        add_field<Computed>(t_fid,strlist_t{"ACCUMULATED","DIVIDE_BY_DT"});
        grid_found = gname;
      }
    }

    EKAT_REQUIRE_MSG (grid_found!="",
        "Error! Could not locate field for tendency request.\n"
        " - atm proc name: " + this->name() + "\n"
        " - tendency request: " + tn + "\n");

    // Add an empty field, to be reset when we set the computed field
    m_proc_tendencies[tname] = {};
    m_tend_to_field[tname] = fn;
  }

  m_compute_proc_tendencies = m_proc_tendencies.size()>0;
}

void AtmosphereProcess::set_required_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_field(f.get_header().get_identifier()),
    "Error! Input field is not required by this atm process.\n"
    "    field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "    atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_fields_in,f)) {
    m_fields_in.emplace_back(f);
  }

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as customer if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as customer to the field
    add_me_as_customer(f);
  }

  set_required_field_impl (f);
}

void AtmosphereProcess::set_computed_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_field(f.get_header().get_identifier()),
    "Error! Input field is not computed by this atm process.\n"
    "   field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_fields_out,f)) {
    m_fields_out.emplace_back(f);
  }

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as provider if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as provider for the field
    add_me_as_provider(f);
  }

  set_computed_field_impl (f);

  if (m_compute_proc_tendencies && m_proc_tendencies.count(f.name())==1) {
    m_proc_tendencies[f.name()] = f;
  }
}

void AtmosphereProcess::set_required_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not require the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_groups_in,group)) {
    m_groups_in.emplace_back(group);
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
  }

  set_required_group_impl(group);
}

void AtmosphereProcess::set_computed_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not compute the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_groups_out,group)) {
    m_groups_out.emplace_back(group);
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
  }

  set_computed_group_impl(group);
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
        for (const auto& f : m_fields_in) {
          if (f.get_header().get_identifier().get_layout().has_tags(tags)) {
            print_field_hyperslab (f,tags,idx,ss);
            ss << " -----------------------------------------------------------------------\n";
          }
        }
        for (const auto& g : m_groups_in) {
          for (const auto& f : g.m_individual_fields) {
            if (f.second->get_header().get_identifier().get_layout().has_tags(tags)) {
              print_field_hyperslab (*f.second,tags,idx,ss);
              ss << " -----------------------------------------------------------------------\n";
            }
          }
        }
        ss << "\n ************************** OUTPUT FIELDS ******************************\n";
        for (const auto& f : m_fields_out) {
          if (f.get_header().get_identifier().get_layout().has_tags(tags)) {
            print_field_hyperslab (f,tags,idx,ss);
            ss << " -----------------------------------------------------------------------\n";
          }
        }
        for (const auto& g : m_groups_out) {
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
  run_property_check(m_column_conservation_check.second,
                     m_column_conservation_check.first,
                     PropertyCheckCategory::Postcondition);
  stop_timer(m_timer_prefix + this->name() + "::run-column-conservation-checks");
  m_atm_logger->debug("[" + this->name() + "] run_column-conservation_checks...done!");
}

void AtmosphereProcess::init_step_tendencies () {
  if (m_compute_proc_tendencies) {
    start_timer(m_timer_prefix + this->name() + "::compute_tendencies");
    for (auto& it : m_start_of_step_fields) {
      const auto& fname = it.first;
      const auto& f     = get_field_out(fname);
            auto& f_beg = it.second;
      f_beg.deep_copy(f);
    }
    stop_timer(m_timer_prefix + this->name() + "::compute_tendencies");
  }
}

void AtmosphereProcess::compute_step_tendencies (const double dt) {
  using namespace ShortFieldTagsNames;
  if (m_compute_proc_tendencies) {
    m_atm_logger->debug("[" + this->name() + "] computing tendencies...");
    start_timer(m_timer_prefix + this->name() + "::compute_tendencies");
    for (auto it : m_proc_tendencies) {
      // Note: f_beg is nonconst, so we can store step tendency in it
      const auto& tname = it.first;
      const auto& fname = m_tend_to_field.at(tname);
      const auto& f     = get_field_out(fname);
            auto& f_beg = m_start_of_step_fields.at(fname);
            auto& tend  = it.second;

      // Compute tend from this atm proc step, then sum into overall atm timestep tendency
      f_beg.update(f,1,-1);
      tend.update(f_beg,1,1);
    }
    stop_timer(m_timer_prefix + this->name() + "::compute_tendencies");
  }
}

bool AtmosphereProcess::has_required_field (const FieldIdentifier& id) const {
  return has_required_field(id.name(),id.get_grid_name());
}

bool AtmosphereProcess::has_required_field (const std::string& name, const std::string& grid_name) const {
  for (const auto& it : m_required_field_requests) {
    if (it.fid.name()==name && it.fid.get_grid_name()==grid_name) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::has_computed_field (const FieldIdentifier& id) const {
  return has_computed_field(id.name(),id.get_grid_name());
}

bool AtmosphereProcess::has_computed_field (const std::string& name, const std::string& grid_name) const {
  for (const auto& it : m_computed_field_requests) {
    if (it.fid.name()==name && it.fid.get_grid_name()==grid_name) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::has_required_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_required_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::has_computed_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_computed_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
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
  for (auto& f : m_fields_out) {
    f.get_header().get_tracking().update_time_stamp(t);
  }
  for (auto& g : m_groups_out) {
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
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

void AtmosphereProcess::add_me_as_customer (const Field& f) {
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void AtmosphereProcess::add_internal_field (const Field& f) {
  m_internal_fields.push_back(f);
}

const Field& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) const {
  return get_field_in_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) {
  return get_field_in_impl(field_name,grid_name);
}

const Field& AtmosphereProcess::
get_field_in(const std::string& field_name) const {
  return get_field_in_impl(field_name);
}

Field& AtmosphereProcess::
get_field_in(const std::string& field_name) {
  return get_field_in_impl(field_name);
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
        s = "Fatal";
        break;
      case CheckFailHandling::Warning:
        s = "Warning";
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
  EKAT_REQUIRE_MSG(m_column_conservation_check.second == nullptr,
                   "Error! Conservation check for process \""+ name() +
                   "\" has already been added.");

  m_column_conservation_check = std::make_pair(cfh,prop_check);
}

void AtmosphereProcess::set_fields_and_groups_pointers () {
  for (auto& f : m_fields_in) {
    const auto& fid = f.get_header().get_identifier();
    m_fields_in_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
  for (auto& f : m_fields_out) {
    const auto& fid = f.get_header().get_identifier();
    m_fields_out_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
  for (auto& g : m_groups_in) {
    const auto& group_name = g.m_info->m_group_name;
    m_groups_in_pointers[group_name][g.grid_name()] = &g;
  }
  for (auto& g : m_groups_out) {
    const auto& group_name = g.m_info->m_group_name;
    m_groups_out_pointers[group_name][g.grid_name()] = &g;
  }
  for (auto& f : m_internal_fields) {
    const auto& fid = f.get_header().get_identifier();
    m_internal_fields_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
}

void AtmosphereProcess::
alias_field_in (const std::string& field_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  try {
    m_fields_in_pointers[alias_name][grid_name] = m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field for aliasing request.\n"
        "    - field name: " + field_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_field_out (const std::string& field_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  try {
    m_fields_out_pointers[alias_name][grid_name] = m_fields_out_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field for aliasing request.\n"
        "    - field name: " + field_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_group_in (const std::string& group_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  try {
    m_groups_in_pointers[alias_name][grid_name] = m_groups_in_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group for aliasing request.\n"
        "    - group name: " + group_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_group_out (const std::string& group_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  try {
    m_groups_out_pointers[alias_name][grid_name] = m_groups_out_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group for aliasing request.\n"
        "    - group name: " + group_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_in_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
  static Field f;
  return f;
}

Field& AtmosphereProcess::
get_field_in_impl(const std::string& field_name) const {
  try {
    auto& copies = m_fields_in_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find input field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
  static Field f;
  return f;
}

Field& AtmosphereProcess::
get_field_out_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_fields_out_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
  static Field f;
  return f;
}

Field& AtmosphereProcess::
get_field_out_impl(const std::string& field_name) const {
  try {
    auto& copies = m_fields_out_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find output field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
  static Field f;
  return f;
}

FieldGroup& AtmosphereProcess::
get_group_in_impl(const std::string& group_name, const std::string& grid_name) const {
  try {
    return *m_groups_in_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
  static FieldGroup g("");
  return g;
}

FieldGroup& AtmosphereProcess::
get_group_in_impl(const std::string& group_name) const {
  try {
    auto& copies = m_groups_in_pointers.at(group_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find input group providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  group name: " + group_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n");
  }
  static FieldGroup g("");
  return g;
}

FieldGroup& AtmosphereProcess::
get_group_out_impl(const std::string& group_name, const std::string& grid_name) const {
  try {
    return *m_groups_out_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
  static FieldGroup g("");
  return g;
}

FieldGroup& AtmosphereProcess::
get_group_out_impl(const std::string& group_name) const {
  try {
    auto& copies = m_groups_out_pointers.at(group_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find output group providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  group name: " + group_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n");
  }
  static FieldGroup g("");
  return g;
}

Field& AtmosphereProcess::
get_internal_field_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_internal_fields_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate internal field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
  static Field f;
  return f;
}

Field& AtmosphereProcess::
get_internal_field_impl(const std::string& field_name) const {
  try {
    auto& copies = m_internal_fields_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find internal field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate internal field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
  static Field f;
  return f;
}

void AtmosphereProcess
::remove_field (const std::string& field_name, const std::string& grid_name) {
  typedef std::list<Field>::iterator It;
  const auto rmf = [&] (std::list<Field>& fields, strmap_t<strmap_t<Field*>>& ptrs) {
    std::vector<It> rm_its;
    for (It it = fields.begin(); it != fields.end(); ++it) {
      const auto& fid = it->get_header().get_identifier();
      if (fid.name() == field_name and fid.get_grid_name() == grid_name) {
        rm_its.push_back(it);
        ptrs[field_name][grid_name] = nullptr;
      }
    }
    for (auto& it : rm_its) fields.erase(it);
  };
  rmf(m_fields_in, m_fields_in_pointers);
  rmf(m_fields_out, m_fields_out_pointers);
  rmf(m_internal_fields, m_internal_fields_pointers);
}

void AtmosphereProcess
::remove_group (const std::string& group_name, const std::string& grid_name) {
  typedef std::list<FieldGroup>::iterator It;
  const auto rmg = [&] (std::list<FieldGroup>& fields, strmap_t<strmap_t<FieldGroup*>>& ptrs) {
    std::vector<It> rm_its;
    for (It it = fields.begin(); it != fields.end(); ++it) {
      if (it->m_info->m_group_name == group_name and it->grid_name() == grid_name) {
        rm_its.push_back(it);
        ptrs[group_name][grid_name] = nullptr;
        for (auto& kv : it->m_individual_fields)
          remove_field(kv.first, grid_name);
      }
    }
    for (auto& it : rm_its) fields.erase(it);
  };
  rmg(m_groups_in, m_groups_in_pointers);
  rmg(m_groups_out, m_groups_out_pointers);
}

void AtmosphereProcess::compute_column_conservation_checks_data (const int dt)
{
  EKAT_REQUIRE_MSG(m_column_conservation_check.second != nullptr,
                   "Error! User set enable_column_conservation_checks=true, "
                   "but no conservation check exists.\n");

  // Set dt and compute current mass and energy.
  const auto& conservation_check =
      std::dynamic_pointer_cast<MassAndEnergyColumnConservationCheck>(m_column_conservation_check.second);
  conservation_check->set_dt(dt);
  conservation_check->compute_current_mass();
  conservation_check->compute_current_energy();
}

} // namespace scream
