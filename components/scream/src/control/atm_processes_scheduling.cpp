#include "atm_processes_scheduling.hpp"

#include "share/error_defs.hpp"
#include "share/util/string_utils.hpp"

namespace scream {

void ProcessesSchedule::check_unique_processes () const {
  std::set<std::string> procs_names;
  check_unique_processes(m_procs_graph,procs_names);
}

bool ProcessesSchedule::check_process_present (const std::string& proc_name) const {
  return check_process_present(m_procs_graph,proc_name);
}

std::vector<std::pair<std::string,Comm>>
ProcessesSchedule::get_processes_pipeline(const Comm& atm_comm) const {
  return get_processes_pipeline (m_procs_graph,atm_comm);
}

void ProcessesSchedule::check_unique_processes (const ProcessList& list, std::set<std::string>& existing_procs) const {
  for (const auto& entry : list.m_list) {
    if (entry.m_is_sublist) {
      // Check the sublist, adding its processes to the set
      check_unique_processes (*entry.sublist, existing_procs);
    } else {
      // Try to insert this process name. If insertion fails, process was already present.
      auto it_bool = existing_procs.insert(*entry.entry);
      error::runtime_check(it_bool.second,"Error! Repeated entry in the processes graph.\n",-1);
    }
  }
}

bool ProcessesSchedule::check_process_present (const ProcessList& list, const std::string& proc_name) const {
  bool present = false;
  std::string proc_name_upper = util::upper_case(proc_name);
  for (const auto& entry : list.m_list) {
    if (entry.m_is_sublist) {
      // Check if present in the sublist
      present |= check_process_present(*entry.sublist,proc_name);
    } else {
      // Check if this process if the desired one
      present |= (util::upper_case(*entry.entry) == proc_name_upper);
    }
  }

  return present;
}

std::vector<std::pair<std::string,Comm>>
ProcessesSchedule::get_processes_pipeline(const ProcessList& list, const Comm& comm) const {
  std::vector<std::pair<std::string,Comm>> procs_names_and_comms;
  if (list.m_schedule==ScheduleType::Parallel && comm.size()>1) {
    error::runtime_abort("Errror! Parallel processes execution not yet implemented.\n",-1);
  } else {
    // Simply grab all the entries
    for (const auto& entry : list.m_list) {
      if (entry.m_is_sublist) {
        // Get the processes in the sublist...
        auto&& tmp = get_processes_pipeline(*entry.sublist,comm);

        // ...and add them to the result
        std::move(tmp.begin(),tmp.end(),std::back_inserter(procs_names_and_comms));
      } else {
        // Simply add this process, with the input comm
        procs_names_and_comms.emplace_back(*entry.entry,comm);
      }
    }
  }

  return procs_names_and_comms;
}

} // namespace scream
