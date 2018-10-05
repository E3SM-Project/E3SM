#ifndef SCREAM_ATM_PROCESSES_SCHEDULING_HPP
#define SCREAM_ATM_PROCESSES_SCHEDULING_HPP

#include "share/mpi/scream_comm.hpp"

#include <vector>   // For std::vector
#include <set>      // For std::set
#include <memory>   // For std::shared_ptr
#include <utility>  // For std::pair

namespace scream {

/*
 *  Data structure needed to create a schedule for the execution of atm processes.
 *
 *  We want to allow (almost) arbitrary scheduling of atmosphere processes.
 *  In particular, processes can run sequentially, in parallel, or in any
 *  combination of the two.
 *  The following diagram illustrate an example of atmosphere processes splitting:
 *  
 *       ---- Dynamics ---
 *      /                 \
 *   -- ---- Radiation ---  --- SurfaceCoupling
 *      \                 /
 *       --- P3 -- SHOC --
 *  
 *  The diagram imposes left-to-right ordering. In particular, SurfaceCoupling
 *  will run after everything else has completed. P3 and SHOC will run sequentially,
 *  but in parallel with Dynamics and Radiation.
 *  To allow this, we need a data structure that can store this kind of diagram
 *  information, and then, for each rank, retrieve the list of atm processes it
 *  is responsible to run. E.g., with 3 ranks, we may have the following pipelines:
 *  
 *   - rank 0: Dynamics, SurfaceCoupling
 *   - rank 1: Radiation, SurfaceCoupling
 *   - rank 2: P3, SHOC, SurfaceCoupling
 *  
 *  In this case, the SurfaceCoupling process will be split among 3 ranks, while
 *  all the other processes will be executed on a single rank.
 */

// Enum to specify whether processes in a list have to be run sequentially or in parallel.
enum class ScheduleType {
  Sequential,
  Parallel
};

// Fwd declaration, needed to break circular dependency between ProcessList and ProcessListEntry,
// since ProcessListEntry could store a ProcessList, and a ProcessList stores ProcessListEntry's.
struct ProcessList;

// We want to create lists that allow recursion. So an entry in a list can be either the name
// of a process or a sublist (possibly with sublists itself). The ideal solution would be
// std::any, but it's in the c++17 standard, which we won't have on some machines (e.g.,Summit).
// An entry of the list is morally a union (proc or list), and a bool, to help
// us see what it is. You may think we could use a 'union' (in the sense of the
// c++ keyword), but that's not possible, since std::shared_ptr has a non-trivial
// destructor.
struct ProcessListEntry {
  std::shared_ptr<std::string>  entry;
  std::shared_ptr<ProcessList>  sublist;
  bool                          m_is_sublist;
};

// Finally, a list of ProcessListEntries
struct ProcessList {
  ScheduleType                    m_schedule;
  std::vector<ProcessListEntry>   m_list;

  // This is needed only if the procs in the list are to be run in parallel
  // The integers should be the number of ranks to be assigned to each
  // entry in the list.
  std::vector<int>                m_resources;
};

// A class that stores meta-data representing a schedule.
class ProcessesSchedule {
public:

  void check_unique_processes () const; 
  bool check_process_present (const std::string& proc_name) const; 
  std::vector<std::pair<std::string,Comm>> get_processes_pipeline (const Comm& atm_comm) const;

private:

  void check_unique_processes (const ProcessList& list, std::set<std::string>& existing_procs) const;
  bool check_process_present (const ProcessList& list, const std::string& proc_name) const; 
  std::vector<std::pair<std::string,Comm>> get_processes_pipeline (const ProcessList& list, const Comm& comm) const;

  ProcessList     m_procs_graph;
};

} // namespace scream

#endif // SCREAM_ATM_PROCESSES_SCHEDULING_HPP
