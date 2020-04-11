#ifndef SCREAM_ATMOSPHERE_PROCESS_DAG_HPP
#define SCREAM_ATMOSPHERE_PROCESS_DAG_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"

namespace scream {

class AtmProcDAG {
public:

  using group_type = AtmosphereProcessGroup;

  void create_dag (const group_type& atm_procs);

  void write_dag (const std::string& fname) const;

  const std::set<std::string>& get_prev_time_step_dependencies () const {
    return m_deps_on_prev_ts_set;
  }

  const std::set<std::string>& get_unmet_dependencies () const {
    return m_unmet_deps_set;
  }

protected:

  void cleanup ();

  void add_nodes (const group_type& atm_procs);

  struct Node {
    std::vector<int>          parents;
    std::vector<int>          children;
    std::string               name;
    int                       id;
    std::vector<std::string>  computed;
    std::vector<std::string>  required;
  };

  // Map each field id to its last provider
  std::map<std::string,int>               m_fid_to_last_provider;

  // Map a node id to a vector of unmet dependencies
  std::map<int,std::vector<std::string>>  m_unmet_deps;

  // A lumped set version of the above, which can be queried to do checks on the AD side.
  std::set<std::string>                   m_unmet_deps_set;

  // Map a node id to a vector of deps met from previous time-step
  std::map<int,std::vector<std::string>>  m_deps_on_prev_ts;

  // A lumped set version of the above, which can be queried to do checks on the AD side.
  std::set<std::string>                   m_deps_on_prev_ts_set;

  using Edge = std::pair<int,int>;

  std::vector<Node>   m_nodes;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_DAG_HPP
