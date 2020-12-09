#ifndef SCREAM_ATMOSPHERE_PROCESS_DAG_HPP
#define SCREAM_ATMOSPHERE_PROCESS_DAG_HPP

#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/field/field_initializer.hpp"

namespace scream {

class AtmProcDAG {
public:

  using group_type = AtmosphereProcessGroup;
  static constexpr int VERB_MAX = 4;

  void create_dag (const group_type& atm_procs,
                   const std::shared_ptr<FieldRepository<Real>> field_repo);

  void add_field_initializer (const FieldInitializer& initializer);

  void write_dag (const std::string& fname, const int verbosity = VERB_MAX) const;

  bool has_unmet_dependencies () const { return m_has_unmet_deps; }
  const std::map<int,std::set<int>>& unmet_deps () const {
    return m_unmet_deps;
  }

protected:

  void cleanup ();

  void add_nodes (const group_type& atm_procs,
                  const std::shared_ptr<FieldRepository<Real>> field_repo);

  // Add fid to list of fields in the dag, and return its position.
  // If already stored, simply return its position
  int add_fid (const FieldIdentifier& fid);

  void update_unmet_deps ();

  struct Node {
    std::vector<int>  children;
    std::string       name;
    int               id;
    std::set<int>     computed;
    std::set<int>     required;
  };

  // Assign an id to each field identifier
  std::vector<FieldIdentifier>    m_fids;

  // Map each field id to its last provider
  std::map<int,int>               m_fid_to_last_provider;

  // Map a node id to a set of unmet field dependencies
  std::map<int,std::set<int>>     m_unmet_deps;
  bool                            m_has_unmet_deps;

  // The nodes in the atm DAG
  std::vector<Node>               m_nodes;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_DAG_HPP
