#ifndef SCREAM_ATMOSPHERE_PROCESS_DAG_HPP
#define SCREAM_ATMOSPHERE_PROCESS_DAG_HPP

#include <memory>
#include <string>
#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/field/field_group.hpp"

namespace scream {

class AtmProcDAG {
public:

  using group_type = AtmosphereProcessGroup;
  static constexpr int VERB_MAX = 4;

  void create_dag (const group_type& atm_procs);

  using grid_field_map = std::map<std::string,std::vector<std::string>>;
  void process_initial_conditions(const grid_field_map &ic_inited);

  void init_atm_proc_nodes(const group_type& atm_procs);

  void add_surface_coupling (const std::set<FieldIdentifier>& imports,
                             const std::set<FieldIdentifier>& exports);

  void write_dag (const std::string& fname, const int verbosity = VERB_MAX) const;

  bool has_unmet_dependencies () const { return m_has_unmet_deps; }
  const std::map<int,std::set<int>>& unmet_deps () const {
    return m_unmet_deps;
  }

protected:

  void cleanup ();

  void add_nodes (const group_type& atm_procs);

  void add_edges ();

  // Add fid to list of fields in the dag, and return its position.
  // If already stored, simply return its position
  int add_fid (const FieldIdentifier& fid);

  // Internally, we store FID's in a vector. This method find the index
  // of the input FID in said vector.
  int get_fid_index (const FieldIdentifier& fid) const;

  void update_unmet_deps ();

  struct Node {
    std::vector<int>  children;
    std::string       name;
    int               id;
    std::set<int>     computed;     // output fields
    std::set<int>     required;     // input  fields
    std::set<int>     gr_computed;  // output groups
    std::set<int>     gr_required;  // input  groups
  };

  // Assign an id to each field identifier
  std::vector<FieldIdentifier>            m_fids;

  // Store groups so we can print info of their members if need be
  std::map<FieldIdentifier,FieldGroup>    m_gr_fid_to_group;

  // Map each field id to its last provider
  std::map<int,int>               m_fid_to_last_provider;

  // Map a node id to a set of unmet field dependencies
  std::map<int,std::set<int>>     m_unmet_deps;
  bool                            m_has_unmet_deps;
  bool                            m_IC_processed;

  // The nodes in the atm DAG
  std::vector<Node>               m_nodes;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_DAG_HPP
