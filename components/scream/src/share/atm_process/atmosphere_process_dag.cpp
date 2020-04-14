#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"

#include <fstream>

namespace scream {

void AtmProcDAG::create_dag(const group_type& atm_procs) {

  cleanup ();

  // begin is a placeholder for the beginning of an atm time step.
  // It comes handy when representing the input from prev. time step.
  Node begin_ts;
  begin_ts.name = "Begin of atm time step";
  begin_ts.id = 0;
  m_nodes.push_back(begin_ts);
  m_unmet_deps[0].resize(0);

  // Create the nodes
  add_nodes(atm_procs);

  // Add the end node. Just like begin, it is a placeholder, this time
  // for the 'next' time step.
  m_nodes.push_back(Node());
  Node& end_ts = m_nodes.back();
  end_ts.name = "End of atm time step";
  end_ts.id = m_nodes.size()-1;
  m_unmet_deps[end_ts.id].resize(0);

  // Next, check if some unmet deps are simply coming from previous time step.
  for (auto& it : m_unmet_deps) {
    int id = it.first;
    auto& unmet = it.second;
    for (size_t k=0; k<unmet.size();) {
      const std::string& fid = unmet[k];
      auto check = m_fid_to_last_provider.find(fid);
      if (check!=m_fid_to_last_provider.end()) {
        // We found this 'unmet' dependency. So a process, which appears
        // later in the time-step dag, will compute it. It's not a 'real'
        // unmet dependency, but rather a dependency on the prev. time step.
        m_nodes[0].computed.push_back(fid);
        m_nodes[0].children.push_back(id);
        m_nodes[id].parents.push_back(0);

        end_ts.required.push_back(check->first);
        end_ts.parents.push_back(check->second);
        m_nodes[check->second].children.push_back(end_ts.id);

        // Erase this field from the unmet deps of this node
        unmet.erase(unmet.begin()+k);
      } else {
        // Modify the name of the 'required' entry, adding ' *** MISSING *** '
        Node& n = m_nodes[id];
        auto req_it = std::find(n.required.begin(),n.required.end(),fid);
        scream_require_msg(req_it!=n.required.end(),
                           "[AtmProcDAG] Internal error! Please, contact developers.\n");
        *req_it += "  *** MISSING ***";

        // Note: we update the iterator only in this branch of the if,
        //       cause in the other case, erasing an entry already makes
        //       fid_it point to the 'next' element.
        ++k;
      }
    }
  }

  // Fill the set m_deps_on_prev_ts_set.
  for (const auto& it : m_deps_on_prev_ts) {
    const auto& deps = it.second;
    for (const auto& fid : deps) {
      m_deps_on_prev_ts_set.insert(fid);
    }
  }

  // Fill the set m_unmet_deps_set
  for (const auto& it : m_unmet_deps) {
    const auto& deps = it.second;
    for (const auto& fid : deps) {
      m_unmet_deps_set.insert(fid);
    }
  }
}

void AtmProcDAG::write_dag (const std::string& fname) const {
  std::ofstream ofile;
  ofile.open (fname.c_str());

  ofile << "digraph G {\n";

  for (const auto& n : m_nodes) {
    std::string color = "";
    // Check if all deps of this node are satisfied.
    if (m_unmet_deps.at(n.id).size()>0) {
      color = "red";
    }

    // Write node, with computed/required fields
    ofile << n.id
          << " [fontcolor=\"" << color
          << "\", label=\"" << n.name;
    if (n.name=="Begin of atm time step") {
      ofile << "\\n Inputs from previous time step:";
    } else if (n.name!="End of atm time step"){
      ofile << "\\n Computed:";
    }
    for (const auto& fid : n.computed) {
      ofile << "\\n   " << fid;
    }
    if (n.name=="End of atm time step") {
      ofile << "\\n Outputs for next time step:";
    } else if (n.name!="Begin of atm time step") {
      ofile << "\\n Required:";
    }
    for (const auto& fid : n.required) {
      ofile << "\\n   " << fid;
    }
    ofile << "\"]\n";

    // Write all outgoing edges
    for (const auto c : n.children) {
      ofile << n.id << "->" << c << "\n";
    }
  }

  // Close the file
  ofile << "}";
  ofile.close();
}

void AtmProcDAG::cleanup () {
  m_nodes.clear();
  m_fid_to_last_provider.clear();
  m_unmet_deps.clear();
  m_deps_on_prev_ts.clear();
  m_deps_on_prev_ts_set.clear();
}

void AtmProcDAG::add_nodes (const group_type& atm_procs) {
  
  const int num_procs = atm_procs.get_num_processes();
  const auto& remappers_in  = atm_procs.get_inputs_remappers ();
  const auto& remappers_out = atm_procs.get_outputs_remappers ();
  const bool sequential = (atm_procs.get_schedule_type()==ScheduleType::Sequential);

  scream_require_msg (sequential, "Error! Parallel splitting dag not yet supported.\n");

  int id = m_nodes.size();
  for (int i=0; i<num_procs; ++i) {
    const auto proc = atm_procs.get_process(i);
    const bool is_group = (proc->type()==AtmosphereProcessType::Group);
    if (is_group) {
      auto group = std::dynamic_pointer_cast<const group_type>(proc);
      scream_require_msg(group, "Error! Unexpected failure in dynamic_pointer_cast.\n"
                                "       Please, contact developers.\n");
      // Add all the stuff in the group.
      // Note: no need to add remappers for this process, because
      //       the sub-group will have its remappers taken care of
      add_nodes(*group);
    } else {
      // Add a node for the remapper(s) in (if needed)
      for (auto r : remappers_in[i]) {
        if (r.second->get_num_fields()>0) {
          const auto& rem = *r.second;
          m_nodes.push_back(Node());
          Node& node = m_nodes.back();;
          node.id = id;
          node.name = proc->name()+" (remap in [" + rem.get_src_grid()->name() + "->" + rem.get_tgt_grid()->name() + "])";
          m_unmet_deps[id].resize(0);
          for (int k=0; k<rem.get_num_fields(); ++k) {
            const auto& fid_in = rem.get_src_field_id(k).get_id_string();
            node.required.push_back(fid_in);
            auto it = m_fid_to_last_provider.find(fid_in);
            if (it==m_fid_to_last_provider.end()) {
              m_unmet_deps[id].push_back(fid_in);
            } else {
              // Establish parent-child relationship
              Node& parent = m_nodes[it->second];
              node.parents.push_back(parent.id);
              parent.children.push_back(node.id);
            }
            const auto& fid_out = rem.get_tgt_field_id(k).get_id_string();
            node.computed.push_back(fid_out);
            m_fid_to_last_provider[fid_out] = id;
          }
          ++id;
        }
      }

      // Note: the braces are just to have this chunk of code in a scope,
      //       and avoid vars clashing
      {
        // Create a node for the process
        // Node& node = m_nodes[proc->name()];
        m_nodes.push_back(Node());
        Node& node = m_nodes.back();;
        node.id = id;
        node.name = proc->name();
        m_unmet_deps[id].resize(0);
        // node.parents.push_back(m_current_head);
        // m_current_head->children.push_back(&node);
        for (auto f : proc->get_required_fields()) {
          const auto& fid = f.get_id_string();
          node.required.push_back(fid);
          auto it = m_fid_to_last_provider.find(fid);
          if (it==m_fid_to_last_provider.end()) {
            m_unmet_deps[id].push_back(fid);
          } else {
            // Establish parent-child relationship
            Node& parent = m_nodes[it->second];
            node.parents.push_back(parent.id);
            parent.children.push_back(node.id);
          }
        }
        for (auto f : proc->get_computed_fields()) {
          const auto& fid = f.get_id_string();
          node.computed.push_back(fid);
          m_fid_to_last_provider[fid] = id;
        }
        ++id;
      }

      // Add a node for the remapper(s) out (if needed)
      for (auto r : remappers_out[i]) {
        if (r.second->get_num_fields()>0) {
          const auto& rem = *r.second;
          m_nodes.push_back(Node());
          Node& node = m_nodes.back();
          node.id = id;
          node.name = proc->name()+" (remap out [" + rem.get_src_grid()->name() + "->" + rem.get_tgt_grid()->name() + "])";
          m_unmet_deps[id].resize(0);
          for (int k=0; k<rem.get_num_fields(); ++k) {
            const auto& fid_in = rem.get_src_field_id(k).get_id_string();
            node.required.push_back(fid_in);
            // A remapper for outputs of proc should *not* have unmet dependencies.
            auto it = m_fid_to_last_provider.find(fid_in);
            scream_require_msg (it!=m_fid_to_last_provider.end(),
                                "Internal error! Something is off with outputs remapper for atm proc '" + proc->name() + "'.\n"
                                "   Please, contact developers.\n");

            // Establish parent-child relationship
            Node& parent = m_nodes[it->second];
            node.parents.push_back(parent.id);
            parent.children.push_back(node.id);

            const auto& fid_out = rem.get_tgt_field_id(k).get_id_string();
            node.computed.push_back(fid_out);
            m_fid_to_last_provider[fid_out] = id;
          }
          ++id;
        }
      }
    }
  }
}

} // namespace scream
