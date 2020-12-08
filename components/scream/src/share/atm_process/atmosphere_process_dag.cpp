#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"

#include <fstream>

namespace scream {

void AtmProcDAG::
create_dag(const group_type& atm_procs,
           const std::shared_ptr<FieldRepository<Real>> field_repo) {

  cleanup ();

  // Create the nodes
  add_nodes(atm_procs,field_repo);

  // Add a 'begin' and 'end' placeholders. While they are not actual
  // nodes of the graph, they come handy when representing inputs
  // coming from previous time step, or output that will be fed to
  // the next time step
  m_nodes.reserve(m_nodes.size()+2);

  m_nodes.push_back(Node());
  Node& begin_ts = m_nodes.back();
  begin_ts.name = "Begin of atm time step";
  begin_ts.id = m_nodes.size()-1;
  m_unmet_deps[begin_ts.id].clear();

  m_nodes.push_back(Node());
  Node& end_ts = m_nodes.back();
  end_ts.name = "End of atm time step";
  end_ts.id = m_nodes.size()-1;
  m_unmet_deps[end_ts.id].clear();

  // Next, check if some unmet deps are simply coming from previous time step.
  for (auto& it : m_unmet_deps) {
    int id = it.first;
    auto& unmet = it.second;
    std::set<int> to_be_erased;
    for (auto fid : unmet) {
      auto check = m_fid_to_last_provider.find(fid);
      if (check!=m_fid_to_last_provider.end()) {
        // We found this 'unmet' dependency. So a process, which appears
        // later in the time-step dag, will compute it. It's not a 'real'
        // unmet dependency, but rather a dependency on the prev. time step.
        begin_ts.computed.insert(fid);
        begin_ts.children.push_back(id);

        end_ts.required.insert(check->first);
        m_nodes[check->second].children.push_back(end_ts.id);

        // This fid is not really an unmet dep. Mark it to be erased from unmet.
        to_be_erased.insert(fid);
      }
    }

    for (auto fid : to_be_erased) {
      ekat::erase(unmet,fid);
    }
  }

  update_unmet_deps ();
}

void AtmProcDAG::add_field_initializer (const FieldInitializer& initializer)
{
  EKAT_REQUIRE_MSG (m_nodes.size()>0,
    "Error! You need to create the dag before adding field initializers.\n");

  const auto& inited_fields = initializer.get_inited_fields();

  // Add a node
  m_nodes.push_back(Node());
  auto& n = m_nodes.back();
  n.id = m_nodes.size()-1;
  n.name = initializer.name() + " (init only)";
  m_unmet_deps[n.id].clear();

  for (const auto& f : inited_fields) {
    auto fid = add_fid(f);

    // Add the fid to the list of 'computed' fields
    n.computed.insert(fid);

    // Now, remove the unmet dependency (if any)
    for (auto& it : m_unmet_deps) {
      // Erase the unmet dependency (if any)
      int erased = it.second.erase(fid);

      if (erased==1) {
        // Establish parent-child relationship
        n.children.push_back(it.first);
      }
    }
  }

  // We need to re-check whether there are unmet deps
  update_unmet_deps();
}

void AtmProcDAG::write_dag (const std::string& fname, const int verbosity) const {

  if (verbosity<=0) {
    return;
  }

  // Handy lambda to print a fid with different degrees of verbosity
  auto print_fid = [] (const FieldIdentifier& fid, const int verbosity) -> std::string {
    // Print will always print the name
    std::string s = fid.name();

    // If verbosity is >= 1, add layout info
    if (verbosity>0) {
      // If verbosity is >= 2, add grid name
      if (verbosity>1) {
        s += " [" + fid.get_grid_name() + "]";
      }
      s += " <";
      for (auto t : fid.get_layout().tags()) {
        s += tag2string(t);
        s += ",";
      }
      // Remove last ',' and add '>'.
      s.back() = '>';
      s += "(";
      for (auto dim : fid.get_layout().dims()) {
        s += std::to_string(dim);
        s += ",";
      }
      // Remove last ',' and add ')'.
      s.back() = ')';

      if (verbosity>2) {
        s += " [" + fid.get_units().get_string() + "]";
      }
    }
    return s;
  };

  auto html_fix = [] (const std::string& s) -> std::string {
    std::string out(s);
    auto pos = out.find('<');
    while (pos!=std::string::npos) {
      out.replace(pos,1,"&lt;");
      pos = out.find('<');
    }
    pos = out.find('>');
    while (pos!=std::string::npos) {
      out.replace(pos,1,"&gt;");
      pos = out.find('>');
    }
    return out;
  };

  std::ofstream ofile;
  ofile.open (fname.c_str());

  ofile << "strict digraph G {\n";

  for (const auto& n : m_nodes) {
    const auto& unmet = m_unmet_deps.at(n.id);

    // Write node, with computed/required fields
    ofile << n.id
          << " [\n"
          << "  shape=box\n"
          << "  label=<\n"
          << "    <table border=\"0\">\n"
          << "      <tr><td><b>" << html_fix(n.name) << "</b></td></tr>";
    if (verbosity>1) {
      // FieldIntentifier prints bare min with verb 0.
      // DAG starts printing fids with verb 2, so fid verb is verb-2;
      int fid_verb = verbosity-2;
      ofile << "<hr/>\n";
      if (n.name=="Begin of atm time step") {
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Inputs from previous time step:</font></td></tr>\n";
      } else if (n.name!="End of atm time step"){
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Computed:</font></td></tr>\n";
      }
      for (const auto& fid : n.computed) {
        ofile << "      <tr><td align=\"left\">  " << html_fix(print_fid(m_fids[fid],fid_verb)) << "</td></tr>\n";
      }
      if (n.name=="End of atm time step") {
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Outputs for next time step:</font></td></tr>\n";
      } else if (n.name!="Begin of atm time step") {
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Required:</font></td></tr>\n";
      }
      for (const auto& fid : n.required) {
        std::string fc = "<font color=\"";
        fc += (ekat::contains(unmet,fid) ? "red" : "black");
        fc += "\">  ";
        ofile << "      <tr><td align=\"left\">" << fc << html_fix(print_fid(m_fids[fid],fid_verb));
        if (ekat::contains(m_unmet_deps.at(n.id),fid)) {
          ofile << "<b>  *** MISSING ***</b>";
        }
        ofile << "</font></td></tr>\n";
      }
    } else {
      ofile << "\n";
    }
    ofile << "    </table>\n"
          << "  >\n];\n";

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
  m_has_unmet_deps = false;
}

void AtmProcDAG::
add_nodes (const group_type& atm_procs,
           const std::shared_ptr<FieldRepository<Real>> field_repo) {
  
  const int num_procs = atm_procs.get_num_processes();
  const auto& remappers_in  = atm_procs.get_inputs_remappers ();
  const auto& remappers_out = atm_procs.get_outputs_remappers ();
  const bool sequential = (atm_procs.get_schedule_type()==ScheduleType::Sequential);

  EKAT_REQUIRE_MSG (sequential, "Error! Parallel splitting dag not yet supported.\n");

  int id = m_nodes.size();
  for (int i=0; i<num_procs; ++i) {
    const auto proc = atm_procs.get_process(i);
    const bool is_group = (proc->type()==AtmosphereProcessType::Group);
    if (is_group) {
      auto group = std::dynamic_pointer_cast<const group_type>(proc);
      EKAT_REQUIRE_MSG(group, "Error! Unexpected failure in dynamic_pointer_cast.\n"
                                "       Please, contact developers.\n");
      // Add all the stuff in the group.
      // Note: no need to add remappers for this process, because
      //       the sub-group will have its remappers taken care of
      add_nodes(*group,field_repo);
    } else {
      // Add a node for the remapper(s) in (if needed)
      for (auto r : remappers_in[i]) {
        const auto& rem = *r.second;
        const int nfields = rem.get_num_registered_fields();
        if (nfields>0) {
          m_nodes.push_back(Node());
          Node& node = m_nodes.back();;
          node.id = id;
          node.name = proc->name()+" (remap in [" + rem.get_src_grid()->name() + "->" + rem.get_tgt_grid()->name() + "])";
          m_unmet_deps[id].clear();
          for (int k=0; k<nfields; ++k) {
            const auto& fid_in = rem.get_src_field_id(k);
            const int fid_in_id = add_fid(fid_in);
            node.required.insert(fid_in_id);
            auto it = m_fid_to_last_provider.find(fid_in_id);
            if (it==m_fid_to_last_provider.end()) {
              m_unmet_deps[id].insert(fid_in_id);
            } else {
              // Establish parent-child relationship
              Node& parent = m_nodes[it->second];
              parent.children.push_back(node.id);
            }

            const auto& fid_out = rem.get_tgt_field_id(k);
            const int fid_out_id = add_fid(fid_out);
            node.computed.insert(fid_out_id);
            m_fid_to_last_provider[fid_out_id] = id;
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
        m_unmet_deps[id].clear();
        for (auto fid : proc->get_required_fields()) {
          const int fid_id = add_fid(fid);
          node.required.insert(fid_id);
          auto it = m_fid_to_last_provider.find(fid_id);
          if (it==m_fid_to_last_provider.end()) {
            m_unmet_deps[id].insert(fid_id);
          } else {
            // Establish parent-child relationship
            Node& parent = m_nodes[it->second];
            parent.children.push_back(node.id);
          }
        }
        for (auto fid : proc->get_computed_fields()) {
          const int fid_id = add_fid(fid);
          node.computed.insert(fid_id);
          m_fid_to_last_provider[fid_id] = id;
        }
        for (auto itg : proc->get_updated_groups()) {
          EKAT_REQUIRE_MSG (field_repo, "Error! Field repo pointer is null.\n");
          auto group = field_repo->get_const_field_group(itg.first,itg.second);
          for (const auto& f : group) {
            const auto& fid = f.get_header().get_identifier();
            const int fid_id = add_fid(fid);
            node.required.insert(fid_id);
            auto it = m_fid_to_last_provider.find(fid_id);
            if (it==m_fid_to_last_provider.end()) {
              m_unmet_deps[id].insert(fid_id);
            } else {
              // Establish parent-child relationship
              Node& parent = m_nodes[it->second];
              parent.children.push_back(node.id);
            }
          }
        }
        for (auto itg : proc->get_required_groups()) {
          EKAT_REQUIRE_MSG (field_repo, "Error! Field repo pointer is null.\n");
          auto group = field_repo->get_const_field_group(itg.first,itg.second);
          for (const auto& f : group) {
            const auto& fid = f.get_header().get_identifier();
            const int fid_id = add_fid(fid);
            node.computed.insert(fid_id);
            m_fid_to_last_provider[fid_id] = id;
          }
        }
        ++id;
      }

      // Add a node for the remapper(s) out (if needed)
      for (auto r : remappers_out[i]) {
        const auto& rem = *r.second;
        const int nfields = rem.get_num_registered_fields();
        if (nfields>0) {
          m_nodes.push_back(Node());
          Node& node = m_nodes.back();
          node.id = id;
          node.name = proc->name()+" (remap out [" + rem.get_src_grid()->name() + "->" + rem.get_tgt_grid()->name() + "])";
          m_unmet_deps[id].clear();
          for (int k=0; k<nfields; ++k) {
            const auto& fid_in = rem.get_src_field_id(k);
            const int fid_in_id = add_fid(fid_in);
            node.required.insert(fid_in_id);
            // A remapper for outputs of proc should *not* have unmet dependencies.
            auto it = m_fid_to_last_provider.find(fid_in_id);
            EKAT_REQUIRE_MSG (it!=m_fid_to_last_provider.end(),
                                "Internal error! Something is off with outputs remapper for atm proc '" + proc->name() + "'.\n"
                                "   Please, contact developers.\n");

            // Establish parent-child relationship
            Node& parent = m_nodes[it->second];
            parent.children.push_back(node.id);

            const auto& fid_out = rem.get_tgt_field_id(k);
            const int fid_out_id = add_fid(fid_out);
            node.computed.insert(fid_out_id);
            m_fid_to_last_provider[fid_out_id] = id;
          }
          ++id;
        }
      }
    }
  }
}

int AtmProcDAG::add_fid (const FieldIdentifier& fid) {
  auto it = ekat::find(m_fids,fid);
  if (it==m_fids.end()) {
    m_fids.push_back(fid);
    return m_fids.size()-1;
  } else {
    return std::distance(m_fids.cbegin(),it);
  }
}

void AtmProcDAG::update_unmet_deps () {
  m_has_unmet_deps = false;

  for (const auto& it : m_unmet_deps) {
    const auto& unmet = it.second;
    m_has_unmet_deps |= (unmet.size()>0);
    if (m_has_unmet_deps) {
      break;
    }
  }
}

} // namespace scream
