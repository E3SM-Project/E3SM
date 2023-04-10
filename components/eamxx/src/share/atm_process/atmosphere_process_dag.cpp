#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"

#include <fstream>

namespace scream {

void AtmProcDAG::
create_dag(const group_type& atm_procs)
{
  cleanup ();

  // Create the nodes
  add_nodes(atm_procs);

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

void AtmProcDAG::add_surface_coupling (const std::set<FieldIdentifier>& imports,
                                       const std::set<FieldIdentifier>& exports)
{
  EKAT_REQUIRE_MSG (m_nodes.size()>0,
    "Error! You need to create the dag before adding surface coupling.\n");

  // Process all imports
  m_nodes.push_back(Node());
  auto& atm_imp = m_nodes.back();
  atm_imp.id = m_nodes.size()-1;
  atm_imp.name = "Atm Import";
  m_unmet_deps[atm_imp.id].clear();
  for (const auto& f : imports) {
    auto fid = add_fid(f);

    // Add the fid to the list of 'computed' fields
    atm_imp.computed.insert(fid);

    // Now, remove the unmet dependency (if any)
    for (auto& it : m_unmet_deps) {
      // Erase the unmet dependency (if any)
      int erased = it.second.erase(fid);

      if (erased==1) {
        // Establish parent-child relationship
        atm_imp.children.push_back(it.first);
      }
    }
  }

  // Process all exports
  m_nodes.push_back(Node());
  auto& atm_exp = m_nodes.back();
  atm_exp.id = m_nodes.size()-1;
  atm_exp.name = "Atm Export";
  m_unmet_deps[atm_exp.id].clear();
  for (const auto& f : exports) {
    auto fid = add_fid(f);

    // Add the fid to the list of 'computed' fields
    atm_exp.required.insert(fid);
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
        s += e2str(t);
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

      // Computed fields
      if (n.name=="Begin of atm time step") {
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Atm input fields from previous time step:</font></td></tr>\n";
      } else if (n.name!="End of atm time step"){
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Computed Fields:</font></td></tr>\n";
      }
      for (const auto& fid : n.computed) {
        std::string fc = "<font color=\"";
        fc += "black";
        fc += "\">  ";
        ofile << "      <tr><td align=\"left\">" << fc << html_fix(print_fid(m_fids[fid],fid_verb)) << "</font></td></tr>\n";
      }

      // Required fields
      if (n.name=="End of atm time step") {
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Atm output fields for next time step:</font></td></tr>\n";
      } else if (n.name!="Begin of atm time step") {
        ofile << "      <tr><td align=\"left\"><font color=\"blue\">Required Fields:</font></td></tr>\n";
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

      // Computed groups
      if (n.gr_computed.size()>0) {
        if (n.name=="Begin of atm time step") {
          ofile << "      <tr><td align=\"left\"><font color=\"blue\">Atm Input groups:</font></td></tr>\n";
        } else if (n.name!="End of atm time step"){
          ofile << "      <tr><td align=\"left\"><font color=\"blue\">Computed Groups:</font></td></tr>\n";
        }
        for (const auto& gr_fid : n.gr_computed) {
          std::string fc = "<font color=\"";
          fc += (ekat::contains(unmet,gr_fid) ? "red" : "black");
          fc += "\">  ";
          ofile << "      <tr><td align=\"left\">" << fc << html_fix(print_fid(m_fids[gr_fid],fid_verb));
          ofile << "</font></td></tr>\n";
          if (verbosity>2) {
            ofile << "      <tr><td align=\"left\">  Members:";
            const auto& group = m_gr_fid_to_group.at(m_fids[gr_fid]);
            const auto& members = group.m_fields;
            const auto& members_names = group.m_info->m_fields_names;
            const size_t max_len = 40;
            size_t len = 0;
            size_t i = 0;
            for (const auto& fn : members_names) {
              const auto f = members.at(fn);
              const auto& mfid = f->get_header().get_identifier();
              const auto mfid_id = get_fid_index(mfid);
              std::string mfc = "<font color=\"";
              mfc += (ekat::contains(unmet,mfid_id) ? "red" : "black");
              mfc += "\">";
              if (len>0) {
                ofile << ",";
                ++len;
              }
              ofile << " " << mfc << fn << "</font>";
              len += fn.size()+1;
              if (len>max_len) {
                if (i<members.size()-1) {
                  ofile << ",";
                }
                ofile << "</td></tr>";
                ofile << "      <tr><td align=\"left\">                  ";
                len = 0;
              }
              ++i;
            }
            ofile <<  "</td></tr>";
          }
        }
      }

      // Required groups
      if (n.gr_required.size()>0) {
        if (n.name=="End of atm time step") {
          ofile << "      <tr><td align=\"left\"><font color=\"blue\">Atm Output Groups:</font></td></tr>\n";
        } else if (n.name!="Begin of atm time step") {
          ofile << "      <tr><td align=\"left\"><font color=\"blue\">Required Groups:</font></td></tr>\n";
        }
        for (const auto& gr_fid : n.gr_required) {
          std::string fc = "<font color=\"";
          fc += (ekat::contains(unmet,gr_fid) ? "red" : "black");
          fc += "\">  ";
          ofile << "      <tr><td align=\"left\">" << fc << html_fix(print_fid(m_fids[gr_fid],fid_verb));
          if (ekat::contains(m_unmet_deps.at(n.id),gr_fid)) {
            ofile << "<b>  *** MISSING ***</b>";
          }
          ofile << "</font></td></tr>\n";
          if (verbosity>2) {
            ofile << "      <tr><td align=\"left\">  Members:";
            const auto& group = m_gr_fid_to_group.at(m_fids[gr_fid]);
            const auto& members = group.m_fields;
            const auto& members_names = group.m_info->m_fields_names;
            const size_t max_len = 40;
            size_t len = 0;
            size_t i = 0;
            for (const auto& fn : members_names) {
              const auto f = members.at(fn);
              const auto& mfid = f->get_header().get_identifier();
              const auto mfid_id = get_fid_index(mfid);
              std::string mfc = "<font color=\"";
              mfc += (ekat::contains(unmet,mfid_id) ? "red" : "black");
              mfc += "\">";
              if (len>0) {
                ofile << ",";
                ++len;
              }
              ofile << " " << mfc << fn << "</font>";
              len += fn.size()+1;
              if (len>max_len) {
                if (i<members.size()-1) {
                  ofile << ",";
                }
                ofile << "</td></tr>";
                ofile << "      <tr><td align=\"left\">                  ";
                len = 0;
              }
              ++i;
            }
            ofile <<  "</td></tr>";
          }
        }
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
add_nodes (const group_type& atm_procs)
{
  const int num_procs = atm_procs.get_num_processes();
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
      add_nodes(*group);
    } else {
      // Create a node for the process
      // Node& node = m_nodes[proc->name()];
      m_nodes.push_back(Node());
      Node& node = m_nodes.back();;
      node.id = id;
      node.name = proc->name();
      m_unmet_deps[id].clear();

      // Input fields
      for (const auto& f : proc->get_fields_in()) {
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

      // Output fields
      for (const auto& f : proc->get_fields_out()) {
        const auto& fid = f.get_header().get_identifier();
        const int fid_id = add_fid(fid);
        node.computed.insert(fid_id);
        m_fid_to_last_provider[fid_id] = id;
      }

      // Input groups
      for (const auto& group : proc->get_groups_in()) {
        if (!group.m_info->m_bundled) {
          // Group is not bundled: process fields individually
          for (const auto& it_f : group.m_fields) {
            const auto& fid = it_f.second->get_header().get_identifier();
            const int fid_id = add_fid(fid);
            node.computed.insert(fid_id);
            m_fid_to_last_provider[fid_id] = id;
          }
        } else {
          // Group is bundled: process the bundled field
          const auto& gr_fid = group.m_bundle->get_header().get_identifier();
          const int gr_fid_id = add_fid(gr_fid);
          node.gr_required.insert(gr_fid_id);
          auto it = m_fid_to_last_provider.find(gr_fid_id);
          if (it==m_fid_to_last_provider.end()) {
            // It might still be ok, as long as there is a provider for all the fields in the group
            bool all_members_have_providers = true;
            for (auto it_f : group.m_fields) {
              const auto& fid = it_f.second->get_header().get_identifier();
              const int fid_id = add_fid(fid);
              auto it_p = m_fid_to_last_provider.find(fid_id);
              if (it_p==m_fid_to_last_provider.end()) {
                m_unmet_deps[id].insert(fid_id);
                all_members_have_providers = false;
              } else {
                Node& parent = m_nodes[it_p->second];
                parent.children.push_back(node.id);
              }
            }
            if (!all_members_have_providers) {
              m_unmet_deps[id].insert(gr_fid_id);
            }
          } else {
            // Establish parent-child relationship
            Node& parent = m_nodes[it->second];
            parent.children.push_back(node.id);
          }
          m_gr_fid_to_group.emplace(gr_fid,group);
        }
      }

      // Output groups
      for (const auto& group : proc->get_groups_out()) {
        if (!group.m_info->m_bundled) {
          // Group is not bundled: process fields in the group individually
          for (const auto& it_f : group.m_fields) {
            const auto& fid = it_f.second->get_header().get_identifier();
            const int fid_id = add_fid(fid);
            node.computed.insert(fid_id);
            m_fid_to_last_provider[fid_id] = id;
          }
        } else {
          // Group is bundled: process the bundled field
          const auto& gr_fid = group.m_bundle->get_header().get_identifier();
          const int gr_fid_id = add_fid(gr_fid);
          node.gr_computed.insert(gr_fid_id);
          m_fid_to_last_provider[gr_fid_id] = id;
          m_gr_fid_to_group.emplace(gr_fid,group);

          // Additionally, each field in the group is implicitly 'computed'
          // by this node, so update their last provider
          for (auto it_f : group.m_fields) {
            const auto& fid = it_f.second->get_header().get_identifier();
            const int fid_id = add_fid(fid);
            m_fid_to_last_provider[fid_id] = id;
          }
        }
      }
      ++id;
    }
  }
}

int AtmProcDAG::add_fid (const FieldIdentifier& fid) {
  auto it = ekat::find(m_fids,fid);
  if (it==m_fids.end()) {
    m_fids.push_back(fid);
    return m_fids.size()-1;
  } else {
    return std::distance(m_fids.begin(), it);
  }
}

int AtmProcDAG::get_fid_index (const FieldIdentifier& fid) const {
  auto it = ekat::find(m_fids,fid);
  if (it==m_fids.end()) {
    return -1;
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
