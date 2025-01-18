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

  // Now all nodes are created, we need to create edges, by checking
  // which is the first node (if any) that computes any of the input
  // fields of each node
  add_edges ();

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
        s += " [" + fid.get_units().to_string() + "]";
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

  ofile << "strict digraph G {\n"
        << "rankdir=\"LR\"\n";

  bool has_IC_field = false;
  for (const auto& n : m_nodes) {
    const auto& unmet = m_unmet_deps.at(n.id);

    std::string box_fmt;
    int id_IC = -1;
    int id_begin = -1;
    int id_end = -1;
    if (n.name == "Begin of atm time step") {
      id_begin = n.id;
      box_fmt = "  color=\"#00667E\"\n  fontcolor=\"#00667E\"\n  style=filled\n"
                "  fillcolor=\"#b9d4dc\"\n";
    } else if (n.name == "Initial Conditions") {
      id_IC = n.id;
      box_fmt = "  color=\"#006219\"\n  fontcolor=\"#006219\"\n  style=filled\n"
                "  fillcolor=\"#b9dcc2\"\n";
    } else if (n.name == "End of atm time step") {
      id_end = n.id;
      box_fmt = "  color=\"#88621e\"\n  fontcolor=\"#88621e\"\n  style=filled\n"
                "  fillcolor=\"#dccfb9\"\n";
    }

    // Write node, with computed/required fields
    ofile << n.id
          << " [\n"
          << "  shape=box\n"
          << box_fmt
          << "  penwidth=4\n"
          << "  fontsize=30\n"
          << "  label=<\n"
          << "    <table border=\"0\">\n"
          << "      <tr><td><b><font point-size=\"40\">" << html_fix(n.name)
          << "</font></b></td></tr>\n";

    int sz_comp = n.computed.size(), sz_req = n.required.size(),
        sz_grcomp = n.gr_computed.size(), sz_grreq = n.gr_required.size();
    int nfield = sz_comp + sz_req + sz_grcomp + sz_grreq;
    if (verbosity > 1 && nfield > 0) {
      // FieldIntentifier prints bare min with verb 0.
      // DAG starts printing fids with verb 2, so fid verb is verb-2;
      int fid_verb = verbosity-2;
      ofile << "      <hr/>\n";

      if (sz_comp > 0) {
        // Computed fields
        if (n.id == id_begin) {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#00667E\">"
                << "Atm input fields from previous time step:</font></b></td></tr>\n";
        } else if (n.id == id_IC) {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#00667E\">"
                << "Initial Fields:</font></b></td></tr>\n";
        } else if (n.id != id_end) {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#88621e\">"
                << "Computed Fields:</font></b></td></tr>\n";
        }

        for (const auto& fid : n.computed) {
          std::string fc = "<font color=\"";
          int fid_out = std::abs(fid);
          fc += "black";
          fc += "\">  ";
          ofile << "      <tr><td align=\"left\">" << fc
                << html_fix(print_fid(m_fids[fid_out], fid_verb))
                << "</font></td></tr>\n";
        }
      }

      if (sz_req > 0) {
        // Required fields
        if (n.id == id_end) {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#88621e\">"
                << "Atm output fields for next time step:</font></b></td></tr>\n";
        } else if (n.id != id_begin && n.id != id_IC) {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#00667E\">"
                << "Required Fields:</font></b></td></tr>\n";
        }
        for (const auto& fid : n.required) {
          std::string fc = "<font color=\"";
          if (ekat::contains(unmet, fid)) {
            fc += "red";
          } else if (ekat::contains(unmet, -fid)) {
            fc +=  "#006219";
          } else {
            fc += "black";
          }
          fc += "\">  ";
          ofile << "      <tr><td align=\"left\">" << fc << html_fix(print_fid(m_fids[fid],fid_verb));
          if (ekat::contains(unmet, fid)) {
            ofile << "<b>  *** MISSING ***</b>";
          } else if (ekat::contains(unmet, -fid)) {
            ofile << "<b>  (Init. Cond.)</b>";
            has_IC_field = true;
          }
          ofile << "</font></td></tr>\n";
        }
      }

      // Computed groups
      if (sz_grcomp > 0) {
        if (n.id == id_begin) {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#00667E\">Atm Input groups:</font></b></td></tr>\n";
        } else if (n.id != id_end){
          ofile << "      <tr><td align=\"left\"><b><font color=\"#88621e\">Computed Groups:</font></b></td></tr>\n";
        }
        for (const auto& gr_fid : n.gr_computed) {
          std::string fc = "<font color=\"";
          fc += "black";
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
              std::string mfc = "<font color=\"";
              mfc += "black";
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
      if (sz_grreq > 0) {
        if (n.name=="End of atm time step") {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#00667E\">Atm Output Groups:</font></b></td></tr>\n";
        } else if (n.name!="Begin of atm time step") {
          ofile << "      <tr><td align=\"left\"><b><font color=\"#00667E\">Required Groups:</font></b></td></tr>\n";
        }
        for (const auto& gr_fid : n.gr_required) {
          std::string fc = "<font color=\"";
          fc += (ekat::contains(unmet,gr_fid) ? "red" : "black");
          fc += "\">  ";
          ofile << "      <tr><td align=\"left\">" << fc << html_fix(print_fid(m_fids[gr_fid],fid_verb));
          if (ekat::contains(unmet, gr_fid)) {
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
      ofile << n.id << "->" << c << "[penwidth=4];\n";
    }
  }

  if (!m_IC_processed && m_has_unmet_deps) {
    int this_node_id = m_nodes.size() + 1;
    ofile << this_node_id << " [\n"
          << "  shape=box\n"
          << "  color=\"#605d57\"\n"
          << "  fontcolor=\"#034a4a\"\n"
          << "  penwidth=8\n"
          << "  fontsize=40\n"
          << "  style=filled\n"
          << "  fillcolor=\"#999999\"\n"
          << "  align=\"center\"\n"
          << "  label=<<b><font color=\"#774006\">NOTE:</font> "
             "Fields marked missing may be<br align=\"center\"/>provided by "
             "the as-yet-unprocessed<br align=\"center\"/>initial condition</b>>\n"
          << "];\n";
  }

  if (m_IC_processed && has_IC_field) {
    int this_node_id = m_nodes.size() + 1;
    ofile << this_node_id << " [\n"
          << "  shape=box\n"
          << "  color=\"#605d57\"\n"
          << "  fontcolor=\"#031576\"\n"
          << "  penwidth=8\n"
          << "  fontsize=40\n"
          << "  style=filled\n"
          << "  fillcolor=\"#cccccc\"\n"
          << "  align=\"center\"\n"
          << "  label=<<b><font color=\"#3d2906\">NOTE:</font> Fields denoted "
             "with <font color=\"#006219\"><b>green text</b></font> "
             "<br align=\"center\"/>indicate the field was provided by the "
             "<br align=\"center\"/>initial conditions and never updated</b>>\n"
          << "];\n";
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
  m_IC_processed = false;
}

void AtmProcDAG::
add_nodes (const group_type& atm_procs)
{
  const int num_procs = atm_procs.get_num_processes();
  const bool sequential = (atm_procs.get_schedule_type()==ScheduleType::Sequential);

  EKAT_REQUIRE_MSG (sequential, "Error! Parallel splitting dag not yet supported.\n");

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
      int id = m_nodes.size();
      m_nodes.push_back(Node());
      Node& node = m_nodes.back();
      node.id = id;
      node.name = proc->name();
      m_unmet_deps[id].clear(); // Ensures an entry for this id is in the map

      // Input fields
      for (const auto& f : proc->get_fields_in()) {
        const auto& fid = f.get_header().get_identifier();
        const int fid_id = add_fid(fid);
        node.required.insert(fid_id);
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
    }
  }
}

void AtmProcDAG::add_edges () {
  for (auto& node : m_nodes) {
    // First individual input fields. Add this node as a children
    // of any *previous* node that computes them. If none provides
    // them, add to the unmet deps list
    for (auto id : node.required) {
      auto it = m_fid_to_last_provider.find(id);
      // Note: check that last provider id is SMALLER than this node id
      if (it!=m_fid_to_last_provider.end() and it->second<node.id) {
        auto parent_id = it->second;
        m_nodes[parent_id].children.push_back(node.id);
      } else {
        m_unmet_deps[node.id].insert(id);
      }
    }
    // Then process groups, looking at both the bundled field and individual fields.
    // NOTE: we don't know if the group as a whole is the last to be updated
    //       OR if each group member is updated after the last "group-update".
    //       So get the id of the last node that updates each field and the group,
    //       and use the most recent one
    for (auto id : node.gr_required) {
      const auto& gr_fid = m_fids[id];
      const auto& group = m_gr_fid_to_group.at(gr_fid);
      const int   size  = group.m_info->size();

      int last_group_update_id = -1;
      std::vector<int> last_members_update_id(size,-1);

      // First check when the group as a whole was last updated
      auto it = m_fid_to_last_provider.find(id);
      // Note: check that last provider id is SMALLER than this node id
      if (it!=m_fid_to_last_provider.end() and it->second<node.id) {
        last_group_update_id = it->second;
      }
      // Then check when each group member was last updated
      int i=0;
      for (auto f_it : group.m_fields) {
        const auto& fid = f_it.second->get_header().get_identifier();
        auto fid_id = std::find(m_fids.begin(),m_fids.end(),fid) - m_fids.begin();
        it = m_fid_to_last_provider.find(fid_id);
        // Note: check that last provider id is SMALLER than this node id
        if (it!=m_fid_to_last_provider.end() and it->second<node.id) {
          last_members_update_id[i] = it->second;
        }
        ++i;
      }

      auto min = *std::min_element(last_members_update_id.begin(),last_members_update_id.end());
      if (min==-1 && last_group_update_id==-1) {
        // Nobody computed the group as a whole and some member was not computed
        m_unmet_deps[node.id].insert(id);
      } else if (min>last_group_update_id) {
        // All members are updated after the group
        for (auto fid_id : last_members_update_id) {
          m_nodes[fid_id].children.push_back(node.id);
        }
      } else {
        // Add the group provider as a parent, but also the provider of each
        // field which is updated after the group
        m_nodes[id].children.push_back(node.id);
        for (auto fid_id : last_members_update_id) {
          if (fid_id>last_group_update_id) {
            m_nodes[fid_id].children.push_back(node.id);
          }
        }
      }
    }
  }
}

void AtmProcDAG::process_initial_conditions(const grid_field_map &ic_inited) {
  // First, add the fields that were determined to come from the previous time
  // step => IC for t = 0

  // Create a node for the ICs by copying the begin_node. Recall that so far
  // m_nodes contains [<processes>, 'beg-of-step', 'end-of-step']
  // WARNING: do NOT get a ref to beg-of-step node, since calls to m_nodes.push_back
  //          may resize the vector, invalidating the reference.
  auto begin_node = m_nodes[m_nodes.size()-2];
  auto& ic_node = m_nodes.emplace_back(begin_node);

  // now set/clear the basic data for the ic_node
  int id = m_nodes.size();
  ic_node.id = id;
  ic_node.name = "Initial Conditions";
  m_unmet_deps[id].clear();
  ic_node.children.clear();
  // now add the begin_node as a child of the ic_node
  ic_node.children.push_back(begin_node.id);
  // return if there's nothing to process in the ic_inited vector
  if (ic_inited.size() == 0) {
    return;
  }
  std::set<int> to_be_marked;
  for (auto &node : m_nodes) {
    if (m_unmet_deps.at(node.id).empty()) {
      continue;
    } else {
      // NOTE: node_unmet_fields is a std::set<int>
      auto &node_unmet_fields = m_unmet_deps.at(node.id);
      // add the current node as a child of the IC node
      ic_node.children.push_back(node.id);
      for (auto &um_fid : node_unmet_fields) {
        for (auto &it1 : ic_inited) {
          const auto &grid_name = it1.first;
          // if this unmet-dependency field's name is in the ic_inited map for
          // the provided grid_name key, we record the field id in to_be_marked
          // (because changing it messes up the iterator)
          if (ekat::contains(ic_inited.at(grid_name), m_fids[um_fid].name())) {
            to_be_marked.insert(um_fid);
            // add the fid of the formerly unmet dep to the initial condition
            // node's computed list
            ic_node.computed.insert(um_fid);
          } else {
            continue;
          }
        }
      }
      if (to_be_marked.empty()) {
        continue;
      } else {
        // change the previously unmet dependency's field id to be negative,
        // indicating that it is now met and provided by the initial condition
        for (auto &fid : to_be_marked) {
          node_unmet_fields.erase(fid);
          node_unmet_fields.insert(-fid);
        }
      }
    }
  }
  m_IC_processed = true;
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
