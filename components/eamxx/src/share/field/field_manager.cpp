#include "share/field/field_manager.hpp"

namespace scream
{

FieldManager::
FieldManager (const std::shared_ptr<const AbstractGrid>& grid)
 : FieldManager (std::make_shared<LibraryGridsManager>(grid))
{
  // Nothing else to do
}

FieldManager::
FieldManager (const std::shared_ptr<const GridsManager>& gm)
 : m_grids_mgr (gm)
{
  EKAT_REQUIRE_MSG (m_grids_mgr!=nullptr,
      "Error! Input grids manager pointer is not valid.");

  // For each grid, initialize maps
  for (auto gname : m_grids_mgr->get_grid_names()) {
    m_fields[gname] = std::map<ci_string,std::shared_ptr<Field>>();
    m_field_groups[gname] = std::map<ci_string, std::shared_ptr<FieldGroup>>();
    m_group_requests[gname] = std::map<std::string, std::set<GroupRequest>>();
  }

  m_repo_state = RepoState::Clean;
}

void FieldManager::register_field (const FieldRequest& req)
{
  const auto& id = req.fid;
  const auto& grid_name = id.get_grid_name();

  // Sanity checks
  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Clean,
      "Error! Repo state is not 'Open' yet. You must call registration_begins() first.\n");
  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Closed,
      "Error! Repo state is not 'Open' anymore. You already called registration_ends().\n");

  // Make sure this FM contains a grid corresponding to the input grid name
  EKAT_REQUIRE_MSG(m_grids_mgr->has_grid(grid_name),
    "Error! Attempting to register field on grid not in the FM's grids manager:\n"
    "  - FieldRequest grid: " + grid_name + "\n"
    "  - Grids stored by FM: " + m_grids_mgr->print_available_grids() + "\n");

  // FieldManager does not allow incomplete requests
  EKAT_REQUIRE_MSG(not req.incomplete,
    "Error! FieldManager does not allow registering incomplete FieldRequest.\n");

  // Get or create the new field
  if (!has_field(id.name(), grid_name)) {

    EKAT_REQUIRE_MSG (id.data_type()==field_valid_data_types().at<Real>(),
        "Error! While refactoring, we only allow the Field data type to be Real.\n"
        "       If you're done with refactoring, go back and fix things.\n");
    m_fields[grid_name][id.name()] = std::make_shared<Field>(id);
  } else {
    // Make sure the input field has the same layout and units as the field already stored.
    // TODO: this is the easiest way to ensure everyone uses the same units.
    //       However, in the future, we *may* allow different units, providing
    //       the users with conversion routines perhaps.
    const auto id0 = m_fields[grid_name][id.name()]->get_header().get_identifier();
    EKAT_REQUIRE_MSG(id.get_units()==id0.get_units(),
        "Error! Field '" + id.name() + "' already registered with different units:\n"
        "         - input field units:  " + id.get_units().to_string() + "\n"
        "         - stored field units: " + id0.get_units().to_string() + "\n"
        "       Please, check and make sure all atmosphere processes use the same units.\n");

    EKAT_REQUIRE_MSG(id.get_layout()==id0.get_layout(),
        "Error! Field '" + id.name() + "' already registered with different layout:\n"
        "         - input id:  " + id.get_id_string() + "\n"
        "         - stored id: " + id0.get_id_string() + "\n"
        "       Please, check and make sure all atmosphere processes use the same layout for a given field.\n");
  }

  // Make sure the field can accommodate the requested value type
  m_fields[grid_name][id.name()]->get_header().get_alloc_properties().request_allocation(req.pack_size);

  // Finally, add the field to the given groups
  // Note: we do *not* set the group info struct in the field header yet.
  //       we will do that when we end the registration phase.
  //       The reason is that at registration_ends we will know *all* the groups
  //       that each field belongs to.
  for (const auto& group_name : req.groups) {
    // Get group (and init ptr, if necessary)
    auto& group_info = m_field_group_info[group_name];
    if (not group_info) {
      group_info = std::make_shared<FieldGroupInfo>(group_name);
    }

    // Add the field name to the list of fields belonging to this group
    if (not ekat::contains(group_info->m_fields_names,id.name())) {
      group_info->m_fields_names.push_back(id.name());
    }

    // Add which grid has registered this field
    if (not ekat::contains(group_info->m_grid_registered[id.name()], grid_name)) {
      group_info->m_grid_registered[id.name()].push_back(grid_name);
    }

    // Ensure that each group in m_field_group_info also appears in the m_group_requests map by
    // adding a "trivial" GroupRequest for this group, meaning no monolithic allocation and pack size 1.
    register_group(GroupRequest(group_name, grid_name));
  }
}

void FieldManager::register_group (const GroupRequest& req)
{
  // Make sure this FM contains a grid corresponding to the input grid name
  EKAT_REQUIRE_MSG(m_grids_mgr->has_grid(req.grid),
    "Error! Attempting to register group on grid not in the FM's grids manager:\n"
    "  - GroupRequest grid:  " + req.grid + "\n"
    "  - Grids stored by FM: " + m_grids_mgr->print_available_grids() + "\n");

  // Groups have to be handled once registration is over, so for now simply store the request,
  // and create an empty group info (if it does not already exist)
  m_group_requests[req.grid][req.name].insert(req);

  auto& group_info = m_field_group_info[req.name];
  if (not group_info) {
    group_info = std::make_shared<FieldGroupInfo>(req.name);
    group_info->m_monolithic_allocation = (req.monolithic_alloc==MonolithicAlloc::Required);
  } else if (not group_info->m_monolithic_allocation) {
    group_info->m_monolithic_allocation = (req.monolithic_alloc==MonolithicAlloc::Required);
  }

  if (not ekat::contains(group_info->m_requested_grids, req.grid)) {
    group_info->m_requested_grids.push_back(req.grid);
  }
}

void FieldManager::
add_to_group (const std::string& field_name, const std::string& grid_name, const std::string& group_name)
{
  EKAT_REQUIRE_MSG (m_repo_state==RepoState::Closed,
      "Error! You cannot call 'add_to_group' until after 'registration_ends' has been called.\n");
  auto& group = m_field_group_info[group_name];
  if (not group) {
    group = std::make_shared<FieldGroupInfo>(group_name);
  }
  EKAT_REQUIRE_MSG (not group->m_monolithic_allocation,
      "Error! Cannot add fields to a group that has a monolithic allocation.\n"
      "   field name: " + field_name + "\n"
      "   group name: " + group_name + "\n");

  EKAT_REQUIRE_MSG (has_field(field_name, grid_name),
      "Error! Cannot add field to group, since the field is not present in this FieldManager.\n"
      "   field name: " + field_name + "\n"
      "   grid name:  " + grid_name + "\n"
      "   group name: " + group_name + "\n");

  if (not ekat::contains(group->m_fields_names, field_name)) {
    group->m_fields_names.push_back(field_name);
  }
  if (not ekat::contains(group->m_grid_registered[field_name], grid_name)) {
    group->m_grid_registered[field_name].push_back(grid_name);
  }

  auto& ft = get_field(field_name, grid_name).get_header().get_tracking();
  ft.add_to_group(group);
}

bool FieldManager::
has_field (const std::string& field_name, const std::string& grid_name) const {
  return m_fields.at(grid_name).find(field_name)!=m_fields.at(grid_name).end();
}

bool FieldManager::
has_group (const std::string& group_name, const std::string& grid_name) const {
  // Return true if there exists a group "name" and group
  // "name" has at least one field registered on "grid_name"
  if (m_field_group_info.count(group_name)>0) {
    auto& group_info = m_field_group_info.at(group_name);
    for (auto& fname : group_info->m_fields_names) {
      if (ekat::contains(group_info->m_grid_registered.at(fname), grid_name)) {
        return true;
      }
    }
  }
  return false;
}

Field FieldManager::get_field (const std::string& name, const std::string& grid_name) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
    "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  EKAT_REQUIRE_MSG(has_field(name, grid_name),
    "Error! Field " + name + " on grid " + grid_name + " not found.\n");

  return *m_fields.at(grid_name).at(name);
}

Field& FieldManager::get_field (const std::string& name, const std::string& grid_name) {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
    "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  EKAT_REQUIRE_MSG(has_field(name, grid_name),
    "Error! Field " + name + " on grid " + grid_name + " not found.\n");

  return *m_fields.at(grid_name).at(name);
}

FieldGroupInfo FieldManager::
get_group_info (const std::string& group_name, const std::string& grid_name) const
{
  // Sanity checks
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
    "Error! Cannot get field groups from the repo while registration has not yet completed.\n");
  EKAT_REQUIRE_MSG (has_group(group_name, grid_name),
    "Error! Field group '" + group_name + "' on grid '" + grid_name + "' not found.\n");

  auto info = *m_field_group_info.at(group_name);

  if (info.m_monolithic_allocation) {
    // All fields in a group with a monolithic allocation should exist on any grid in that group
    for (const auto& fname : info.m_fields_names) {
      EKAT_REQUIRE_MSG(has_field(fname, grid_name),
        "Internal FieldManager Error! Groups with monolithic allocation must contain all "
        "fields in m_fields_names on any grid from m_grids_in_group.\n"
        "The following field should, but does not, exist:\n"
        "  - Group: " + group_name + "\n"
        "  - Field: " + fname + "\n"
        "  - Grid:  " + grid_name + "\n");
    }
  } else {
    // For all other groups, remove fields in the group
    // that do not exist on the requested grid
    std::list<ci_string> remove_fnames;
    for (const auto& fname : info.m_fields_names) {
      if (not ekat::contains(info.m_grid_registered.at(fname), grid_name)) {
        remove_fnames.push_back(fname);
      }
    }
    for (auto fname : remove_fnames) {
      info.m_fields_names.remove(fname);
      info.m_grid_registered.erase(fname);
      info.m_subview_idx.erase(fname);
    }
  }

  return info;
}

FieldGroup FieldManager::
get_field_group (const std::string& group_name, const std::string& grid_name)
{
  // Sanity checks
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");
  EKAT_REQUIRE_MSG (has_group(group_name, grid_name),
      "Error! Field group '" + group_name + "' on grid '" + grid_name + "' not found.\n");

  // If FieldGroup has already been built, return that
  if (m_field_groups.at(grid_name).find(group_name) != m_field_groups.at(grid_name).end()) {
    return *m_field_groups.at(grid_name).at(group_name);
  }

  // FieldGroup does not yet exist, create entry in m_field_groups
  auto& group = m_field_groups[grid_name][group_name];
  group = std::make_shared<FieldGroup>(group_name);

  // Set the info in the group
  group->m_info = std::make_shared<FieldGroupInfo>(get_group_info(group_name, grid_name));

  // Find all the fields registered on given grid and set them in the group
  for (const auto& fname : group->m_info->m_fields_names) {
    group->m_individual_fields[fname] = m_fields.at(grid_name).at(fname);
  }

  // Fetch the monolithic field (if applicable)
  if (group->m_info->m_monolithic_allocation) {
    // All fields in a group have the same parent, get the parent from the 1st one
    const auto& parent_header = group->m_individual_fields.begin()->second->get_header().get_parent();

    EKAT_REQUIRE_MSG(parent_header!=nullptr,
      "Error! A field belonging to a field group with monolithic allocation is missing its 'parent'.\n");

    const auto& parent_id = parent_header->get_identifier();

    EKAT_REQUIRE_MSG(has_field(parent_id.name(), grid_name),
      "Internal FieldManager Error! Parent field must exist in the FieldManager:\n"
      "  - Group: " + group_name + "\n"
      "  - Bundeled field: " + parent_id.name() + "\n"
      "  - Grid:  " + grid_name + "\n");

    const auto& parent_field = m_fields.at(grid_name).at(parent_id.name());

    // Get the parent field of one of the group fields and check if
    // all it's children are in the group.
    bool all_children_in_group = true;
    for (auto child : parent_header->get_children()) {
      bool in_group = false;
      for (auto f : group->m_individual_fields) {
        if (child.lock()->get_identifier() == f.second->get_header().get_identifier()) {
          in_group = true;
          break;
        }
      }
      if (not in_group) {
        all_children_in_group = false;
        break;
      }
    }

    // Set m_monolithic_field field
    if (all_children_in_group) {
      // If all children are in the group, this is the parent group and
      // we can just query the monolithic field and set here
      group->m_monolithic_field = parent_field;
    } else {
      // If children exist that aren't in this group, we need to get a
      // subfield of the parent group containing indices of group fields
      std::vector<int> ordered_subview_indices;
      for (auto fn : group->m_info->m_fields_names) {
        const auto parent_child_index =
          m_field_group_info.at(parent_field->name())->m_subview_idx.at(fn);
        ordered_subview_indices.push_back(parent_child_index);
      }
      std::sort(ordered_subview_indices.begin(),ordered_subview_indices.end());

      // Check that all subview indices are contiguous
      size_t span = ordered_subview_indices.back() - ordered_subview_indices.front() + 1;
      EKAT_REQUIRE_MSG (ordered_subview_indices.size()==span,
        "Error! Non-contiguous subview indices found in group \""+group_name+"\"\n");

      group->m_monolithic_field = std::make_shared<Field>(
        parent_field->subfield(
          group_name, parent_header->get_identifier().get_units(),
          group->m_info->m_subview_dim,
          ordered_subview_indices.front(), ordered_subview_indices.back()+1));

    }
  }

  return *group;
}

void FieldManager::
init_fields_time_stamp (const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot set initial time stamp until registration has completed.\n");

  for (auto& [grid_name, field_repo] : m_fields) {
    for (auto& field_it : field_repo) {
      field_it.second->get_header().get_tracking().update_time_stamp(t0);
    }
  }
}

void FieldManager::registration_begins ()
{
  // Update the state of the repo
  m_repo_state = RepoState::Open;
}

void FieldManager::registration_ends ()
{
  // This method is responsible of allocating the fields in the repo. The most delicate part is
  // the allocation of fields group, in the case where a monolithic allocation is required. If we
  // are not able to honor requests, we will error out. An example of a scenario where we can't
  // honor all requests is given by the three groups G1=(A,B), G2=(B,C), and G3=(A,C). Clearly,
  // only two of these groups can have contiguous allocation.

  // To understand how we can parse the groups to figure out if/how to accommodate all requests,
  // consider the following GR:  G1=(A,B,C), G2=(A,B,C,D,E), G3=(C,D), G4=(C,D,E,F), G5=((D,E,F,G).
  // The ordering of fields (A,B,C,D,E,F,G) clearly accommodates all GRs, though it might not
  // be obvious if one scrambled the fields within each GR.

  // For group requests on different grids, we consider the group members to be the union of all
  // fields registered to a particular group, over all grids. So if fields A and B are registered
  // as part of group G on grid1, and field C is registered as a part of G on grid2, both G on grid1
  // and G on grid2 will contain fields A,B,C, and we ensure that the subview index matches over grids.

  // This is the sequence of operations that allows us to establish if (and how) we can honor
  // all the requests:
  //  1) ensure all groups with a monolithic allocation contain the desired members. This means that
  //     we need to loop over GroupRequest (GR), and make sure there are fields registered in those
  //     groups (querying m_field_group_info info structs).
  //  2) Focus only on GR that require a monolithic allocation, discarding others.
  //     All the remaining group can simply "grab" individual fields later (and they
  //     can even grab some "individual" fields, and some fields that are slices of
  //     another group).
  //  3) Divide the groups found at step 2 into 'clusters'. A cluster is a collection
  //     of field groups such that each group in it share at least one field with
  //     at least another group in the cluster. If two groups are not in the same
  //     cluster, they cannot share any field. Each cluster of field groups can be
  //     processed individually, since allocations are completely independent.
  //     This allows to reduce the combination of groups to checks later.
  //     Note: to build a cluster, take any group not yet in a cluster. Then iteratively
  //     add any group that intersects the cluster, until no intersections are found.
  //  4) If there is no group in the cluster that contains all the fields of the cluster,
  //     create the cluster group C. In order to accommodate all requests for groups,
  //     we must be able to allocate C contiguously.
  //  5) For each cluster, call the function contiguous_superset from eamxx_utils.hpp
  //     (see that file for details). If the fcn fails to find an ordering of the cluster's
  //     field that accommodate all requests, it will return an empty list.
  //     Otherwise it will return the ordering of all fields in the cluster that allows all
  //     groups of the cluster to be a contiguous subset of C.
  //  Note: this can be done independent of grid since each grid will contain the same fields

  // Start by processing group request. This function checks that all group fields are properly
  // registered on the appropriate grid and the FieldGroupInfo is up to date. If group G is
  // requested, we make sure to register the union of all fields in group G on any grids they appear.
  pre_process_monolithic_group_requests();

  // Gather a list of groups whose fields require monolithic allocation
  std::list<std::string> groups_with_monolithic_allocation;
  for (auto& grid_it : m_grids_mgr->get_repo()) {
    for (const auto& greqs : m_group_requests.at(grid_it.second->name())) {
      for (const auto& r : greqs.second) {
        if (r.monolithic_alloc==MonolithicAlloc::Required) {
          // There's at least one request for this group to be allocated.
          groups_with_monolithic_allocation.push_back(greqs.first);
          break;
        }
      }
    }
  }
  // TODO: Is this still needed? Probably no.
  ::scream::sort(groups_with_monolithic_allocation);
  groups_with_monolithic_allocation.unique();

  // Homme currently wants qv to be the first tracer. We should be able to
  // modify Homme, to use something like qv_idx. However, that requires
  // extensive changes in Homme. Instead, we hack our way around this limitatin
  // (for now), and rearrange groups/fields so that we can expect qv to be the
  // first tracer.
  // Note: there are potentially two groups of tracers, "tracers" (which
  // homme advects) and "turbulence_advected_tracers" (which shoc advects).
  // "turbulence_advected_tracers" is guarenteed to be a subset of "tracers", so
  // moving qv first in "tracers" will do the same for "turbulence_advected_tracers".
  // Only homme requires qv first, so we only need to move qv if "tracers" exist.
  bool qv_must_come_first = false;
  if (ekat::contains(groups_with_monolithic_allocation,"tracers")
      and
      ekat::contains(m_field_group_info["tracers"]->m_fields_names, "qv")) {
    // Bring tracers to the front, so it will be processed first
    auto it = ekat::find(groups_with_monolithic_allocation,"tracers");
    std::swap(*it,groups_with_monolithic_allocation.front());

    // Adding the 'fake' group G=(qv) at the front of groups_with_monolithic_allocation ensures qv won't be put
    // in the middle of the tracers group. We use a highly unlikely group name, to avoid clashing
    // with a real group name. Later, after having found an global ordering for the tracers fields,
    // we will remove this group.
    groups_with_monolithic_allocation.push_front("__qv__");
    m_field_group_info.emplace("__qv__",std::make_shared<FieldGroupInfo>("__qv__"));
    m_field_group_info.at("__qv__")->m_fields_names.push_back("qv");
    qv_must_come_first = true;
  }

  // Do all the monolithic allocation stuff only if there are groups which require it.
  if (groups_with_monolithic_allocation.size()>0) {
    using namespace ShortFieldTagsNames;

    // A cluster is a list of names of groups in the cluster
    using cluster_type = std::list<std::string>;

    // Determine if two lists have elements in common (does not compute the intersection)
    auto intersect = [] (const std::list<ci_string>& lhs,
                        const std::list<ci_string>& rhs) -> bool {
      for (const auto& s : lhs) {
        if (ekat::contains(rhs,s)) {
          return true;
        }
      }
      return false;
    };

    std::list<cluster_type> clusters;
    std::list<std::string> added_to_a_cluster;
    while (added_to_a_cluster.size()<groups_with_monolithic_allocation.size()) {
      cluster_type c;
      auto first = groups_with_monolithic_allocation.begin();
      c.push_front(*first);
      groups_with_monolithic_allocation.erase(first);

      for (const auto& gn : groups_with_monolithic_allocation) {
        if (ekat::contains(added_to_a_cluster,gn)) {
          // This group has been added to a cluster already
          continue;
        }

        // Get the fields of this group
        const auto& fnames = m_field_group_info.at(gn)->m_fields_names;
        for (const auto& c_gn : c) {
          const auto& c_fnames = m_field_group_info.at(c_gn)->m_fields_names;
          if (intersect(fnames,c_fnames)) {
            // Ok, group gn intersects the cluser in at least one group (c_gn).
            // We add gn to the cluster, then break
            c.push_back(gn);
            added_to_a_cluster.push_back(gn);
            break;
          }
        }
      }

      clusters.emplace_back(c);
    }

    // Now we have clusters. For each cluster, build the list of lists, and call
    // the contiguous_superset utility.
    for (auto& cluster : clusters) {
      using LOL_t = std::list<std::list<ci_string>>;

      LOL_t groups_fields;
      for (const auto& gn : cluster) {
        groups_fields.push_back(m_field_group_info.at(gn)->m_fields_names);
        ::scream::sort(groups_fields.back());
      }

      auto cluster_ordered_fields = contiguous_superset(groups_fields);

      if (cluster_ordered_fields.size()==0) {
        // We were not able to accommodate all the allocation requests for the groups
        // in this cluster. We have to error out. But first, let's print some
        // information, so the developers/users can have a shot at fixing this
        // (e.g., by changing the request for monolithic allocation in some Atm Proc).
        std::cout << "Error! Field manager was not able to accommodate the following\n"
                  << "       requests for monolithically allocated groups:\n";
        for (const auto& gn : cluster) {
          std::cout << "   - " << gn << "\n";
        }
        std::cout << "       Consdier modifying the Atm Procs where these groups are requested.\n";

        EKAT_ERROR_MSG (" -- ERROR --\n");
      }

      // Ok, if we got here, it means we can allocate the cluster as a field, and then subview all the groups
      // Steps:
      //  - check if there's a group in the cluster containing all the fields. If yes, use that group
      //    name for the grouped field, otherwise make one up from the names of all groups in the cluster.
      //  - allocate the cluster field F.
      //  - loop over the groups in the cluster, and subview F at the proper (contiguous) indices.

      // WARNING: this lines should be removed if we move away from Homme's requirement that
      //          qv be the first tracer
      if (qv_must_come_first) {
        auto qv_it = ekat::find(cluster_ordered_fields,"qv");
        if (qv_it!=cluster_ordered_fields.end()) {
          // Check that qv comes first or last (if last, reverse the list). If not, error out.
          // NOTE: I *think* this should not happen, unless 'tracers' is a subgroup of a bigger group.
          EKAT_REQUIRE_MSG(qv_it==cluster_ordered_fields.begin() || std::next(qv_it,1)==cluster_ordered_fields.end(),
              "Error! The water vapor field has to be the first tracer, but it is not.\n");

          if (qv_it!=cluster_ordered_fields.begin()) {
            // Note: reversing the order of the fields preserves subgroups contiguity
            cluster_ordered_fields.reverse();
          }

          // Remove __qv__ from the cluster (if present)
          auto qv_gr_it = ekat::find(cluster,"__qv__");
          if (qv_gr_it!=cluster.end()) {
            cluster.erase(qv_gr_it);
          }
        }
      }

      // Check if there is a group with all the fields. Notice that it is enough to check
      // if any list in the LOL has the same length as cluster_ordered_fields.
      // If not, we will set cluster_name = $group1_name | $group2_name | ...
      // Note: cluster_name will be the name of the field allocating all fields in the cluster's groups
      std::string cluster_name;
      for (const auto& gn : cluster) {
        // Start building cluster_name by "or-ing" all gn's.
        if (gn==cluster.front()) {
          cluster_name = gn;
        } else {
          cluster_name += " | " + gn;
        }

        const auto& fnames = m_field_group_info.at(gn)->m_fields_names;
        if (fnames.size()==cluster_ordered_fields.size()) {
          // Found a group in the cluster that contains all fields.
          cluster_name = gn;
          break;
        }
      }

      // We must allocate this cluster field on every grid that
      // exists in the group requests contained in the cluster
      std::list<std::string> grids_in_cluster;
      for (const auto& gn : cluster) {
        for (auto grid_name: m_grids_mgr->get_grid_names()) {
          if (m_group_requests.at(grid_name).find(gn)!=m_group_requests.at(grid_name).end()) {
            grids_in_cluster.push_back(grid_name);
          }
        }
      }
      ::scream::sort(grids_in_cluster);
      grids_in_cluster.unique();

      for (auto cluster_grid_name : grids_in_cluster) {
        const auto& cluster_grid = m_grids_mgr->get_grid(cluster_grid_name);

        // Figure out the layout of the fields in this cluster,
        // and make sure they all have the same layout
        LayoutType lt = LayoutType::Invalid;
        FieldLayout f_layout = FieldLayout::invalid();
        for (const auto& fname : cluster_ordered_fields) {
          const auto& f = m_fields.at(cluster_grid_name).at(fname);
          const auto& id = f->get_header().get_identifier();
          if (lt==LayoutType::Invalid) {
          f_layout = id.get_layout();
          lt = f_layout.type();
          } else {
            EKAT_REQUIRE_MSG (lt==id.get_layout().type(),
                "Error! Found a group containing fields with different layouts.\n"
                "       Group name: " + cluster_name + "\n"
                "       Layout 1: " + e2str(lt) + "\n"
                "       Layout 2: " + e2str(id.get_layout().type()) + "\n");
          }
        }

        EKAT_REQUIRE_MSG(lt==LayoutType::Scalar2D || lt==LayoutType::Scalar3D,
            "Error! We can only monolithically allocate scalar fields. Found " + e2str(lt) + " fields instead.\n");

        FieldLayout c_layout = FieldLayout::invalid();
        if (lt==LayoutType::Scalar2D) {
          c_layout = cluster_grid->get_2d_vector_layout(cluster_ordered_fields.size());
        } else {
          c_layout = cluster_grid->get_3d_vector_layout(f_layout.tags().back()==LEV,cluster_ordered_fields.size());
        }

        // The units for the monolithic field are nondimensional, cause checking whether
        // all fields in the group have the same units so we can use those is too long and pointless.
        auto nondim = ekat::units::Units::nondimensional();

        // Allocate cluster field
        FieldIdentifier c_fid(cluster_name,c_layout,nondim,cluster_grid->name());
        register_field(c_fid);
        const auto& C = m_fields.at(cluster_grid_name).at(c_fid.name());

        // Scan all fields in this cluster, get their alloc props, and make sure
        // C can accommodate all of them. Also, make sure C can accommodate all
        // pack sizes of all GroupRequest of groups in this cluster.
        auto& C_ap = C->get_header().get_alloc_properties();
        for (const auto& fn : cluster_ordered_fields) {
          const auto& f = m_fields.at(cluster_grid_name).at(fn);
          C_ap.request_allocation(f->get_header().get_alloc_properties());
        }
        for (const auto& gn : cluster) {
          if (m_group_requests.at(cluster_grid_name).find(gn)!=m_group_requests.at(cluster_grid_name).end()) {
            for (const auto& req : m_group_requests.at(cluster_grid_name).at(gn)) {
              C_ap.request_allocation(req.pack_size);
            }
          }
        }

        // Allocate
        C->allocate_view();

        // Note: as of 02/2021, idim should *always* be 1, but we store it just in case,
        //       to avoid bugs in the future.
        const auto& C_tags = C->get_header().get_identifier().get_layout().tags();
        const int idim = std::distance(C_tags.begin(),ekat::find(C_tags,CMP));
        // See below where we assume idim is 1
        EKAT_REQUIRE_MSG(idim==1, "Error! idim is assumed to be 1 in FieldManager::registration_ends().\n");

        // Create all individual subfields
        for (const auto& fn : cluster_ordered_fields) {
          const auto pos = ekat::find(cluster_ordered_fields,fn);
          const auto idx = std::distance(cluster_ordered_fields.begin(),pos);

          const auto& f = m_fields.at(cluster_grid_name).at(fn);
          const auto& fid = f->get_header().get_identifier();
          EKAT_REQUIRE_MSG (fid.get_units()!=ekat::units::Units::invalid(),
              "Error! A field was registered without providing valid units.\n"
              "  - field id: " + fid.get_id_string() + "\n");
          auto fi = C->subfield(fn,fid.get_units(),idim,idx);

          // Overwrite existing field with subfield
          *f = fi;
        }
      }

      if (ekat::contains(cluster,"__qv__")) {
        // Erase the 'fake group' we added (to guarantee qv would be first/last in the tracers)
        m_field_group_info.erase(m_field_group_info.find("__qv__"));
        cluster.erase(ekat::find(cluster,"__qv__"));
      }

      // Now, update the group info of all the field groups in the cluster
      for (const auto& gn : cluster) {
        auto& info = *m_field_group_info.at(gn);
        const auto n = info.size();

        // Find the first field of this group in the ordered cluster names.
        auto first = std::find_first_of(cluster_ordered_fields.begin(),cluster_ordered_fields.end(),
                                        info.m_fields_names.begin(),info.m_fields_names.end());
        auto last = std::next(first,n);

        // Some sanity checks: info.m_fields_names should be a rearrangement of what is
        // in cluster_odered_fields[first,last)
        EKAT_REQUIRE_MSG(first!=cluster_ordered_fields.end(),
            "Error! Could not find any field of this group in the ordered cluster.\n"
            "       Group name: " + gn + "\n");
        EKAT_REQUIRE_MSG(std::distance(last,cluster_ordered_fields.end())>=0,
            "Error! Something went wrong while looking for fields of this group in the ordered cluster.\n"
            "       Group name: " + gn + "\n");
        EKAT_REQUIRE_MSG(std::is_permutation(first,last,info.m_fields_names.begin(),info.m_fields_names.end()),
            "Error! Something went wrong while looking for fields of this group in the ordered cluster.\n"
            "       Group name: " + gn + "\n");

        // Update the list of fields in the group info, mark it as monolithically
        // allocated, and update the subfield indices too.
        info.m_fields_names.clear();
        for (auto it=first; it!=last; ++it) {
          info.m_fields_names.push_back(*it);
          info.m_subview_dim = 1; // Assumption is checked above
          info.m_subview_idx [*it] = std::distance(cluster_ordered_fields.begin(),it);
        }
        info.m_monolithic_allocation = true;

        // The subview indices will need to be corrected in the case of this group being
        // itself a subfile of the cluster field. We are guarenteed that the indices
        // are contiguous, so just calculate the min idx and subtract from all fields indices.
        int min_idx = std::numeric_limits<int>::max();
        for (auto& it : info.m_subview_idx) {
          if (it.second<min_idx) {
            min_idx = it.second;
          }
        }
        for (auto& it : info.m_subview_idx) {
          it.second -= min_idx;
        }
      }
    }
  }

  for (auto grid_name : m_grids_mgr->get_grid_names()) {
    for (auto& it : m_fields.at(grid_name)) {
      if (it.second->is_allocated()) {
        // If the field has been already allocated, then it was in a bunlded group, so skip it.
        continue;
      }
      // A brand new field. Allocate it
      it.second->allocate_view();
    }

    for (const auto& it : m_field_group_info) {
      if (m_group_requests.at(grid_name).find(it.first)!=m_group_requests.at(grid_name).end()) {
        // Get fields in this group
        auto group_info = it.second;
        const auto& fnames = group_info->m_fields_names;
        for (const auto& fn : fnames) {
          // Update the field tracking
          m_fields.at(grid_name).at(fn)->get_header().get_tracking().add_to_group(group_info);
        }
      }
    }
  }

  // Prohibit further registration of fields to this repo
  m_repo_state = RepoState::Closed;
}

void FieldManager::clean_up() {
  // Clean field map
  m_fields.clear();

  // Clean group info
  m_field_group_info.clear();

  // Reset repo state
  m_repo_state = RepoState::Clean;
}

void FieldManager::clean_up(const std::string& grid_name) {
  // Clear the maps
  m_fields[grid_name].clear();
}

void FieldManager::add_field (const Field& f) {
  const auto& grid_name = f.get_header().get_identifier().get_grid_name();

  // This method has a few restrictions on the input field.
  EKAT_REQUIRE_MSG(m_grids_mgr->has_grid(grid_name),
    "Error! Attempting to add_field on grid not in the FM's grids manager:\n"
    "  - Field name: " + f.name() + "\n"
    "  - Grid:       " + grid_name + "\n"
    "  - Grids stored by FM: " + m_grids_mgr->print_available_grids() + "\n");
  EKAT_REQUIRE_MSG (not (m_repo_state==RepoState::Open),
      "Error! The method 'add_field' can only be called on a closed repo.\n");
  EKAT_REQUIRE_MSG (f.is_allocated(),
      "Error! The method 'add_field' requires the input field to be already allocated.\n");
  EKAT_REQUIRE_MSG (m_grids_mgr->get_grid(grid_name)->is_valid_layout(f.get_header().get_identifier().get_layout()),
      "Error! Input field to 'add_field' has a layout not compatible with the stored grid.\n"
      "  - input field name : " + f.name() + "\n"
      "  - field manager grid: " + grid_name + "\n"
      "  - input field layout:   " + f.get_header().get_identifier().get_layout().to_string() + "\n");
  EKAT_REQUIRE_MSG (not has_field(f.name(), grid_name),
      "Error! The method 'add_field' requires the input field to not be already existing.\n"
      "  - field name: " + f.get_header().get_identifier().name() + "\n");
  EKAT_REQUIRE_MSG (f.get_header().get_tracking().get_groups_info().size()==0 ||
                    m_group_requests.at(grid_name).size()==0,
      "Error! When calling 'add_field', one of the following must be true:\n"
      "  - the input field is not be part of any group,\n"
      "  - there were no group requests for this field manager.\n"
      "The reason for this is that otherwise we *might* have missed some inclusion dependency\n"
      "when we allocated the fields for one of those groups.\n");

  // All good, add the field to the repo
  m_fields[grid_name][f.get_header().get_identifier().name()] = std::make_shared<Field>(f);

  if (m_repo_state==RepoState::Clean) m_repo_state = RepoState::Closed;
}

void FieldManager::pre_process_monolithic_group_requests () {
  // For each group, loop over all fields in the group and register
  // on to each grid the group is requested on (if necessary)
  for (auto [group_name, group_info] : m_field_group_info) {
    if (not group_info->m_monolithic_allocation) continue;
    // Gather all grids in this group. We need to check
    // both grids that already exist in the group and grids that
    // at one point were requested by the group.
    std::list<ci_string> grids_in_group;
    for (auto field_name : group_info->m_fields_names) {
      for (auto grid_name : group_info->m_grid_registered.at(field_name)) {
        if (not ekat::contains(grids_in_group, grid_name)) grids_in_group.push_back(grid_name);
      }
    }
    for (auto grid_name : group_info->m_requested_grids) {
      if (not ekat::contains(grids_in_group, grid_name)) grids_in_group.push_back(grid_name);
    }

    // Register fields on all grids
    for (auto field_name : group_info->m_fields_names) {
      for (auto grid_name : grids_in_group) {
        // Field already registered on grid, nothing to do
        if (has_field(field_name, grid_name)) continue;

        // Find the grid which registered this field.
        // Note: We require *exactly one* grid has registered the field
        std::string registered_grid = "";
        for (auto other_grid_name : group_info->m_grid_registered.at(field_name)) {
          if (other_grid_name == grid_name) continue;

          // Check that we haven't found 2 other grids which registered this field
          EKAT_REQUIRE_MSG(registered_grid == "" or registered_grid == other_grid_name,
            "Error! Before registration_end(), FieldManager requires a field within a group to be registered on a single grid.\n"
            "  - Field name: " + field_name + "\n"
            "  - Group name: " + group_name + "\n"
            "  - Grids registered: " + other_grid_name + "," + registered_grid + "\n");

            registered_grid = other_grid_name;
        }

        // Get a FID for this grid by asking for an equivalent layout
        // to the layout on the src grid
        const auto src_fid = m_fields.at(registered_grid).at(field_name)->get_header().get_identifier();
        const auto fl = m_grids_mgr->get_grid(grid_name)->equivalent_layout(src_fid.get_layout());
        FieldIdentifier fid(field_name, fl, src_fid.get_units(), grid_name);

        // Register the field for each group req
        for (auto greq : m_group_requests.at(grid_name).at(group_name)) {
          FieldRequest req(fid, greq.name, greq.pack_size);
          register_field(req);
        }
      }
    }
  }
}

} // namespace scream
