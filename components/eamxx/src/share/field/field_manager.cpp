#include "share/field/field_manager.hpp"

namespace scream
{

FieldManager::
FieldManager (const grid_ptr_type& grid)
  : m_repo_state (RepoState::Clean)
  , m_grid       (grid)
{
  EKAT_REQUIRE_MSG (m_grid!=nullptr,
      "Error! Input grid pointer is not valid.");
}

void FieldManager::register_field (const FieldRequest& req)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Clean,
      "Error! Repo state is not 'Open' yet. You must call registration_begins() first.\n");
  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Closed,
      "Error! Repo state is not 'Open' anymore. You already called registration_ends().\n");

  const auto& id = req.fid;

  // Make sure the grid name from the id matches the name of m_grid
  EKAT_REQUIRE_MSG(id.get_grid_name()==m_grid->name(),
      "Error! Input field request identifier stores a different grid name than the FM:\n"
      "         - input id grid name: " + id.get_grid_name() + "\n"
      "         - stored grid name:   " + m_grid->name() + "\n");

  // Get or create the new field
  if (!has_field(id.name())) {
    EKAT_REQUIRE_MSG (id.data_type()==field_valid_data_types().at<Real>(),
        "Error! While refactoring, we only allow the Field data type to be Real.\n"
        "       If you're done with refactoring, go back and fix things.\n");
    m_fields[id.name()] = std::make_shared<Field>(id);
  } else {
    // Make sure the input field has the same layout and units as the field already stored.
    // TODO: this is the easiest way to ensure everyone uses the same units.
    //       However, in the future, we *may* allow different units, providing
    //       the users with conversion routines perhaps.
    const auto id0 = m_fields[id.name()]->get_header().get_identifier();
    EKAT_REQUIRE_MSG(id.get_units()==id0.get_units(),
        "Error! Field '" + id.name() + "' already registered with different units:\n"
        "         - input field units:  " + to_string(id.get_units()) + "\n"
        "         - stored field units: " + to_string(id0.get_units()) + "\n"
        "       Please, check and make sure all atmosphere processes use the same units.\n");

    EKAT_REQUIRE_MSG(id.get_layout()==id0.get_layout(),
        "Error! Field '" + id.name() + "' already registered with different layout:\n"
        "         - input id:  " + id.get_id_string() + "\n"
        "         - stored id: " + id0.get_id_string() + "\n"
        "       Please, check and make sure all atmosphere processes use the same layout for a given field.\n");
  }

  if (req.subview_info.dim_idx>=0) {
    // This is a request for a subfield. Store request info, so we can correctly set up
    // the subfield at the end of registration_ends() call
    m_subfield_requests.emplace(id.name(),req);
  } else {
    // Make sure the field can accommodate the requested value type
    m_fields[id.name()]->get_header().get_alloc_properties().request_allocation(req.pack_size);
  }

  // Finally, add the field to the given groups
  // Note: we do *not* set the group info struct in the field header yet.
  //       we will do that when we end the registration phase.
  //       The reason is that at registration_ends we will know *all* the groups
  //       that each field belongs to.
  for (const auto& group_name : req.groups) {
    // Get group (and init ptr, if necessary)
    auto& group_info = m_field_groups[group_name];
    if (group_info==nullptr) {
      group_info = std::make_shared<group_info_type>(group_name);
    }

    // Add the field name to the list of fields belonging to this group
    if (ekat::find(group_info->m_fields_names,id.name())==group_info->m_fields_names.end()) {
      group_info->m_fields_names.push_back(id.name());
    }

    // In order to simplify the logic in pre_process_group_requests,
    // add a "trivial" GroupRequest for this group, meaning no bundling,
    // pack size 1, and no derivation from another group. This ensure that each group
    // in m_field_groups also appears in the m_group_requests map.
    register_group(GroupRequest(group_name,m_grid->name()));
  }
}

void FieldManager::register_group (const GroupRequest& req)
{
  EKAT_REQUIRE_MSG (req.grid==m_grid->name(),
      "Error! Input GroupRequest grid does not match the grid of this FieldManager.\n"
      "       Request grid: " + req.grid + "\n"
      "       Stored grid:  " + m_grid->name() + "\n");

  // Groups have to be handled once registration is over, so for now simply store the request,
  // and create an empty group info
  m_group_requests[req.name].insert(req);

  m_field_groups.emplace(req.name,std::make_shared<group_info_type>(req.name));
}

void FieldManager::
add_to_group (const std::string& field_name, const std::string& group_name)
{
  EKAT_REQUIRE_MSG (m_repo_state==RepoState::Closed,
      "Error! You cannot call 'add_to_group' until after 'registration_ends' has been called.\n");
  auto& group = m_field_groups[group_name];
  if (not group) {
    group = std::make_shared<FieldGroupInfo>(group_name);
  }
  EKAT_REQUIRE_MSG (not group->m_bundled,
      "Error! Cannot add fields to a group that is bundled.\n"
      "   group name: " + field_name + "\n");

  EKAT_REQUIRE_MSG (has_field(field_name),
      "Error! Cannot add field to group, since the field is not present in this FieldManager.\n"
      "   field name: " + field_name + "\n"
      "   group name: " + group_name + "\n");

  group->m_fields_names.push_back(field_name);
  auto& ft = get_field(field_name).get_header().get_tracking();
  ft.add_to_group(group);
}

bool FieldManager::has_field (const identifier_type& id) const
{
  return has_field(id.name()) && m_fields.at(id.name())->get_header().get_identifier()==id;
}

Field FieldManager::get_field (const identifier_type& id) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(id);
  EKAT_REQUIRE_MSG(ptr!=nullptr,
      "Error! Field identifier '" + id.get_id_string() + "' not found.\n");
  return *ptr;
}

Field FieldManager::get_field (const std::string& name) const {

  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(name);
  EKAT_REQUIRE_MSG(ptr!=nullptr, "Error! Field " + name + " not found.\n");
  return *ptr;
}

FieldGroup FieldManager::
get_field_group (const std::string& group_name) const
{
  // Sanity checks
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");
  EKAT_REQUIRE_MSG (has_group(group_name),
      "Error! Field group '" + group_name + "' not found.\n");

  // Create an empty group
  FieldGroup group(group_name);

  // Set the info in the group
  group.m_info = m_field_groups.at(group_name);

  // IMPORTANT: this group might have been created as a "Copy" (see field_request.hpp)
  //            of another group. In this case, we only allocate the bundled group,
  //            and we don't create individual subfields. Users can still extract
  //            subfields on the fly, if need be, by using the names/idx data in m_info.
  auto it = m_group_requests.find(group_name);
  bool allocate_subfields = true;
  if (it!=m_group_requests.end()) {
    for (const auto& req : it->second) {
      if (req.derived_type==DerivationType::Copy) {
        allocate_subfields = false;
        break;
      }
    }
  }

  if (allocate_subfields) {
    // Find all the fields and set them in the group
    for (const auto& fname : group.m_info->m_fields_names) {
      auto f = get_field_ptr(fname);
      group.m_fields[fname] = f;
    }
  }

  // Fetch the bundle field (if bundled)
  if (group.m_info->m_bundled) {
    if (allocate_subfields) {
      // We can get the parent from any of the fields in the group.
      auto p = group.m_fields.begin()->second->get_header().get_parent().lock();
      EKAT_REQUIRE_MSG(p!=nullptr,
          "Error! A field belonging to a bundled field group is missing its 'parent'.\n");

      const auto& p_id = p->get_identifier();
      group.m_bundle = get_field_ptr(p_id);
    } else {
      // This group was copied, and therefore we know we created the bundled field with
      // the same name as the group itself. So simply fetch that field
      group.m_bundle = get_field_ptr(group_name);
    }
  }

  return group;
}

void FieldManager::
init_fields_time_stamp (const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot set initial time stamp until registration has completed.\n");

  for (auto it : m_fields) {
    it.second->get_header().get_tracking().update_time_stamp(t0);
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
  // the allocation of fields group, in the case where bundling is requested. In particular,
  // we want to try to honor as many requests for bundling as possible. If we can't accommodate
  // all the GroupRequest's that have bundling either Required or Preferred, we will try again
  // by considering only those with Required bundling (see GroupRequest for detail on those
  // values). If we are still not able to honor requests, we will error out. An example of a
  // scenario where we can't honor all requests is given by the three groups G1=(A,B), G2=(B,C),
  // and G3=(A,C). Clearly, only two of these groups can have contiguous allocation.

  // To understand how we can parse the groups to figure out if/how to accommodate all requests,
  // consider the following GR:  G1=(A,B,C), G2=(A,B,C,D,E), G3=(C,D), G4=(C,D,E,F), G5=((D,E,F,G).
  // The ordering of fields (A,B,C,D,E,F,G) clearly accommodates all GRs, though it might not
  // be obvious if one scrambled the fields within each GR.

  // This is the sequence of operations that allows us to establish if (and how) we can honor
  // all the requests:
  //  1) ensure all groups contain the desired members. This means that we need to
  //     loop over GroupRequest (GR), and make sure there are fields registered in those
  //     groups (querying m_field_groups info structs). If a GR is derived from another
  //     GR, like a 'Subset' group (see GR header for details), we make sure the group
  //     is in m_field_groups (if not, add it), and contains all the proper fields.
  //  2) Focus only on GR that require (or prefer) a bundled group, discarding others.
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
  //     create the cluster group C. In order to accommodate all requests for bundled groups,
  //     we must be able to allocate C bundled.
  //  5) For each cluster, call the function contiguous_superset from scream_utils.hpp
  //     (see that file for details). If the fcn fails to find an ordering of the cluster's
  //     field that accommodate all bundled requests, it will return an empty list.
  //     Otherwise it will return the ordering of all fields in the cluster that allows all
  //     groups of the cluster to be a contiguous subset of C.
  //  6) If step 4 fails for a cluster, remove from the cluster the groups whose bundling
  //     is only "Preferred", and re-try.
  //     Notes:
  //      - eliminating 1+ groups from the cluster may actually "disconnect" the cluster,
  //        into 2+ separate clusters, which can be treated separately. For simplicity,
  //        we won't re build the cluster. It might make calling the function
  //        contiguous_superset a bit more expensive, but it won't alter the existence
  //        of a contiguous ordering.
  //      - it is possible that, say, removing G1 still doesn't yield a "bundle-able"
  //        cluster, but removing G2 does. In general, we should try to remove groups
  //        one at a time, then two at a time, then three at a time,... until we
  //        reach a point where the remaining groups can all be bundled. This is overly
  //        complicated, so if 4 fails, we simply start removing groups with "Preferred"
  //        bundling until 4 succeeds or we run out of groups with bundling=Preferred
  //        to remove.
  //


  // Start by processing group request. This function will ensure that, if there's a
  // request for group A that depends on the content of group B, the FieldGroupInfo
  // for group A is updated to contain the correct fields, based on the content
  // of group B. It also checks that the requests are not inconsistent.
  pre_process_group_requests ();

  // Gather a list of groups to be bundled
  // NOTE: copied groups are always created bundled, but in a second phase,
  //       without creating individual subfields.
  std::list<std::string> groups_to_bundle, copied_groups;
  for (const auto& greqs : m_group_requests) {
    for (const auto& r : greqs.second) {
      if (r.derived_type!=DerivationType::Copy && r.bundling!=Bundling::NotNeeded) {
        // There's at least one request for this group to be bunlded.
        groups_to_bundle.push_back(greqs.first);
      } else if (r.derived_type==DerivationType::Copy) {
        copied_groups.push_back(greqs.first);
      }
    }
  }
  ::scream::sort(groups_to_bundle);
  ::scream::sort(copied_groups);
  groups_to_bundle.unique();
  copied_groups.unique();

  // Homme currently wants qv to be the first tracer. We should be able to
  // modify Homme, to use something like qv_idx. However, that requires
  // extensive changes in Homme. Instead, we hack our way around this limitatin
  // (for now), and rearrange groups/fields so that we can expect qv to be the
  // first tracer.
  bool qv_must_come_first = false;
  if (has_field("qv")) {
    auto it = ekat::find(groups_to_bundle,"tracers");
    if (it!=groups_to_bundle.end()) {
      // Bring tracers to the front, so it will be processed first
      std::swap(*it,groups_to_bundle.front());

      // Adding the 'fake' group G=(qv) at the front of groups_to_bundle ensures qv won't be put
      // in the middle of the tracers group. We use a highly unlikely group name, to avoid clashing
      // with a real group name. Later, after having found an global ordering for the tracers fields,
      // we will remove this group.
      groups_to_bundle.push_front("__qv__");
      m_field_groups.emplace("__qv__",std::make_shared<group_info_type>("__qv__"));
      m_field_groups.at("__qv__")->m_fields_names.push_back("qv");
      qv_must_come_first = true;
    }
  }

  // Do all the bundling stuff only if there are groups do bundle at all.
  if (groups_to_bundle.size()>0) {
    using namespace ShortFieldTagsNames;

    // A cluster is a pair <cluster_name,list of names of groups in the cluster>
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
    while (added_to_a_cluster.size()<groups_to_bundle.size()) {
      cluster_type c;
      auto first = groups_to_bundle.begin();
      c.push_front(*first);
      groups_to_bundle.erase(first);

      for (const auto& gn : groups_to_bundle) {
        if (ekat::contains(added_to_a_cluster,gn)) {
          // This group has been added to a cluster already
          continue;
        }

        // Get the fields of this group
        const auto& fnames = m_field_groups.at(gn)->m_fields_names;
        for (const auto& c_gn : c) {
          const auto& c_fnames = m_field_groups.at(c_gn)->m_fields_names;
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
    // the contiguous_superset method.
    for (auto& cluster : clusters) {
      using LOL_t = std::list<std::list<ci_string>>;

      LOL_t groups_fields;
      for (const auto& gn : cluster) {
        groups_fields.push_back(m_field_groups.at(gn)->m_fields_names);
        ::scream::sort(groups_fields.back());
      }

      auto cluster_ordered_fields = contiguous_superset(groups_fields);
      while (cluster_ordered_fields.size()==0) {
        // Try to see if there's a group we can remove, that is, a group
        // for which bundling is only Preferred. If there's no such group,
        // we break the loop, and we will crap out.

        std::list<std::string>::iterator it = cluster.end();
        for (const auto& gn : cluster) {
          const auto& reqs = m_group_requests.at(gn);
          bool remove_this = true;
          for (const auto& r : reqs) {
            if (r.bundling==Bundling::Required) {
              // Can't remove this group, cause at least one request "requires" bundling
              remove_this = false;
              break;
            }
          }

          if (remove_this) {
            // It's ok to try and remove this group from the cluster
            it = ekat::find(cluster,gn);
            break;
          }
        }

        if (it==cluster.end()) {
          // We were not able to find a group that can be removed from the cluster,
          // meaning that all groups in the cluster are "Required" to be bunlded.
          // Time to quit.
          break;
        }

        // Ok, let's remove group *it, and try again
        // Note: we have to remove the group name from gnames, but we also have to
        //       remove its list of fields from group_fields
        auto pos = std::distance(cluster.begin(),it);
        cluster.erase(std::next(cluster.begin(),pos));
        groups_fields.erase(std::next(groups_fields.begin(),pos));
        cluster_ordered_fields = contiguous_superset(groups_fields);
      }

      if (cluster_ordered_fields.size()==0) {
        // We were not able to accommodate all the requests bundling the groups
        // in this cluster. We have to error out. But first, let's print some
        // information, so the developers/users can have a shot at fixing this
        // (e.g., by changing the request for bundling in some Atm Proc).
        std::cout << "Error! Field manager on grid " << m_grid->name() <<  " was not able to accommodate\n"
                  << "       the following requests for bundled groups:\n";
        for (const auto& gn : cluster) {
          std::cout << "   - " << gn << "\n";
        }
        std::cout << "       Consdier modifying the Atm Procs where these groups are requested.\n"
                  << "       For instance, you may ask for 'Preferred' bundling, rather than 'Required'\n";

        EKAT_ERROR_MSG (" -- ERROR --\n");
      }

      // Ok, if we got here, it means we can allocate the cluster as a field, and then subview all the groups
      // Steps:
      //  - check if there's a group in the cluster containing all the fields. If yes, use that group
      //    name for the bundled field, otherwise make one up from the names of all groups in the cluster.
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
      // Note: cluster_name will be the name of the field bundling all fields in the cluster's groups
      std::string cluster_name;
      for (const auto& gn : cluster) {
        // Start building cluster_name by "or-ing" all gn's.
        if (gn==cluster.front()) {
          cluster_name = gn;
        } else {
          cluster_name += " | " + gn;
        }

        const auto& fnames = m_field_groups.at(gn)->m_fields_names;
        if (fnames.size()==cluster_ordered_fields.size()) {
          // Found a group in the cluster that contains all fields.
          cluster_name = gn;
          break;
        }
      }

      // Figure out the layout of the fields in this cluster,
      // and make sure they all have the same layout
      LayoutType lt = LayoutType::Invalid;
      std::shared_ptr<const FieldLayout> f_layout;
      for (const auto& fname : cluster_ordered_fields) {
        const auto& f = m_fields.at(fname);
        const auto& id = f->get_header().get_identifier();
        if (lt==LayoutType::Invalid) {
         f_layout = id.get_layout_ptr();
         lt = get_layout_type(f_layout->tags());
        } else {
          EKAT_REQUIRE_MSG (lt==get_layout_type(id.get_layout().tags()),
              "Error! Found a group to bundle containing fields with different layouts.\n"
              "       Group name: " + cluster_name + "\n"
              "       Layout 1: " + e2str(lt) + "\n"
              "       Layout 2: " + e2str(get_layout_type(id.get_layout().tags())) + "\n");
        }
      }

      EKAT_REQUIRE_MSG(lt==LayoutType::Scalar2D || lt==LayoutType::Scalar3D,
          "Error! We can only bundle scalar fields. Found " + e2str(lt) + " fields instead.\n");

      FieldLayout c_layout = FieldLayout::invalid();
      if (lt==LayoutType::Scalar2D) {
        c_layout = m_grid->get_2d_vector_layout(CMP,cluster_ordered_fields.size());
      } else {
        c_layout = m_grid->get_3d_vector_layout(f_layout->tags().back()==LEV,CMP,cluster_ordered_fields.size());
      }

      // The units for the bundled field are nondimensional, cause checking whether
      // all fields in the bundle have the same units so we can use those is too long and pointless.
      auto nondim = ekat::units::Units::nondimensional();

      // Allocate cluster field
      FieldIdentifier c_fid(cluster_name,c_layout,nondim,m_grid->name());
      register_field(c_fid);
      const auto& C = m_fields.at(c_fid.name());

      // Scan all fields in this cluster, get their alloc props, and make sure
      // C can accommodate all of them. Also, make sure C can accommodate all
      // pack sizes of all GroupRequest of groups in this cluster.
      auto& C_ap = C->get_header().get_alloc_properties();
      for (const auto& fn : cluster_ordered_fields) {
        const auto& f = m_fields.at(fn);
        C_ap.request_allocation(f->get_header().get_alloc_properties());
      }
      for (const auto& gn : cluster) {
        for (const auto& req : m_group_requests.at(gn)) {
          C_ap.request_allocation(req.pack_size);
        }
      }

      // Allocate
      C->allocate_view();

      // Note: as of 02/2021, idim should *always* be 1, but we store it just in case,
      //       to avoid bugs in the future.
      const auto& C_tags = C->get_header().get_identifier().get_layout().tags();
      const int idim = std::distance(C_tags.begin(),ekat::find(C_tags,CMP));

      if (ekat::contains(cluster,"__qv__")) {
        // Erase the 'fake group' we added (to guarantee qv would be first/last in the tracers)
        m_field_groups.erase(m_field_groups.find("__qv__"));
        cluster.erase(ekat::find(cluster,"__qv__"));
      }

      // Create all individual subfields
      for (const auto& fn : cluster_ordered_fields) {
        const auto pos = ekat::find(cluster_ordered_fields,fn);
        const auto idx = std::distance(cluster_ordered_fields.begin(),pos);

        const auto& f = m_fields.at(fn);
        const auto& fid = f->get_header().get_identifier();
        EKAT_REQUIRE_MSG (fid.get_units()!=ekat::units::Units::invalid(),
            "Error! A field was registered without providing valid units.\n"
            "  - field id: " + fid.get_id_string() + "\n");
        auto fi = C->subfield(fn,fid.get_units(),idim,idx);

        // Overwrite existing field with subfield
        *f = fi;
      }

      // Now, update the group info of all the field groups in the cluster
      for (const auto& gn : cluster) {
        auto& info = *m_field_groups.at(gn);
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

        // Update the list of fields in the group info, mark it as bundled,
        // and update the subfield indices too.
        info.m_fields_names.clear();
        for (auto it=first; it!=last; ++it) {
          info.m_fields_names.push_back(*it);
          info.m_subview_dim = idim;
          info.m_subview_idx [*it] = std::distance(cluster_ordered_fields.begin(),it);
        }
        info.m_bundled = true;
      }
    }
  }

  // Now that all bundled fields have been taken care of, we can
  //  - create copied groups (as a bundled field) and update their info
  //  - allocate remaining fields (not part of a bundled field)
  //  - update the tracking of each field, by adding the group info of all
  //    the groups the field belongs to.
  for (const auto& gname : copied_groups) {
    using namespace ShortFieldTagsNames;
    auto& info  = *m_field_groups.at(gname);
    const int   size  = info.size();

    // The units for the bundled field are nondimensional, cause checking whether
    // all fields in the bundle have the same units so we can use those is too long and pointless.
    auto nondim = ekat::units::Units::nondimensional();

    // Take any of the fields in the group to be copied, and determine
    // whether they are 2d or 3d.
    auto f1 = m_fields.at(info.m_fields_names.front());
    auto f1_layout = f1->get_header().get_identifier().get_layout();
    auto lt = get_layout_type(f1_layout.tags());
    FieldLayout g_layout = FieldLayout::invalid();
    if (lt==LayoutType::Scalar2D) {
      g_layout = m_grid->get_2d_vector_layout(CMP,size);
    } else {
      bool mid = f1_layout.tags().back()==LEV;
      g_layout = m_grid->get_3d_vector_layout(mid,CMP,size);
    }

    FieldIdentifier g_fid(gname,g_layout,nondim,m_grid->name());
    register_field(g_fid);
    const auto& G = m_fields.at(g_fid.name());

    // Loop over the requests, and set the alloc props
    auto& G_ap = G->get_header().get_alloc_properties();
    for (const auto& req : m_group_requests.at(gname)) {
      G_ap.request_allocation(req.pack_size);
    }
    G->allocate_view();

    // Now, update the group info of the copied group, by setting the
    // correct subview_idx, in case the user wants to extract the
    // individual fields as Field structs.
    // NOTE: individual fields are *not* registered in the repo,
    //       but the user *can* create subfields, if needed.
    //       E.g., user might create subfields for remapping purposes.
    // Note: as of 08/2021, idim should *always* be 1, but we store it just in case,
    //       to avoid bugs in the future.
    const auto& G_tags = g_layout.tags();
    const int idim = std::distance(G_tags.begin(),ekat::find(G_tags,CMP));
    for (int i=0; i<size; ++i) {
      const auto& fn = *std::next(info.m_fields_names.begin(),i);
      info.m_subview_dim = idim;
      info.m_subview_idx [fn] = i;
    }
    info.m_bundled = true;
  }

  for (auto& it : m_fields) {
    if (it.second->is_allocated()) {
      // If the field has been already allocated, then it was in a bunlded group, so skip it.
      continue;
    }
    // A brand new field. Allocate it
    it.second->allocate_view();
  }

  for (const auto& it : m_field_groups) {
    if (ekat::contains(copied_groups,it.first)) {
      // Copied groups use the same field names as the group they copied,
      // but those fields do not belong to the copy group.
      continue;
    }
    // Get fields in this group
    auto group_info = it.second;
    const auto& fnames = group_info->m_fields_names;
    for (const auto& fn : fnames) {
      // Update the field tracking
      m_fields.at(fn)->get_header().get_tracking().add_to_group(group_info);
    }
  }

  // Prohibit further registration of fields
  m_repo_state = RepoState::Closed;
}

void FieldManager::clean_up() {
  // Clear the maps
  m_fields.clear();
  m_field_groups.clear();

  // Reset repo state
  m_repo_state = RepoState::Clean;
}

void FieldManager::add_field (const Field& f) {
  // This method has a few restrictions on the input field.
  EKAT_REQUIRE_MSG (m_repo_state==RepoState::Closed,
      "Error! The method 'add_field' can only be called on a closed repo.\n");
  EKAT_REQUIRE_MSG (f.is_allocated(),
      "Error! The method 'add_field' requires the input field to be already allocated.\n");
  EKAT_REQUIRE_MSG (f.get_header().get_identifier().get_grid_name()==m_grid->name(),
      "Error! Input field to 'add_field' is defined on a grid different from the one stored.\n"
      "  - field manager grid: " + m_grid->name() + "\n"
      "  - input field grid:   " + f.get_header().get_identifier().get_grid_name() + "\n");
  EKAT_REQUIRE_MSG (not has_field(f.get_header().get_identifier().name()),
      "Error! The method 'add_field' requires the input field to not be already existing.\n"
      "  - field name: " + f.get_header().get_identifier().name() + "\n");
  EKAT_REQUIRE_MSG (f.get_header().get_tracking().get_groups_info().size()==0,
      "Error! The method 'add_field' requires the input field to not be part of any group.\n");

  // All good, add the field to the repo
  m_fields[f.get_header().get_identifier().name()] = std::make_shared<Field>(f);
}

std::shared_ptr<Field>
FieldManager::get_field_ptr (const identifier_type& id) const {
  auto it = m_fields.find(id.name());
  if (it!=m_fields.end() && it->second->get_header().get_identifier()==id) {
    return it->second;
  }
  return nullptr;
}

std::shared_ptr<Field>
FieldManager::get_field_ptr (const std::string& name) const {
  auto it = m_fields.find(name);
  return it==m_fields.end() ? nullptr : it->second;
}

void FieldManager::pre_process_group_requests () {
  // GroupRequests that are formulated in terms of another group (i.e., where
  // req.derived_type!=DerivationType::None), need a preprocessing step.
  // This step needs to do two things: 1) make sure the group contains all
  // the fields it needs, and 2) make sure the different group requests are
  // consistent and valid. Point 1 is necessary, since nobody may have added
  // any field to the group. E.g., if group B is a subset of group A={f1,f2,f3}
  // with f1 removed, we need to make sure B={f2,f3} (might be nobody registered
  // f2 and f3 as part of B). Point 2 is to make sure we don't ask for things
  // that are either too complicated to handle (and likely never going to be needed)
  // or inconsistent. E.g., if group A is requested as a "Copy" of group B,
  // and group B is requested as a "Copy" of group A (whish is copying which?).
  // Note: it is allowed for requests to be "chained". E.g., group A can be
  // listed as subset of B={f1,f2,f3,f4}, with f4 removed, and group B can be
  // listed as subset of A, with f2 removed.
  // Since requests can be "chained" (in the sense explained above), we must process
  // them in the "right" order. That is, if we have a req for group A that
  // depends on group B, and there are reqs for group B that depend on other groups,
  // we must process group B first. So first, let's sort group requests
  // so that they are in the correct order

  // A list storing the order in which we process the requests.
  std::list<std::string> ordered_groups;
  while (ordered_groups.size()<m_group_requests.size()) {
    // If this iteration of the while loop does not increase the size of ordered_groups,
    // then there are circular deps, and we need to error out.
    size_t curr_size = ordered_groups.size();

    for (const auto& greqs : m_group_requests) {
      if (ekat::contains(ordered_groups,greqs.first)) {
        // We already processed this group
        continue;
      }

      // In order for this group to be able to be processed,
      // it must only have deps on groups already added to ordered_groups.
      bool can_process_this_group = true;
      for (const auto& req : greqs.second) {
        if (req.derived_type!=DerivationType::None &&
            not ekat::contains(ordered_groups,req.src_name)) {
          can_process_this_group = false;
          break;
        }
      }

      if (can_process_this_group) {
        // Ok, this group does not depend on any other group, or it depends
        // on groups we already processed. Append it to ordered_groups.
        ordered_groups.push_back(greqs.first);
      }
    }

    EKAT_REQUIRE_MSG (curr_size<ordered_groups.size(),
        "Error! There are circular dependencies between group requests.\n"
        "       The FieldManager cannot handle them.\n");

  }

  // Now we have an order in which we can process group requests. Go in order.
  for (const auto& gname : ordered_groups) {
    const auto& reqs = m_group_requests.at(gname);

    // First, we add any field needed to the corresponding FieldGroupInfo
    // Note: we only need to parse Superset, Subset, and Copy requests
    auto& info = m_field_groups.at(gname);
    for (const auto& r : reqs) {
      switch(r.derived_type) {
        case DerivationType::None: // [[fallthrough]] // c++17 keyword
        case DerivationType::Import:
          break;
        case DerivationType::Copy:  // [[fallthrough]] // c++17 keyword
        case DerivationType::Superset:
          for (const auto& fname : m_field_groups.at(r.src_name)->m_fields_names) {
            // Since we don't want to sort field names lists, adding duplciates
            // would make removing them hard, so might as well check before adding names.
            if (not ekat::contains(info->m_fields_names,fname)) {
              info->m_fields_names.push_back(fname);
            }
          }
          break;
        case DerivationType::Subset:
          for (const auto& fname : m_field_groups.at(r.src_name)->m_fields_names) {
            // Since we don't want to sort field names lists, adding duplciates
            // would make removing them hard, so might as well check before adding names.
            if (not ekat::contains(info->m_fields_names,fname) &&
                not ekat::contains(r.exclude,fname)) {
              info->m_fields_names.push_back(fname);
            }
          }
          break;
      }
    }

    // Now that all groups contain all fields, we can check consistency.
    for (const auto& req : reqs) {
      // Depending on the derivation type (if any), do different checks
      switch (req.derived_type) {
        case DerivationType::None:
          // Nothing to do for this request.
          break;
        case DerivationType::Import:
          // An imported group cannot depend on other groups. That is, we don't allow
          // to do something like G1 = import(G2,from_grid_blah) + G3.
          // We allow, however, to have multiple 'Import' requests, simply with different
          // pack sizes or bundling requests
          for (const auto& r : reqs) {
            if (r.derived_type!=DerivationType::None) {
              EKAT_REQUIRE_MSG (r.derived_type==DerivationType::Import,
                  "Error! An Impor request cannot be mixed with non-Import requests.\n"
                  "  group name: " + req.name + "\n"
                  "  second req type: " + e2str(req.derived_type) + "\n");
              EKAT_REQUIRE_MSG (r.src_name==req.src_name &&
                                r.src_grid==req.src_grid,
                  "Error! An Impor request cannot be mixed with other Import requests\n"
                  "       of different groups or from different grids.\n"
                  "  group name: " + req.name + "\n"
                  "  current req source group: " + req.src_name + "\n"
                  "  current req source grid: " + req.src_grid + "\n"
                  "  second req source group: " + r.src_name + "\n"
                  "  second req source grid: " + r.src_grid + "\n");
            }
          }

          break;
        case DerivationType::Copy:
          // We can only copy a group from this grid
          EKAT_REQUIRE_MSG (req.src_grid==m_grid->name(),
              "Error! Use 'Copy' group requests to copy a group from the same grid.\n");
          // Make sure we're not mixing 'Copy' with other requests types other than Copy or None
          for (const auto& r : reqs) {
            if (r.derived_type!=DerivationType::None) {
              EKAT_REQUIRE_MSG (r.derived_type==DerivationType::Copy,
                  "Error! A Copy request cannot be mixed with non-Copy requests.\n"
                  "  group name: " + req.name + "\n"
                  "  second req type: " + e2str(req.derived_type) + "\n");
              EKAT_REQUIRE_MSG (r.src_name==req.src_name,
                  "Error! A Copy request cannot be mixed with other Copy requests of different groups.\n"
                  "  group name: " + req.name + "\n"
                  "  current req source group: " + req.src_name + "\n"
                  "  second req source group: " + r.src_name + "\n");
            }
          }

          break;
        case DerivationType::Subset:
        {
          // Say group A is a subset of B={f1,f2,f3}, with f2 removed. If someone registered
          // field f2 as part of A, then we have an inconsistency. We also have an inconsistency
          // if someone asked for A to be a superset of C={f2,f3}, since f2 was supposed to be removed.
          const auto& fnames = m_field_groups.at(gname)->m_fields_names;
          for (const auto& excl : req.exclude) {
            EKAT_REQUIRE_MSG (not ekat::contains(fnames,excl),
                "Error! Found a request for group A to be a subset of group B,\n"
                "       but one of the excluded fields was explicitly (via register_field)\n"
                "       or implicitly (via another group request) added to group A.\n"
                "  current group (A): " + gname + "\n"
                "  superset group (B): " + req.src_name + "\n"
                "  field name: " + excl + "\n");
          }
          // Also check that all fields of this group are in the source group
          const auto& rel_fnames = m_field_groups.at(req.src_name)->m_fields_names;
          for (const auto& f : fnames) {
            EKAT_REQUIRE_MSG (ekat::contains(rel_fnames,f),
                "Error! Found a request for group A to be a subset of group B,\n"
                "       but it contains a field not in the group B.\n"
                "       Most likely reason: the was also a request for group A\n"
                "       to be a superset of group C, which was not a subset of B\n"
                "  current group (A): " + f + "\n"
                "  superset group (B): " + req.src_name + "\n"
                "  field name: " + f + "\n");
          }
          break;
        }
        case DerivationType::Superset:
        {
          const auto& fnames = m_field_groups.at(gname)->m_fields_names;
          const auto& rel_fnames = m_field_groups.at(req.src_name)->m_fields_names;
          for (const auto& f : rel_fnames) {
            EKAT_REQUIRE_MSG (ekat::contains(fnames,f),
                "Error! Found a request for group A to be a superset of group B,\n"
                "       but group B contains a field not in group A.\n"
                "       Most likely reason: the was also a request for group A\n"
                "       to be a subset of group C, which was not a superset of B\n"
                "  current group (A): " + f + "\n"
                "  subset group (B): " + req.src_name + "\n"
                "  field name: " + f + "\n");
          }
          break;
        }
      }
    }
  }
}

} // namespace scream
