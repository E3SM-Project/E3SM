#ifndef SCREAM_FIELD_MANAGER_HPP
#define SCREAM_FIELD_MANAGER_HPP

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/field/field.hpp"
#include "share/field/field_request.hpp"
#include "share/util/map_key_iterator.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <algorithm>
#include <initializer_list>
#include <map>
#include <memory>
#include <set>

namespace scream
{

 /*
  *  A database for all the persistent fields needed in an atm time step
  *  We template a field manager over the field's (real) value type.
  *  This is enough to fully deduce the type of the stored views. All views
  *  are stored on the default device.
  *
  *  The FieldManager is associated with a specific grid. While there
  *  *may* be multiple fields with the same name (e.g., a 3d scalar at
  *  level midpoints and a 3d scalar at level interfaces), we decide to
  *  enforece a *single* copy for each field. That means that, in the
  *  example above, the 3d scalar at midpoints and interfaces would have
  *  to have different names (e.g., X_mid, X_int).
  */

template<typename RealType>
class FieldManager {
public:

  // Public types
  using RT               = typename std::remove_const<RealType>::type;
  using const_RT         = typename std::add_const<RT>::type;
  using field_type       = Field<RT>;
  using const_field_type = typename Field<RT>::const_field_type;
  using header_type      = typename field_type::header_type;
  using identifier_type  = typename field_type::identifier_type;
  using ci_string        = typename identifier_type::ci_string;
  using repo_type        = std::map<ci_string,std::shared_ptr<field_type>>;
  using group_info_type  = FieldGroupInfo;
  using group_info_map   = std::map<ci_string,std::shared_ptr<group_info_type>>;
  using grid_ptr_type    = std::shared_ptr<const AbstractGrid>;

  // Constructor(s)
  explicit FieldManager (const grid_ptr_type& grid);

  // No copies, cause the internal database is not a shared_ptr.
  // NOTE: you can change this if you find that copies are needed/useful.
  FieldManager (const FieldManager&) = delete;
  FieldManager& operator= (const FieldManager&) = delete;

  // Change the state of the database
  void registration_begins ();
  void register_field (const FieldIdentifier& fid) { register_field(FieldRequest(fid)); }
  void register_field (const FieldRequest& req);
  void register_group (const GroupRequest& req);
  void registration_ends ();
  void clean_up ();


  // Get information about the state of the repo
  int size () const { return m_fields.size(); }
  RepoState repository_state () const { return m_repo_state; }

  // Return the grid associated to this FieldManager
  grid_ptr_type get_grid () const { return m_grid; }

  // Get the group_name->group_info map of all stored groups
  const group_info_map& get_groups_info () const { return m_field_groups; }

  // Query for a particular field or group of fields
  bool has_field (const std::string& name) const { return m_fields.find(name)!=m_fields.end(); }
  bool has_field (const identifier_type& id) const;
  bool has_group (const std::string& name) const { return m_field_groups.find(name)!=m_field_groups.end(); }

  field_type get_field (const std::string& name) const;
  field_type get_field (const identifier_type& id) const;

  // Unlike the previous two, these are allowed even if registration is ongoing
  std::shared_ptr<field_type> get_field_ptr(const std::string& name) const;
  std::shared_ptr<field_type> get_field_ptr(const identifier_type& id) const;

  FieldGroup<RT> get_field_group (const std::string& name) const;
  FieldGroup<const_RT> get_const_field_group (const std::string& group_name) const;

  // Set the time stamp of all fields
  // TODO: I think I want to remove this. We don't want to blanket-init
  //       the time stamp. IC reader can init the ts of IC fields, then
  //       let the different atm procs update ts when needed.
  void init_fields_time_stamp (const util::TimeStamp& t0);

protected:

  // The state of the repository
  RepoState           m_repo_state;

  // The actual repo.
  repo_type           m_fields;

  // The map group_name -> FieldGroupInfo
  group_info_map      m_field_groups;

  // Groups need to be created after all fields have been registered,
  // since we may need to rearrange fields inside them. Also, we
  // need all groups to be registered in order to start analyzing
  // the groups request (due to groups possibly having 'parents').
  // So store GroupRequest objects during registration phase.
  std::map<std::string,std::set<GroupRequest>> m_group_requests;

  // The grid where the fields in this FM live
  std::shared_ptr<const AbstractGrid> m_grid;
};

// ============================== IMPLEMENTATION ============================= //

template<typename RealType>
FieldManager<RealType>::
FieldManager (const grid_ptr_type& grid)
  : m_repo_state (RepoState::Clean)
  , m_grid       (grid)
{
  // m_field_groups["tracers"] = std::make_shared<group_info_type>("tracers");

  EKAT_REQUIRE_MSG (m_grid!=nullptr,
      "Error! Input grid pointer is not valid.");
}

template<typename RealType>
void FieldManager<RealType>::register_field (const FieldRequest& req)
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
    m_fields[id.name()] = std::make_shared<field_type>(id);
  } else {
    using ekat::units::to_string;
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

  // Make sure the field can accommodate the requested value type
  constexpr int real_size = sizeof(RealType);
  m_fields[id.name()]->get_header().get_alloc_properties().request_allocation(real_size,req.pack_size);

  // Finally, add the field to the given groups
  // Note: we do *not* set the group info struct in the field header yet.
  //       we will do that when we end the registration phase.
  //       The reason is that at registration_ends we will know *all* the groups
  //       that each field belongs to.
  for (const auto& group_name : req.groups) {
    // Get group (and init ptr, if necessary)
    auto& group = m_field_groups[group_name];
    if (group==nullptr) {
      group = std::make_shared<group_info_type>(group_name);
    }

    // Add the field name to the list of fields belonging to this group
    if (ekat::find(group->m_fields_names,id.name())==group->m_fields_names.end()) {
      group->m_fields_names.push_back(id.name());
    }
  }
}

template<typename RealType>
void FieldManager<RealType>::register_group (const GroupRequest& req)
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

template<typename RealType>
bool FieldManager<RealType>::has_field (const identifier_type& id) const {
  return has_field(id.name()) && m_fields.at(id.name())->get_header()->get_identifier()==id;
}

template<typename RealType>
typename FieldManager<RealType>::field_type
FieldManager<RealType>::get_field (const identifier_type& id) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(id);
  EKAT_REQUIRE_MSG(ptr!=nullptr,
      "Error! Field identifier '" + id.get_id_string() + "' not found.\n");
  return *ptr;
}

template<typename RealType>
typename FieldManager<RealType>::field_type
FieldManager<RealType>::get_field (const std::string& name) const {

  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(name);
  EKAT_REQUIRE_MSG(ptr!=nullptr, "Error! Field " + name + " not found.\n");
  return *ptr;
}

template<typename RealType>
FieldGroup<typename FieldManager<RealType>::RT>
FieldManager<RealType>::
get_field_group (const std::string& group_name) const
{
  // Sanity checks
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");
  EKAT_REQUIRE_MSG (has_group(group_name),
      "Error! Field group '" + group_name + "' not found.\n");

  // Create an empty group
  FieldGroup<RT> group(group_name);

  // Set the info in the group
  group.m_info = m_field_groups.at(group_name);

  // Find all the fields and set them in the group
  for (const auto& fname : group.m_info->m_fields_names) {
    auto f = get_field_ptr(fname);
    group.m_fields[fname] = f;
  }

  // Fetch the bundle field (if bundled)
  if (group.m_info->m_bundled) {
    // We can get the parent from any of the fields in the group.
    auto p = group.m_fields.begin()->second->get_header().get_parent().lock();
    EKAT_REQUIRE_MSG(p!=nullptr,
        "Error! A field belonging to a bundled field group is missing its 'parent'.\n");

    const auto& p_id = p->get_identifier();
    group.m_bundle = get_field_ptr(p_id);
  }

  return group;
}

template<typename RealType>
FieldGroup<typename FieldManager<RealType>::const_RT>
FieldManager<RealType>::
get_const_field_group (const std::string& group_name) const
{
  // Saniti checks
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");
  EKAT_REQUIRE_MSG (has_group(group_name),
      "Error! Field group '" + group_name + "' not found.\n");

  // Create an empty group
  FieldGroup<const_RT> group(group_name);

  // Set the info in the group
  group.m_info = m_field_groups.at(group_name);

  // Find all the fields and set them in the group
  for (const auto& fname : group.m_info->m_fields_names) {
    auto f = get_field_ptr(fname);
    auto cf = std::make_shared<const_field_type>(f->get_const());
    group.m_fields[fname] = cf;
  }

  // Fetch the bundle field (if bundled)
  if (group.m_info->m_bundled) {
    // We can get the parent from any of the fields in the group.
    auto p = group.m_fields.begin()->second->get_header().get_parent().lock();
    EKAT_REQUIRE_MSG(p!=nullptr,
        "Error! A field belonging to a bundled field group is missing its 'parent'.\n");

    const auto& p_id = p->get_identifier();
    group.m_bundle = std::make_shared<const_field_type>(get_field_ptr(p_id)->get_const());
  }

  return group;
}

template<typename RealType>
void FieldManager<RealType>::
init_fields_time_stamp (const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot set initial time stamp until registration has completed.\n");

  for (auto it : m_fields) {
    it.second->get_header().get_tracking().update_time_stamp(t0);
  }
}

template<typename RealType>
void FieldManager<RealType>::registration_begins ()
{
  // Update the state of the repo
  m_repo_state = RepoState::Open;
}

template<typename RealType>
void FieldManager<RealType>::
registration_ends ()
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
  //     groups (querying m_field_groups info structs). If a GR is 'relative' to another
  //     GR, like a 'Child' group (see GR header for details), we make sure the group
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
  //     field that accommodate all bundler request, it will return an empty list.
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
  //        one at a time, then two at a time, then three at a time,... untile we
  //        reach a point where the remaining groups can all be bundled. This is overly
  //        complicated, so if 4 fails, we simply start removing groups with "Preferred"
  //        bundling until 4 succeeds or we run out of groups with bundling=preferred
  //        to remove.
  //

  // This lambda process a group request, and ensures the group contains only what it should contain
  auto ensure_group_members_correctness = [&] (const GroupRequest& r) {
    if (r.relative==nullptr || r.relative_type==Relationship::Parent) {
      // This request is either not related to another group, or it is a superset
      // of another group. Either way, the group of this request must exist.
      // TODO: should this be removed? What if a group is 'optional'?
      EKAT_REQUIRE_MSG (m_field_groups.find(r.name)!=m_field_groups.end(),
          "Error! Found a requested group that has no fields associated to it.\n"
          "    -  Group name: " + r.name + "\n"
          "    -  Grid name:  " + r.grid + "\n");
    } else {
      // If this request is a relative (alias or subset) of a relative group, then the
      // relative group must exist.
      EKAT_REQUIRE_MSG(m_field_groups.find(r.relative->name)!=m_field_groups.end(),
          "Error! Found a requested group with a relative that has no fields associated to it.\n"
          "    -  Group name:    " + r.name + "\n"
          "    -  Relative name: " + r.relative->name + "\n"
          "    -  Relative type: " + e2str(r.relative_type) + "\n"
          "    -  Grid name:     " + r.grid + "\n");
    }

    if (r.relative!=nullptr) {
      if (r.relative_type==Relationship::Child || r.relative_type==Relationship::Alias) {
        // All fields in the child/alias group must be added to this group
        auto& members = m_field_groups.at(r.name)->m_fields_names;
        const auto& relatives = m_field_groups.at(r.relative->name)->m_fields_names;
        for (const auto& n : relatives) {
          if (!ekat::contains(members,n)) {
            members.push_back(n);
          }
        }
      } else if (r.relative_type==Relationship::Parent) {
        // We take all fields in the parent group and add them to this group,
        // provided that they are not in the 'exclude' list.
        auto& members = m_field_groups.at(r.name)->m_fields_names;
        const auto& relatives = m_field_groups.at(r.relative->name)->m_fields_names;
        for (const auto& n : relatives) {
          if (!ekat::contains(members,n) && !ekat::contains(r.exclude,n)) {
            members.push_back(n);
          }
        }
      }
    }

    // Additional check. Say group G1 is a subset of G2=(f1,f2,f3), excluding f2.
    // If by chance, someone registered field f2 as part of G1, we have inconsistent
    // requests. Rather than allow sneaky bugs to go in, let's error out.
    const auto& members = m_field_groups.at(r.name)->m_fields_names;
    for (const auto& e : r.exclude) {
      EKAT_REQUIRE_MSG (not ekat::contains(members,e),
          "Error! Found group containing fields that were supposed to be excluded from relative group fields list.\n"
          "       Most likely cause is that someone individually registered a field as part of this group, while\n"
          "       someone else requested this group as subset of another group, excluding this field.\n"
          "    -  Group name: " + r.name + "\n"
          "    -  Grid name:  " + r.grid + "\n"
          "    -  Field name: " + e + "\n");

    }
  };

  // First, check groups are well defined, and contain the correct members
  for (const auto& greqs : m_group_requests) {
    for (const auto& r : greqs.second) {
      ensure_group_members_correctness (r);
    }
  }
  for (const auto& it : m_field_groups) {
    EKAT_REQUIRE_MSG (it.second->size()>0,
        "Error! We have a group in this Field Manager that has no members.\n"
        "       Group name: " + it.first + "\n"
        "       Grid name:  " + m_grid->name() + "\n");
  }

  // Gather a list of groups to be bundled
  std::list<std::string> groups_to_bundle;
  for (const auto& greqs : m_group_requests) {
    for (const auto& r : greqs.second) {
      if (r.bundling!=Bundling::NotNeeded) {
        // There's at least one request for this group to be bunlded.
        groups_to_bundle.push_back(r.name);
        break;
      }
    }
  }

  // Homme currently wants qv to be the first tracer. We should be able to
  // modify Homme, to use something like qv_idx. However, that requires
  // extensive changes in Homme. Instead, we hack our way around this limitatin
  // (for now), and rearrange groups/fields so that we can expect qv to be the
  // first tracer.
  if (has_field("qv")) {
    auto it = ekat::find(groups_to_bundle,"tracers");
    if (it!=groups_to_bundle.end()) {
      // Bring tracers to the front
      std::swap(*it,groups_to_bundle.front());

      // Adding the 'fake' group G=(qv) at the front of groups_to_bundle ensures qv won't be put
      // in the middle of the tracers group. We use a highly unlikely group name, to avoid clashing
      // with a real group name. Later, after having found an global ordering for the tracers fields,
      // we will remove this group.
      groups_to_bundle.push_front("__qv__");
      m_field_groups.emplace("__qv__",std::make_shared<group_info_type>("__qv__"));
      m_field_groups.at("__qv__")->m_fields_names.push_back("qv");
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
        groups_fields.back().sort();
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
      //  - check if there's a group in the cluster containing all the fields. If yes, use that
      //    group name for the bundled field, otherwise make one up from the names of all groups.
      //  - allocate the cluster field F.
      //  - loop over the groups in the cluster, and subview F at the proper (contiguous) indices.

      // WARNING: this lines should be removed if we move away from Homme's requirement that
      //          qv be the first tracer
      auto qv_it = ekat::find(cluster_ordered_fields,"qv");
      if (qv_it!=cluster_ordered_fields.end()) {
        // Check that qv comes first or last (if last, reverse the list). If not, error out.
        // NOTE: I *think* this should not happen, unless 'tracers' is a subgroup of a bigger group.
        EKAT_REQUIRE_MSG(qv_it==cluster_ordered_fields.begin() || qv_it==(cluster_ordered_fields.end()--),
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
          C_ap.request_allocation(sizeof(RealType),req.pack_size);
        }
      }

      // Allocate
      C->allocate_view();

      // Note: as of 02/2021, idim should *always* be 1, but we store it just in case,
      //       to avoid bugs in the future.
      const auto& C_tags = C->get_header().get_identifier().get_layout().tags();
      const int idim = std::distance(C_tags.begin(),ekat::find(C_tags,CMP));

      if (auto it = ekat::contains(cluster,"__qv__")) {
        // Erase the 'fake group' we added (to guarantee qv would be first/last in the tracers)
        m_field_groups.erase(m_field_groups.find("__qv__"));
        cluster.erase(ekat::find(cluster,"__qv__"));
      }

      // Create all individual subfields
      for (const auto& fn : cluster_ordered_fields) {
        const auto pos = ekat::find(cluster_ordered_fields,fn);
        const auto idx = std::distance(cluster_ordered_fields.begin(),pos);

        const auto& f = m_fields.at(fn);
        const auto& fh = f->get_header();
        auto fi = C->subfield(fn,fh.get_identifier().get_units(),idim,idx);

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
  //  - allocate remaining fields
  //  - update the tracking of each field, by adding the group info of all
  //    the groups the field belongs to.

  for (auto& it : m_fields) {
    if (it.second->is_allocated()) {
      // If the field has been already allocated, then it was in a bunlded group, so skip it.
      continue;
    }

    // A brand new field. Allocate it
    it.second->allocate_view();
  }

  for (const auto& it : m_field_groups) {
    auto fgi = it.second;

    // Get fields in this group
    const auto& fnames = fgi->m_fields_names;
    for (const auto& fn : fnames) {
      // Update the field tracking
      m_fields.at(fn)->get_header().get_tracking().add_to_group(fgi);
    }
  }

  // Prohibit further registration of fields
  m_repo_state = RepoState::Closed;
}

template<typename RealType>
void FieldManager<RealType>::clean_up() {
  // Clear the maps
  m_fields.clear();
  m_field_groups.clear();

  // Reset repo state
  m_repo_state = RepoState::Clean;
}

template<typename RealType>
std::shared_ptr<typename FieldManager<RealType>::field_type>
FieldManager<RealType>::get_field_ptr (const identifier_type& id) const {
  auto it = m_fields.find(id.name());
  if (it!=m_fields.end() && it->second->get_header().get_identifier()==id) {
    return it->second;
  }
  return nullptr;
}

template<typename RealType>
std::shared_ptr<typename FieldManager<RealType>::field_type>
FieldManager<RealType>::get_field_ptr (const std::string& name) const {
  auto it = m_fields.find(name);
  return it==m_fields.end() ? nullptr : it->second;
}

} // namespace scream

#endif // SCREAM_FIELD_MANAGER_HPP
