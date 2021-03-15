#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include "share/grid/grids_manager.hpp"
#include "share/field/field.hpp"
#include "share/util/map_key_iterator.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <map>
#include <set>

namespace scream
{

 /*
  *  A database for all the persistent fields needed in an atm time step
  *  We template a field repository over the field's (real) value type.
  *  This is enough to fully deduce the type of the stored views. All views
  *  are stored on the default device.
  *
  *  The fields are internally organized by name. Within each name,
  *  there can be multiple fields, which differ by their layout.
  *  For instance, we could have two version of 'temperature',
  *  one on the "physics" grid (tags: Column, Level), and one on
  *  the "dynamics" grid (tags: Element, GaussPoint, GaussPoint, Level).
  *
  *  When you query the repo for its 'size', you get the number
  *  of different *names*. So if the repo is storing the two
  *  versions of 'temperature' above and nothing else, its size
  *  will be 1. To get the number of different Field objects
  *  stored, you need to use the 'internal_size' method.
  *
  */

template<typename RealType>
class FieldRepository {
public:

  // Public types
  using RT               = typename std::remove_const<RealType>::type;
  using const_RT         = typename std::add_const<RT>::type;
  using field_type       = Field<RT>;
  using const_field_type = typename Field<RT>::const_field_type;
  using header_type      = typename field_type::header_type;
  using identifier_type  = typename field_type::identifier_type;
  using ci_string        = typename identifier_type::ci_string;
  using alias_map_type   = std::map<identifier_type,std::shared_ptr<field_type>>;
  using repo_type        = std::map<ci_string,alias_map_type>;
  using group_info_type  = FieldGroupInfo;
  using group_info_map   = std::map<ci_string,std::shared_ptr<group_info_type>>;

  // Constructor(s)
  FieldRepository ();

  // No copies, cause the internal database is not a shared_ptr.
  // NOTE: you can change this if you find that copies are needed/useful.
  FieldRepository (const FieldRepository&) = delete;
  FieldRepository& operator= (const FieldRepository&) = delete;

  // Change the state of the database
  // Note: the GridsManager is needed *if and only if* the field group TRACERS
  //       is non-empty. If so, the GM is used to retrieve a FieldLayout of
  //       a vector field on every grid on which tracers are declared.
  //       If you are writing a unit test with no use of 'TRACERS' group,
  //       then you do *not* need to pass a GM pointer.
  // Note: there is also the 'TRACERS TENDENCY' group, with some requirements
  //       on his specs. In particular, for every 'blah_tendency' field in the
  //       group 'TRACERS TENDENCY' there must be a 'blah' field in 'TRACERS'.
  void registration_begins ();
  void registration_ends (const std::shared_ptr<const GridsManager>& gm = nullptr);
  void clean_up ();

  // Deduce the pack size from the scalar type (which must be of type Pack<RealType,N>, for some int N>0, or RealType)
  template<typename RequestedValueType = RT>
  void register_field (const identifier_type& identifier, const std::initializer_list<std::string>& groups_names);

  template<typename RequestedValueType = RT>
  void register_field (const identifier_type& identifier, const std::set<std::string>& groups_names);

  template<typename RequestedValueType = RT>
  void register_field (const identifier_type& identifier, const std::string& field_group);

  template<typename RequestedValueType = RT>
  void register_field (const identifier_type& identifier);

  void register_field (const identifier_type& identifier, const int pack_size, const std::string& field_group);

  // Get information about the state of the repo
  int size () const { return m_fields.size(); }
  int internal_size () const;
  RepoState repository_state () const { return m_repo_state; }

  // Query for a particular field or group of fields
  bool has_field (const std::string& name) const;
  bool has_field (const identifier_type& identifier) const;
  field_type get_field (const identifier_type& identifier) const;
  field_type get_field(const std::string& name,const std::string& grid) const;

  // Get iterators to the keys (i.e., identifier_type) of all fields with a given name
  map_key_const_iterator<alias_map_type> cbegin (const std::string& name) const;
  map_key_const_iterator<alias_map_type> cend   (const std::string& name) const;

  // Get the names of the groups, with the names of all fields belonging to each group
  const group_info_map& get_groups_info () const { return m_field_groups; }

  // Note: when you request a group of fields, you must also specify the grid on which you need them.
  //       If you need the group G on grids foo and bar, then you must issue two separate requests to
  //       the field repo, for the group/grid pairs (G,foo) and (G,bar).
  // Note: it is *ASSUMED* that for each field in the group, there is only ONE field with such name
  //       on the requested grid.
  FieldGroup<RT> get_field_group (const std::string& name, const std::string& grid_name) const;
  FieldGroup<const_RT> get_const_field_group (const std::string& group_name, const std::string& grid_name) const;

  // Set the time stamp of all fields
  void init_fields_time_stamp (const util::TimeStamp& t0);

protected:

  std::shared_ptr<field_type> get_field_ptr(const std::string& name,const std::string& grid) const;
  std::shared_ptr<field_type> get_field_ptr(const identifier_type& id) const;

  // The state of the repository
  RepoState           m_repo_state;

  // The actual repo.
  repo_type           m_fields;

  // The map group_name -> FieldGroupInfo
  group_info_map      m_field_groups;

  // For each field, store the names of the groups it belongs to
  std::map<ci_string,std::set<ci_string>> m_field_to_groups;

  // These groups are allocated as a big bundled field, and each
  // field belonging to this group is created as a subview of that field
  // The map works m_bundled_groups[group_name] = bundled_field_name
  std::map<ci_string,ci_string>  m_bundled_groups;
};

// ============================== IMPLEMENTATION ============================= //

template<typename RealType>
FieldRepository<RealType>::
FieldRepository ()
 : m_repo_state (RepoState::Clean)
{
  m_bundled_groups["TRACERS"] = "Q";
  m_bundled_groups["TRACERS TENDENCY"] = "FQ";
  m_field_groups["TRACERS"] = std::make_shared<group_info_type>("TRACERS");
  m_field_groups["TRACERS TENDENCY"] = std::make_shared<group_info_type>("TRACERS TENDENCY");
}

template<typename RealType>
template<typename RequestedValueType>
void FieldRepository<RealType>::register_field (const identifier_type& id) {
  std::set<std::string> empty_set;
  register_field<RequestedValueType>(id,empty_set);
}

template<typename RealType>
template<typename RequestedValueType>
void FieldRepository<RealType>::
register_field (const identifier_type& id, const std::string& group_name) {
  std::set<std::string> group_name_set;
  group_name_set.insert(group_name);
  register_field<RequestedValueType>(id,group_name_set);
}

template<typename RealType>
template<typename RequestedValueType>
void FieldRepository<RealType>::
register_field (const identifier_type& id, const std::initializer_list<std::string>& groups_names) {
  register_field<RequestedValueType>(id,std::set<std::string>(groups_names));
}

template<typename RealType>
template<typename RequestedValueType>
void FieldRepository<RealType>::
register_field (const identifier_type& id, const std::set<std::string>& groups_names) {

  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Clean,
      "Error! Repo state is not 'Open' yet. You must call registration_begins() first.\n");
  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Closed,
      "Error! Repo state is not 'Open' anymore. You must register field before calling registration_ends().\n");

  using ekat::ScalarTraits;

  // Check that the requested value type is either RealType or Pack<RealType,N>, for some N>0.
  static_assert(std::is_same<RealType,RequestedValueType>::value ||
                std::is_same<RealType,typename ScalarTraits<RequestedValueType>::scalar_type>::value,
                "Error! The template argument 'RequestedValueType' of this function must either match "
                "the template argument 'RealType' of this class or be a Pack type based on RealType.\n");

  // Sanity checks
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Open,"Error! Registration of new fields not started or no longer allowed.\n");

  // Get the map of all fields with this name
  auto& map = m_fields[id.name()];

  if (map.size()>0) {
    using ekat::units::to_string;
    // Make sure the new field id stores the same units as all the other ones.
    // TODO: this is the easiest way to ensure everyone uses the same units.
    //       However, in the future, we *may* allow different units, providing
    //       the users with conversion routines perhaps.
    EKAT_REQUIRE_MSG(id.get_units()==map.begin()->first.get_units(),
                       "Error! Request to register field '" + id.name() + "' with units '" +
                       to_string(id.get_units()) + "',\n"
                       "       but there is already a request for the same field with units '" +
                       to_string(map.begin()->first.get_units()) + "'\n"
                       "       Please, check and make sure all atmosphere processes use the same units.\n");
  }

  if (map.find(id)==map.end()) {
    map[id] = std::make_shared<field_type>(id);
  }

  // Make sure the field can accommodate the requested value type
  auto& f = *map[id];
  f.get_header().get_alloc_properties().template request_allocation<RequestedValueType>();

  // Finally, add the field to the given groups
  // Note: we do *not* set the group info struct in the field header yet.
  //       we will do that when we end the registration phase.
  //       The reason is that at registration_ends we will know *all* the groups
  //       that each field belongs to.
  for (const auto& group_name : groups_names) {
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
void FieldRepository<RealType>::
register_field (const identifier_type& identifier, const int pack_size, const std::string& field_group)
{
  register_field(identifier,field_group);
  auto f = get_field_ptr(identifier);
  f->get_header().get_alloc_properties().template request_allocation<RealType>(pack_size);
}

template<typename RealType>
int FieldRepository<RealType>::
internal_size () const {
  int s = 0;
  for (auto x : m_fields) {
    s+= x.second.size();
  }
  return s;
}

template<typename RealType>
bool FieldRepository<RealType>::
has_field (const std::string& name) const {
  auto it = m_fields.find(name);
  return it!=m_fields.end();
}

template<typename RealType>
bool FieldRepository<RealType>::
has_field (const identifier_type& identifier) const {
  auto it = m_fields.find(identifier.name());
  return it!=m_fields.end() && it->second.find(identifier)!=it->second.end();
}

template<typename RealType>
typename FieldRepository<RealType>::field_type
FieldRepository<RealType>::get_field (const identifier_type& id) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(id);
  EKAT_REQUIRE_MSG(ptr!=nullptr, "Error! Field not found.\n");
  return *ptr;
}

template<typename RealType>
typename FieldRepository<RealType>::field_type
FieldRepository<RealType>::get_field (const std::string& name, const std::string& grid) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(name,grid);
  EKAT_REQUIRE_MSG(ptr!=nullptr, "Error! Field not found.\n");
  return *ptr;
}

template<typename RealType>
FieldGroup<typename FieldRepository<RealType>::RT>
FieldRepository<RealType>::
get_field_group (const std::string& group_name, const std::string& grid_name) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");

  // Create an empty group
  FieldGroup<RT> group(group_name,grid_name);

  // Allow returning an empty group
  auto it = m_field_groups.find(group_name);
  if (it!=m_field_groups.end()) {
    group.m_info = it->second;;
    for (const auto& fname : group.m_info->m_fields_names) {
      auto f = get_field_ptr(fname,grid_name);
      group.m_fields[fname] = f;
      // Fetch the bundle fields (if bundled) just once
      if (group.m_info->m_bundled && group.m_bundle==nullptr) {
        auto p = f->get_header().get_parent().lock();
        EKAT_REQUIRE_MSG(p!=nullptr,
            "Error! A field belonging to a bundled field group is missing its 'parent'.\n");

        const auto& p_id = p->get_identifier();
        group.m_bundle = get_field_ptr(p_id);
      }
    }
  }

  return group;
}

template<typename RealType>
FieldGroup<typename FieldRepository<RealType>::const_RT>
FieldRepository<RealType>::
get_const_field_group (const std::string& group_name, const std::string& grid_name) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");

  // Create an empty group
  FieldGroup<const_RT> group(group_name,grid_name);

  // Allow returning an empty group
  auto it = m_field_groups.find(group_name);
  if (it!=m_field_groups.end()) {
    group.m_info = it->second;;
    for (const auto& fname : group.m_info->m_fields_names) {
      auto f = get_field_ptr(fname,grid_name);
      auto cf = std::make_shared<const_field_type>(f->get_const());
      group.m_fields[fname] = cf;
      // Fetch the bundle fields (if bundled) just once
      if (group.m_info->m_bundled && group.m_bundle==nullptr) {
        auto p = f->get_header().get_parent().lock();
        EKAT_REQUIRE_MSG(p!=nullptr, "Error! Something is amiss with a bundled field group.\n");

        const auto& p_id = p->get_identifier();
        group.m_bundle = std::make_shared<const_field_type>(get_field_ptr(p_id)->get_const());
      }
    }
  }

  return group;
}

template<typename RealType>
void FieldRepository<RealType>::
init_fields_time_stamp (const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot set initial time stamp until registration has completed.\n");

  for (auto it : m_fields) {
    for (auto f : it.second) {
      f.second->get_header().get_tracking().update_time_stamp(t0);
    }
  }
}

template<typename RealType>
map_key_const_iterator<typename FieldRepository<RealType>::alias_map_type>
FieldRepository<RealType>::cbegin (const std::string& name) const
{
  EKAT_REQUIRE_MSG (m_fields.find(name)!=m_fields.end(),
      "Error! No field called '" << name << "' registered in the repo.\n");
  const auto& aliases = m_fields.at(name);
  return map_key_const_iterator<alias_map_type>(aliases.cbegin());
}

template<typename RealType>
map_key_const_iterator<typename FieldRepository<RealType>::alias_map_type>
FieldRepository<RealType>::cend (const std::string& name) const
{
  EKAT_REQUIRE_MSG (m_fields.find(name)!=m_fields.end(),
      "Error! No field called '" << name << "' registered in the repo.\n");
  const auto& aliases = m_fields.at(name);
  return map_key_const_iterator<alias_map_type>(aliases.cend());
}

template<typename RealType>
void FieldRepository<RealType>::registration_begins () {
  // Update the state of the repo
  m_repo_state = RepoState::Open;
}

template<typename RealType>
void FieldRepository<RealType>::
registration_ends (const std::shared_ptr<const GridsManager>& gm) {

  // Count fields in the 'TRACERS' group
  FieldGroupInfo& tr_gr = *m_field_groups["TRACERS"];
  const int nt = tr_gr.m_fields_names.size();
  std::set<std::string> tr_grids;
  for (auto& fn : tr_gr.m_fields_names) {
    // Scan all the fields with this name, and get all grids
    const auto& map = m_fields.at(fn);
    for (const auto& it_f : map) {
      tr_grids.insert(it_f.second->get_header().get_identifier().get_grid_name());
    }
  }

  // Homme (and possibly other places in EAM) hard-codes water vapor as the 1st
  // tracer in Q. Therfore, find qv in the TRACERS group, and make sure it's the 1st
  auto qv_pos = ekat::find(tr_gr.m_fields_names,"qv");
  if (qv_pos!=tr_gr.m_fields_names.end()) {
    std::iter_swap(tr_gr.m_fields_names.begin(),qv_pos);
  }

  // Inspect the TRACERS TENDENCY group, and make sure that:
  //  - has the same size as the TRACERS group
  //  - its fields are called "blah_tendency", and "blah" is a name
  //    that appears in the TRACERS group
  //  - "blah" and "blah_tendency" have the same index in the two groups.
  // At this moment, it is possible the tendency group is smaller than the
  // tracers group (if there's a tracer that no atm proc is interested in
  // its specific tendency). But it *can't* happen that the tendency group
  // is larger than the tracers group.

  auto& tr_tend_gr = *m_field_groups["TRACERS TENDENCY"];
  EKAT_REQUIRE_MSG (tr_tend_gr.size()<=tr_gr.size(),
      "Error! Tracers tendency group is larger than the tracers group.\n"
      "       We don't think this makes sense, so we are erroring out.\n"
      "       If you think this CAN happen, remove this check (and add others).\n");

  // Make sure tend is of the form "blah_tendency", and that "blah"
  // is in the tracers group
  for (const auto& tend : tr_tend_gr.m_fields_names) {
    bool ok = false;
    for (const auto& tr : tr_gr.m_fields_names) {
      if (tend==(tr+"_tendency")) {
        ok = true;
        break;
      }
    }
    EKAT_REQUIRE_MSG (ok,
        "Error! Tracer tendency '" << tend << "' has an invalid name.\n"
        "       Tendency names must be 'blah_tendency', with 'blah' being a valid tracer name.\n");
  }

  // Empty the tracers tendencies group. When we create tracers, we will re-fill it,
  // making sure the tracer tendencies are in the same order as the tracers
  tr_tend_gr = FieldGroupInfo(tr_tend_gr.m_group_name);

  // If there are tracers, we need gm to be valid
  EKAT_REQUIRE_MSG (nt==0 || gm!=nullptr,
      "Error! Tracers allocation requires a valid grid manager pointer.\n");

  const std::string& Q_name = m_bundled_groups.at("TRACERS");
  const std::string& FQ_name = m_bundled_groups.at("TRACERS TENDENCY");
  constexpr auto VAR = ShortFieldTagsNames::VAR;
  const auto q_units = ekat::units::Units::nondimensional();

  // Create "bundled" tracers and tracers forcing
  for (const auto& gn : tr_grids) {
    // Create field id for Q
    auto layout = gm->get_grid(gn)->get_3d_vector_layout(true,VAR,nt);
    FieldIdentifier fid_Q (Q_name,layout,q_units,gn);

    // Create Q field
    register_field(fid_Q);
    auto Q = get_field_ptr(fid_Q);

    // Scan all tracers on this grid, get their alloc prop,
    // and make sure Q can accommodate all of them
    auto& Q_ap = Q->get_header().get_alloc_properties();
    for (const auto& fn : tr_gr.m_fields_names) {
      auto q = get_field_ptr (fn, gn);
      if (q!=nullptr) {
        Q_ap.request_allocation(q->get_header().get_alloc_properties());
      }
    }

    // Allocate
    Q->allocate_view();

    // Register tracers forcing. Alloc props must at least accommodate Q's ones.
    FieldIdentifier fid_FQ (FQ_name,layout, q_units/ekat::units::s, gn);
    register_field(fid_FQ);
    auto FQ = get_field_ptr(fid_FQ);
    FQ->get_header().get_alloc_properties().request_allocation(Q_ap);
    FQ->allocate_view();
  }

  // Helper lambdas to detect if a field name corresponds to a tracer or
  // tracer tendency.
  const ci_string tend_tail = "_tendency";
  auto is_tracer = [&](const std::string& name) -> bool {
    return ekat::contains(tr_gr.m_fields_names,name);
  };
  auto is_tendency = [&](const std::string& name) -> bool {
    return name.size()>tend_tail.size() &&
           (name.substr(name.size()-tend_tail.size()) == tend_tail);
  };
  auto get_tr_name_from_tend_name = [&](const std::string& name) -> std::string {
    return name.substr(0,name.size()-tend_tail.size());
  };
  auto is_tr_tendency = [&](const std::string& name) -> bool {
    return is_tendency(name) && is_tracer(name.substr(0,name.size()-tend_tail.size()));
  };

  // Proceed to allocate other fields, and subview tracers
  for (auto& it : m_fields) {
    const auto& fname = it.first;
    if (fname==Q_name || fname==FQ_name) {
      // We already allocated Q and FQ, so skip it.
      continue;
    }

    for (auto& it_f : it.second) {
      auto f = it_f.second;
      if (is_tracer(fname)) {
        // A tracer must be a subview of the big Q field.
        auto pos = ekat::find(tr_gr.m_fields_names,fname);
        const int iq = std::distance(tr_gr.m_fields_names.begin(),pos);
        EKAT_REQUIRE_MSG (iq>=0 && iq<tr_gr.size(),
            "Error! Field '" << fname << "' is a tracer, but could not locate it in the tracer group.\n");

        const auto& fh = f->get_header();
        const auto  Q = get_field_ptr(Q_name,fh.get_identifier().get_grid_name());
        const auto& Q_tags = Q->get_header().get_identifier().get_layout().tags();
        // Note: as of 02/2021, idim should *always* be 1, but we store it just in case,
        //       to avoid bugs in the future.
        const int idim = std::distance(Q_tags.begin(),ekat::find(Q_tags,VAR));
        auto q = Q->subfield(fname,fh.get_identifier().get_units(),idim,iq);

        // Either this is the first tracer we set in the group (m_subview_dim still -1),
        // or idim should match what was already in the group info.
        EKAT_REQUIRE_MSG (tr_gr.m_subview_dim==-1 || tr_gr.m_subview_dim==idim,
            "Error! Something is amiss with the creation of tracers subviews.\n");
        tr_gr.m_subview_idx[fname] = iq;
        tr_gr.m_subview_dim = idim;
        tr_gr.m_bundled = true;

        // Overwrite f with q.
        *f = q;
      } else if (is_tr_tendency(fname)) {
        const auto tr_name = get_tr_name_from_tend_name(fname);
        auto pos = ekat::find(tr_gr.m_fields_names,tr_name);
        const int iq = std::distance(tr_gr.m_fields_names.begin(),pos);
        EKAT_REQUIRE_MSG (iq>=0 && iq<tr_gr.size(),
            "Error! Field '" << fname << "' is a tracer-tendency, but could not locate\n"
            "       the corresponding tracer name '" << tr_name << "' in the tracer group.\n");

        const auto& fh = f->get_header();
        const auto  FQ = get_field_ptr(FQ_name,fh.get_identifier().get_grid_name());
        const auto& FQ_tags = FQ->get_header().get_identifier().get_layout().tags();
        // Note: as of 02/2021, idim should *always* be 1, but we store it just in case,
        //       to avoid bugs in the future.
        const int idim = std::distance(FQ_tags.begin(),ekat::find(FQ_tags,VAR));
        auto fq = FQ->subfield(fname,fh.get_identifier().get_units(),idim,iq);

        // Either this is the first tracer we set in the group (m_subview_dim still -1),
        // or idim should match what was already in the group info.
        EKAT_REQUIRE_MSG (tr_tend_gr.m_subview_dim==-1 || tr_tend_gr.m_subview_dim==idim,
            "Error! Something is amiss with the creation of tracers subviews.\n");
        tr_tend_gr.m_subview_idx[fname] = iq;
        tr_tend_gr.m_subview_dim = idim;
        tr_tend_gr.m_bundled = true;

        // Overwrite f with fq.
        *f = fq;
      } else {
        // A completely independent field. Allocate it.
        f->allocate_view();
      }
    }
  }

  // Update the tracking of all fields with the group info of all
  // the groups they belong to
  for (const auto& fgi_it : m_field_groups) {
    auto fgi = fgi_it.second;
    // Get fields in this group
    const auto& fnames = fgi->m_fields_names;
    for (const auto& fn : fnames) {
      // Get all fields with this name
      auto& map = m_fields.at(fn);
      for (auto& it : map) {
        // Update the field tracking
        it.second->get_header().get_tracking().add_to_group(fgi);
      }
    }
  }

  // Prohibit further registration of fields
  m_repo_state = RepoState::Closed;
}

template<typename RealType>
void FieldRepository<RealType>::clean_up() {
  m_fields.clear();
  m_repo_state = RepoState::Clean;
}

template<typename RealType>
std::shared_ptr<typename FieldRepository<RealType>::field_type>
FieldRepository<RealType>::get_field_ptr (const identifier_type& id) const {
  if (!has_field(id)) {
    return nullptr;
  }
  return m_fields.at(id.name()).at(id);
}

template<typename RealType>
std::shared_ptr<typename FieldRepository<RealType>::field_type>
FieldRepository<RealType>::get_field_ptr (const std::string& name, const std::string& grid) const {

  if (!has_field(name)) {
    return nullptr;
  }

  // Keep track of the number of fields found for this name/grid combo
  std::vector<identifier_type> f_matches;

  //  Search subset of field repo for matching gridname.
  for (const auto& it : m_fields.at(name)) {
    if ( it.first.get_grid_name()==grid) {
      f_matches.push_back(it.first);
    }
  }

  // Check to make sure a) the field was found on this grid, and b) only once.
  EKAT_REQUIRE_MSG(f_matches.size()==1, "Error! get_field: " + name + " found " + std::to_string(f_matches.size()) + " matches on grid " + grid + ".\n");

  // Use this field id to grab field itself.
  return get_field_ptr(f_matches[0]);
}

} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
