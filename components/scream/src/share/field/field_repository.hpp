#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include "share/grid/grids_manager.hpp"
#include "share/field/field.hpp"
#include "share/util/map_key_iterator.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <initializer_list>
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
  *  The FieldRepository is associated with a specific grid. While there
  *  *may* be multiple fields with the same name (e.g., a 3d scalar at
  *  level midpoints and a 3d scalar at level interfaces), we decide to
  *  enforece a *single* copy for each field. That means that, in the
  *  example above, the 3d scalar at midpoints and interfaces would have
  *  to have different names (e.g., X_mid, X_int).
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
  using repo_type        = std::map<ci_string,std::shared_ptr<field_type>>;
  using group_info_type  = FieldGroupInfo;
  using group_info_map   = std::map<ci_string,std::shared_ptr<group_info_type>>;
  using grid_ptr_type    = std::shared_ptr<const AbstractGrid>;

  // Constructor(s)
  explicit FieldRepository (const grid_ptr_type& grid);

  // No copies, cause the internal database is not a shared_ptr.
  // NOTE: you can change this if you find that copies are needed/useful.
  FieldRepository (const FieldRepository&) = delete;
  FieldRepository& operator= (const FieldRepository&) = delete;

  // Change the state of the database
  void registration_begins ();
  void registration_ends ();
  void clean_up ();

  // --- Registration methods --- //
  template<typename RequestedValueType = RT>
  void register_field (const identifier_type& id, const std::string& group_name) {
    register_field<RequestedValueType>(id,{group_name});
  }

  template<typename RequestedValueType = RT>
  void register_field (const identifier_type& id,
                       const std::initializer_list<std::string>& groups_names = {}) {
    register_field<RequestedValueType>(id,std::set<std::string>(groups_names));
  }

  template<typename RequestedValueType = RT>
  void register_field (const identifier_type& id,
                       const std::set<std::string>& groups_names);

  void register_field (const identifier_type& id, const int pack_size,
                       const std::string& group_name) {
    register_field(id,pack_size,{group_name});
  }

  void register_field (const identifier_type& id, const int pack_size,
                       const std::initializer_list<std::string>& groups_names);

  // Get information about the state of the repo
  int size () const { return m_fields.size(); }
  RepoState repository_state () const { return m_repo_state; }

  // Return the grid associated to this FieldRepository
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

  // For each field, store the names of the groups it belongs to
  std::map<ci_string,std::set<ci_string>> m_field_to_groups;

  // These groups are allocated as a big bundled field, and each
  // field belonging to this group is created as a subview of that field
  // The map works m_bundled_groups[group_name] = bundled_field_name
  std::map<ci_string,ci_string>  m_bundled_groups;

  // The grid where the fields in this repo live
  std::shared_ptr<const AbstractGrid> m_grid;
};

// ============================== IMPLEMENTATION ============================= //

template<typename RealType>
FieldRepository<RealType>::
FieldRepository (const grid_ptr_type& grid)
  : m_repo_state (RepoState::Clean)
  , m_grid       (grid)
{
  m_bundled_groups["TRACERS"] = "Q";
  m_field_groups["TRACERS"] = std::make_shared<group_info_type>("TRACERS");

  EKAT_REQUIRE_MSG (m_grid!=nullptr,
      "Error! Input grid pointer is not valid.");
}

template<typename RealType>
template<typename RequestedValueType>
void FieldRepository<RealType>::
register_field (const identifier_type& id, const std::set<std::string>& groups_names) {
  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Clean,
      "Error! Repo state is not 'Open' yet. You must call registration_begins() first.\n");
  EKAT_REQUIRE_MSG (m_repo_state!=RepoState::Closed,
      "Error! Repo state is not 'Open' anymore. You already called registration_ends().\n");

  using ekat::ScalarTraits;

  // Check that the requested value type is either RealType or Pack<RealType,N>, for some N>0.
  static_assert(std::is_same<RealType,RequestedValueType>::value ||
                std::is_same<RealType,typename ScalarTraits<RequestedValueType>::scalar_type>::value,
                "Error! The template argument 'RequestedValueType' of this function must either match "
                "the template argument 'RealType' of this class or be a Pack type based on RealType.\n");

  // Sanity checks
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Open,"Error! Registration of new fields not started or no longer allowed.\n");

  // Make sure the grid name from the id matches the name of m_grid
  EKAT_REQUIRE_MSG(id.get_grid_name()==m_grid->name(),
      "Error! Input identifier stores a different grid name than the repo:\n"
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
  m_fields[id.name()]->get_header().get_alloc_properties().template request_allocation<RequestedValueType>();

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
register_field (const identifier_type& identifier, const int pack_size,
                const std::initializer_list<std::string>& groups_names) {
  register_field(identifier,groups_names);
  auto f = get_field_ptr(identifier);
  f->get_header().get_alloc_properties().template request_allocation<RealType>(pack_size);
}

template<typename RealType>
bool FieldRepository<RealType>::has_field (const identifier_type& id) const {
  return has_field(id.name()) && m_fields.at(id.name())->get_header()->get_identifier()==id;
}

template<typename RealType>
typename FieldRepository<RealType>::field_type
FieldRepository<RealType>::get_field (const identifier_type& id) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(id);
  EKAT_REQUIRE_MSG(ptr!=nullptr,
      "Error! Field identifier '" + id.get_id_string() + "' not found.\n");
  return *ptr;
}

template<typename RealType>
typename FieldRepository<RealType>::field_type
FieldRepository<RealType>::get_field (const std::string& name) const {

  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(name);
  EKAT_REQUIRE_MSG(ptr!=nullptr, "Error! Field " + name + " not found.\n");
  return *ptr;
}

template<typename RealType>
FieldGroup<typename FieldRepository<RealType>::RT>
FieldRepository<RealType>::
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
FieldGroup<typename FieldRepository<RealType>::const_RT>
FieldRepository<RealType>::
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
void FieldRepository<RealType>::
init_fields_time_stamp (const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot set initial time stamp until registration has completed.\n");

  for (auto it : m_fields) {
    it.second->get_header().get_tracking().update_time_stamp(t0);
  }
}

template<typename RealType>
void FieldRepository<RealType>::registration_begins ()
{
  // Update the state of the repo
  m_repo_state = RepoState::Open;
}

template<typename RealType>
void FieldRepository<RealType>::
registration_ends ()
{
  // Count fields in the 'TRACERS' group
  FieldGroupInfo& tracers = *m_field_groups["TRACERS"];
  const int nt = tracers.m_fields_names.size();

  // Homme (and possibly other places in EAM) hard-codes water vapor as the 1st
  // tracer in Q. Therfore, find qv in the TRACERS group, and make sure it's the 1st
  auto qv_pos = ekat::find(tracers.m_fields_names,"qv");
  if (qv_pos!=tracers.m_fields_names.end()) {
    std::iter_swap(tracers.m_fields_names.begin(),qv_pos);
  }

  const std::string& Q_name = m_bundled_groups.at("TRACERS");
  constexpr auto CMP = ShortFieldTagsNames::CMP;
  const auto q_units = ekat::units::Units::nondimensional();

  // Create "bundled" tracers
  if (nt>0) {
    // If there are tracers, we need the grid ptr to be valid
    EKAT_REQUIRE_MSG (nt==0 || m_grid!=nullptr,
        "Error! Tracers allocation requires a valid grid pointer.\n");

    // Create field id for Q
    auto layout = m_grid->get_3d_vector_layout(true,CMP,nt);
    FieldIdentifier fid_Q (Q_name,layout,q_units,m_grid->name());

    // Create Q field
    register_field(fid_Q);
    auto Q = get_field_ptr(fid_Q);

    // Scan all tracers on this grid, get their alloc prop,
    // and make sure Q can accommodate all of them
    auto& Q_ap = Q->get_header().get_alloc_properties();
    for (const auto& fn : tracers.m_fields_names) {
      auto q = get_field_ptr (fn);
      if (q!=nullptr) {
        Q_ap.request_allocation(q->get_header().get_alloc_properties());
      }
    }

    // Allocate
    Q->allocate_view();
  }

  // Helper lambda to detect if a field name corresponds to a tracer
  auto is_tracer = [&](const std::string& name) -> bool {
    return ekat::contains(tracers.m_fields_names,name);
  };

  // Proceed to allocate other fields, and subview tracers
  for (auto& it : m_fields) {
    const auto& fname = it.first;
    if (fname==Q_name) {
      // We already allocated Q , so skip it.
      continue;
    }

    auto f = it.second;
    if (is_tracer(fname)) {
      // A tracer must be a subview of the big Q field.
      auto pos = ekat::find(tracers.m_fields_names,fname);
      const int iq = std::distance(tracers.m_fields_names.begin(),pos);
      EKAT_REQUIRE_MSG (iq>=0 && iq<tracers.size(),
          "Error! Field '" << fname << "' is a tracer, but could not locate it in the tracer group.\n");

      const auto& fh = f->get_header();
      const auto  Q = get_field_ptr(Q_name);
      const auto& Q_tags = Q->get_header().get_identifier().get_layout().tags();
      // Note: as of 02/2021, idim should *always* be 1, but we store it just in case,
      //       to avoid bugs in the future.
      const int idim = std::distance(Q_tags.begin(),ekat::find(Q_tags,CMP));
      auto qi = Q->subfield(fname,fh.get_identifier().get_units(),idim,iq);

      // Either this is the first tracer we set in the group (m_subview_dim still -1),
      // or idim should match what was already in the group info.
      EKAT_REQUIRE_MSG (tracers.m_subview_dim==-1 || tracers.m_subview_dim==idim,
          "Error! Something is amiss with the creation of tracers subviews.\n");
      tracers.m_subview_idx[fname] = iq;
      tracers.m_subview_dim = idim;
      tracers.m_bundled = true;

      // Overwrite f with qi.
      *f = qi;
    } else {
      // A completely independent field. Allocate it.
      f->allocate_view();
    }
  }

  // Update the tracking of each field, by adding the group info of all
  // the groups the field belongs to.
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
void FieldRepository<RealType>::clean_up() {
  // Clear the maps
  m_fields.clear();
  m_field_groups.clear();
  m_field_to_groups.clear();
  m_bundled_groups.clear();

  // Reset repo state
  m_repo_state = RepoState::Clean;
}

template<typename RealType>
std::shared_ptr<typename FieldRepository<RealType>::field_type>
FieldRepository<RealType>::get_field_ptr (const identifier_type& id) const {
  auto it = m_fields.find(id.name());
  if (it!=m_fields.end() && it->second->get_header().get_identifier()==id) {
    return it->second;
  }
  return nullptr;
}

template<typename RealType>
std::shared_ptr<typename FieldRepository<RealType>::field_type>
FieldRepository<RealType>::get_field_ptr (const std::string& name) const {
  auto it = m_fields.find(name);
  return it==m_fields.end() ? nullptr : it->second;
}

} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
