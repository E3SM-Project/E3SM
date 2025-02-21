#ifndef SCREAM_FIELD_MANAGER_HPP
#define SCREAM_FIELD_MANAGER_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/field/field.hpp"
#include "share/field/field_group.hpp"
#include "share/field/field_request.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/eamxx_types.hpp"

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

class FieldManager {
public:

  // Public types
  using header_type      = typename Field::header_type;
  using identifier_type  = typename Field::identifier_type;
  using ci_string        = typename identifier_type::ci_string;
  using repo_type        = std::map<ci_string,std::shared_ptr<Field>>;
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
  void register_field (const FieldRequest& req);
  void register_group (const GroupRequest& req);
  void registration_ends ();
  void clean_up ();

  // Adds an externally-constructed field to the FieldManager. Allows the FM
  // to make the field available as if it had been built with the usual
  // registration procedures.
  // This can be used to allow atm procs to create some helper fields internally,
  // but still leverage the FM class for certain features (e.g., I/O).
  // NOTE: the repo must be in closed state, and the FieldManager must not already
  //       store a field with the same name.
  void add_field (const Field& f);

  // Get information about the state of the repo
  int size () const { return m_fields.size(); }
  RepoState repository_state () const { return m_repo_state; }

  // Return the grid associated to this FieldManager
  grid_ptr_type get_grid () const { return m_grid; }

  // Get the group_name->group_info map of all stored groups
  const group_info_map& get_groups_info () const { return m_field_groups; }

  // Adds $field_name to group $group_name (creating the group, if necessary).
  // NOTE: if $group_name is allocated as a bundled field, this throws.
  // NOTE: must be called after registration ends
  void add_to_group (const std::string& field_name, const std::string& group_name);

  // Query for a particular field or group of fields
  bool has_field (const std::string& name) const { return m_fields.find(name)!=m_fields.end(); }
  bool has_group (const std::string& name) const { return m_field_groups.find(name)!=m_field_groups.end(); }

  const FieldIdentifier& get_field_id (const std::string& name) const;
  Field get_field (const std::string& name) const;
  Field get_field (const identifier_type& id) const { return get_field(id.name()); }
  Field& get_field (const std::string& name);
  Field& get_field (const identifier_type& id) { return get_field(id.name()); }

  FieldGroup get_field_group (const std::string& name) const;

  repo_type::const_iterator begin () const { return m_fields.cbegin(); }
  repo_type::const_iterator end   () const { return m_fields.cend();   }

  // Set the time stamp of all fields
  // TODO: I think I want to remove this. We don't want to blanket-init
  //       the time stamp. IC reader can init the ts of IC fields, then
  //       let the different atm procs update ts when needed.
  void init_fields_time_stamp (const util::TimeStamp& t0);

protected:

  // These are allowed even if registration is ongoing
  std::shared_ptr<Field> get_field_ptr(const std::string& name) const;
  std::shared_ptr<Field> get_field_ptr(const identifier_type& id) const;

  void pre_process_group_requests ();

  // The state of the repository
  RepoState           m_repo_state;

  // The actual repo.
  repo_type           m_fields;

  // When registering subfields, we might end up registering the subfield before
  // the parent field. So at registration time, simply keep track of the subfields,
  // and create them at registration_ends() time, after all other fields.
  std::map<std::string,FieldRequest> m_subfield_requests;

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

  // If some fields are registered with incomplete FID (just name and grid),
  // we 'skip' them, hoping that some other request will contain the right specs.
  // If no complete request is given for that field, we need to error out
  std::list<std::pair<std::string,std::string>> m_incomplete_requests;
};

} // namespace scream

#endif // SCREAM_FIELD_MANAGER_HPP
