#ifndef SCREAM_FIELD_MANAGER_HPP
#define SCREAM_FIELD_MANAGER_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/grid/library_grids_manager.hpp"
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
  using header_type         = typename Field::header_type;
  using identifier_type     = typename Field::identifier_type;
  using ci_string           = typename identifier_type::ci_string;
  using repo_type           = std::map<ci_string,std::map<ci_string,std::shared_ptr<Field>>>;
  using field_group_type    = std::map<ci_string,std::map<ci_string,std::shared_ptr<FieldGroup>>>;
  using group_info_map      = std::map<ci_string, std::shared_ptr<FieldGroupInfo>>;

  // Constructor(s)
  explicit FieldManager (const std::shared_ptr<const AbstractGrid>& grid);
  explicit FieldManager (const std::shared_ptr<const GridsManager>& grid);

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
  void clean_up (const std::string& grid_name);

  // Adds an externally-constructed field to the FieldManager. Allows the FM
  // to make the field available as if it had been built with the usual
  // registration procedures.
  // This can be used to allow atm procs to create some helper fields internally,
  // but still leverage the FM class for certain features (e.g., I/O).
  // NOTE: the repo must be in closed state, and the FieldManager must not already
  //       store a field with the same name.
  void add_field (const Field& f);

  // Adds $field_name on $grid_name to group $group_name (creating the group, if necessary).
  // NOTE: if $group_name is allocated as a monolithic field, this throws.
  // NOTE: must be called after registration ends
  void add_to_group (const std::string& field_name, const std::string& grid_name, const std::string& group_name);
  void add_to_group (const identifier_type& id, const std::string& group_name) { add_to_group(id.name(), id.get_grid_name(), group_name); }

  // Query for a particular field or group of fields
  bool has_field (const std::string& field_name, const std::string& grid_name) const;
  bool has_field (const std::string& field_name) const {
    EKAT_ASSERT_MSG(m_grids_mgr->size() == 1,
      "Error! More than one grid exists for FieldManager, must specify grid name to query for field.\n"
      "  - Field name: " + field_name + "\n"
      "  - Grids in FM: " + m_grids_mgr->print_available_grids() + "\n");
    return has_field(field_name, m_grids_mgr->get_repo().begin()->second->name());
  }
  bool has_group (const std::string& group_name, const std::string& grid_name) const;
  bool has_group (const std::string& group_name) const {
    EKAT_ASSERT_MSG(m_grids_mgr->size() == 1,
      "Error! More than one grid exists for FieldManager, must specify grid name to query for group.\n"
      "  - Group name: " + group_name + "\n"
      "  - Grids in FM: " + m_grids_mgr->print_available_grids() + "\n");
    return has_group(group_name, m_grids_mgr->get_repo().begin()->second->name());
  }

  Field get_field (const std::string& name, const std::string& grid_name) const;
  Field get_field (const identifier_type& id) const { return get_field(id.name(), id.get_grid_name()); }
  Field get_field (const std::string& name) const {
    EKAT_ASSERT_MSG(m_grids_mgr->size() == 1,
      "Error! More than one grid exists for FieldManager, must specify grid name to get field.\n"
      "  - Field name: " + name + "\n"
      "  - Grids in FM: " + m_grids_mgr->print_available_grids() + "\n");
    return get_field(name, m_grids_mgr->get_repo().begin()->second->name());
  }

  Field& get_field (const std::string& name, const std::string& grid_name);
  Field& get_field (const identifier_type& id) { return get_field(id.name(), id.get_grid_name()); }
  Field& get_field (const std::string& name) {
    EKAT_ASSERT_MSG(m_grids_mgr->size() == 1,
      "Error! More than one grid exists for FieldManager, must specify grid name to get field.\n"
      "  - Field name: " + name + "\n"
      "  - Grids in FM: " + m_grids_mgr->print_available_grids() + "\n");
    return get_field(name, m_grids_mgr->get_repo().begin()->second->name());
  }

  FieldGroupInfo get_group_info (const std::string& group_name) const {
    EKAT_ASSERT_MSG(m_grids_mgr->size() == 1,
      "Error! More than one grid exists for FieldManager, must specify grid name to query for FieldGroupInfo.\n"
      "  - Group name: " + group_name + "\n"
      "  - Grids in FM: " + m_grids_mgr->print_available_grids() + "\n");
    return get_group_info(group_name, m_grids_mgr->get_repo().begin()->second->name());
  }
  FieldGroupInfo get_group_info (const std::string& group_name, const std::string& grid_name) const;

  FieldGroup get_field_group (const std::string& name, const std::string& grid_name);

  const std::map<ci_string,std::shared_ptr<Field>>&
  get_repo () const {
    EKAT_ASSERT_MSG(m_grids_mgr->size() == 1,
      "Error! More than one grid exists for FieldManager, must specify grid name to get field repo.\n"
      "  - Grids in FM: " + m_grids_mgr->print_available_grids() + "\n");
    return m_fields.at(m_grids_mgr->get_repo().begin()->second->name());
  }
  const std::map<ci_string,std::shared_ptr<Field>>&
  get_repo (const std::string& grid_name) const {
    EKAT_REQUIRE_MSG(m_grids_mgr->has_grid(grid_name),
      "Error! This field manager does not contain data on grid \""+grid_name+"\"\n");
    return m_fields.at(grid_name);
  }

  const std::shared_ptr<const AbstractGrid> get_grid () const {
    EKAT_ASSERT_MSG(m_grids_mgr->size() == 1,
      "Error! More than one grid exists for FieldManager, must access grid through grids manager.\n"
      "  - Grids in FM: " + m_grids_mgr->print_available_grids() + "\n");
    return m_grids_mgr->get_repo().begin()->second;
  }

  const std::shared_ptr<const GridsManager>&
  get_grids_manager () const { return m_grids_mgr; }

  // Set the time stamp of all fields on given grid
  // TODO: I think I want to remove this. We don't want to blanket-init
  //       the time stamp. IC reader can init the ts of IC fields, then
  //       let the different atm procs update ts when needed.
  void init_fields_time_stamp (const util::TimeStamp& t0);

protected:

  void pre_process_monolithic_group_requests ();

  // The state of the repository
  RepoState m_repo_state;

  // The actual repo.
  repo_type           m_fields;

  // Preprocessed field groups
  field_group_type    m_field_groups;

  // The map group_name -> FieldGroupInfo
  group_info_map      m_field_group_info;

  // Groups need to be created after all fields have been registered,
  // since we may need to rearrange fields inside them. Also, we
  // need all groups to be registered in order to start analyzing
  // the groups request (due to groups possibly having 'parents').
  // So store GroupRequest objects during registration phase.
  std::map<std::string, std::map<std::string, std::set<GroupRequest>>> m_group_requests;

  // The grids where the fields in this FM live
  std::shared_ptr<const GridsManager> m_grids_mgr;

  // If some fields are registered with incomplete FID (just name and grid),
  // we 'skip' them, hoping that some other request will contain the right specs.
  // If no complete request is given for that field, we need to error out
  std::list<std::pair<std::string,std::string>> m_incomplete_requests;
};

} // namespace scream

#endif // SCREAM_FIELD_MANAGER_HPP
