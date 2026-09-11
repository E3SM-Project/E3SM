#ifndef EAMXX_MODEL_INIT_HPP
#define EAMXX_MODEL_INIT_HPP

#include "share/data_managers/field_manager.hpp"

#include <ekat_parameter_list.hpp>

namespace scream
{

// Handles initialization of a FieldManager's fields, either from a startup
// (initial run) or restart run. Fields are inited, in order of precedence,
// from
//   1) a constant value (scalar or per-component) specified in the input
//      parameter list (startup run only),
//   2) a copy of another field, also specified in the input parameter list
//      (startup run only),
//   3) a startup/restart/topography file.
// The parameter list uses the same schema as the (historical) 'initial_conditions'
// input parameter list: an entry '$field_name: $value' provides a constant
// IC for field $field_name (a scalar, or, for vector fields, either a scalar
// or an array with one value per component), while an entry
// '$field_name: "$other_field_name"' means $field_name should be
// initialized as a copy of $other_field_name. The 'filename' entry gives the
// file to read fields from: for a startup run, this is the IC file; for a
// restart run, the caller (which needs access to the case's rpointer file
// to resolve the actual restart file name) is responsible for resolving the
// restart file name and passing it as 'filename' instead.
class ModelInit {
public:
  ModelInit (const ekat::ParameterList& params);

  virtual ~ModelInit () = default;

  virtual void run (const std::shared_ptr<FieldManager>& fm,
                    const util::TimeStamp& t0,
                    const RunType run_type);

protected:

  template<typename T>
  using strmap_t = std::map<std::string,T>;

  using strvec_t = std::vector<std::string>;

  // Fields in group_name (on grid_name) that still need to be initialized.
  //  - For the STARTUP group, a field that is the parent of other fields
  //    (e.g., a group's monolithic field, or a field with convenience
  //    component subfields, like U/V for horiz_winds) is expanded into its
  //    (not yet inited) children, recursively, since IC files only store
  //    leaf fields.
  //  - For every other group, a field is skipped if its parent is also
  //    part of group_name, since restart/topography files store the parent
  //    only, and updating a parent automatically updates all its children.
  // Virtual, since the pg2 model init may change the field names
  virtual std::vector<Field>
  get_fields (const std::shared_ptr<FieldManager>& fm,
              const std::string& group_name,
              const std::string& grid_name);

  // Startup (initial) run: init fields in the STARTUP group of grid, from
  // (in order of precedence) a constant value, a copy of another field, or
  // the 'filename' input file.
  void init_startup_fields (const std::shared_ptr<FieldManager>& fm,
                            const std::shared_ptr<const AbstractGrid>& grid,
                            const util::TimeStamp& t0);

  // Restart run: init fields in the RESTART group of grid, from the
  // 'filename' input file (which the caller must have set to the resolved
  // restart file name).
  void init_restart_fields (const std::shared_ptr<FieldManager>& fm,
                            const std::shared_ptr<const AbstractGrid>& grid,
                            const util::TimeStamp& t0);

  // Init fields in the TOPOGRAPHY group of grid, from the
  // 'topography_filename' input file.
  void init_topography_fields (const std::shared_ptr<FieldManager>& fm,
                               const std::shared_ptr<const AbstractGrid>& grid,
                               const util::TimeStamp& t0);

  // After initializing the (leaf) fields of the STARTUP group, some parent
  // fields (whose time stamp is NOT automatically updated when their
  // children's is) may now have all their children inited: propagate the
  // time stamp to them too.
  void fixup_parents_time_stamp (const std::shared_ptr<FieldManager>& fm,
                                 const std::string& group_name,
                                 const std::string& grid_name,
                                 const util::TimeStamp& t0);

  // If m_params has an entry '$name: $value' (a number, or an array of
  // numbers, for vector fields), assign $value to f, and stamp its time.
  // Returns true if such an entry was found (and used).
  bool set_constant_field (Field& f, const util::TimeStamp& t0);

  // If m_params has an entry '$name: "$other_name"', return $other_name.
  // Otherwise, return an empty string.
  std::string get_copy_source (const std::string& name) const;

  // Read fields from filename, optionally renaming some layout tags
  // (e.g., to account for a different dimension name on file).
  void read_fields (const std::string& filename,
                    std::vector<Field>& fields,
                    const std::shared_ptr<const AbstractGrid>& grid,
                    const strmap_t<std::string>& tag_rename = {});

  // By default, only the GLL grid's topography fields need a rename (the
  // topography file always uses 'ncol_d' for the GLL 'ncol' dimension, so
  // that a single topography file can serve both GLL and PG2 runs).
  // Virtual, since the pg2 model init may need to remap more dimension
  // names to the CGLL grid ones.
  virtual strmap_t<std::string>
  get_tag_rename (const std::string& group_name,
                  const std::string& grid_name) const;

  ekat::ParameterList       m_params;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
