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
//
// This class is meant to be subclassed (e.g., for the PG2 physics grid,
// where dynamics owns the fields on that grid, and the base class's
// straight FieldManager/file interaction is not enough). Everything that
// may need to change in a derived class is a separate virtual method:
//  - get_leaf_fields/get_fields decide WHICH fields need to be inited;
//  - get_tag_rename/get_topography_file_names decide how a field/dimension
//    is named on file, when that differs from its eamxx name.
// The rest of the class (the init_*_fields methods, driving the above hooks)
// is not virtual, and a derived class is expected to reuse it as-is.
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

  // Leaf fields (still needing initialization) of group_name on grid_name.
  // IC files only store leaf fields, so a composite field is expanded into
  // its (not yet inited) leaves, recursively. A field is a composite either
  // because it is itself a registered FieldManager group (e.g., a group's
  // own monolithic field, expanded into that group's declared members --
  // which may not be its header children, if the group is an overlapping
  // subset of a bigger one, like SHOC's "turbulence_advected_tracers" is of
  // "tracers"), or because it has header children of its own without being
  // a group (e.g., "horiz_winds", whose U/V component subfields are its
  // header children). Used for the STARTUP group. Virtual, since a derived
  // class may need to change which fields are considered leaves (e.g., to
  // exclude fields owned by dynamics).
  virtual std::vector<Field>
  get_leaf_fields (const std::shared_ptr<FieldManager>& fm,
                   const std::string& group_name,
                   const std::string& grid_name);

  // Fields (still needing initialization) of group_name on grid_name, as
  // they are stored in the group, with no leaf expansion: a field is
  // skipped only if its parent is also part of group_name, since a
  // restart/topography file stores the parent as a whole, and updating a
  // parent automatically updates all its children. Used for the RESTART
  // and TOPOGRAPHY groups. Virtual, for the same reason as get_leaf_fields.
  virtual std::vector<Field>
  get_fields (const std::shared_ptr<FieldManager>& fm,
              const std::string& group_name,
              const std::string& grid_name);

  // After initializing the (leaf) fields of the STARTUP group, some parent
  // fields (whose time stamp is NOT automatically updated when their
  // children's is) may now have all their children inited: propagate the
  // time stamp to them too.
  void fixup_parents_time_stamp (const std::shared_ptr<FieldManager>& fm,
                                 const std::string& group_name,
                                 const std::string& grid_name,
                                 const util::TimeStamp& t0);

  // If m_params has a (non-empty) 'perturbed_fields' entry, apply a random
  // (relative) perturbation to each of those (already inited) GLL fields,
  // at levels below 'perturbation_minimum_pressure'. A no-op otherwise.
  void perturb_fields (const std::shared_ptr<FieldManager>& fm);

  // The random seed to use for perturb_fields, resolved from the
  // 'generate_perturbation_random_seed'/'perturbation_random_seed' params.
  int get_perturbation_seed ();

  // A (LEV-layout, int) field, used by perturb_fields to decide which
  // levels of gll_grid to perturb: level k is perturbed iff the field's
  // value at k is nonzero. A level is included iff its reference pressure
  // (computed from gll_grid's hyam/hybm geometry data) exceeds the
  // 'perturbation_minimum_pressure' param.
  Field build_perturbation_level_mask (const std::shared_ptr<const AbstractGrid>& gll_grid);

  // If m_params has an entry '$name: $value' (a number, or an array of
  // numbers, for vector fields), assign $value to f, and stamp its time.
  // Returns true if such an entry was found (and used).
  bool set_constant_field (Field& f, const util::TimeStamp& t0);

  // If m_params has an entry '$name: "$other_name"', return $other_name.
  // Otherwise, return an empty string.
  std::string get_copy_source (const std::string& name) const;

  // Read fields from filename, and stamp their time to t0. Optionally,
  // rename some layout tags (e.g., to account for a different dimension
  // name on file).
  void read_fields (const std::string& filename,
                    std::vector<Field>& fields,
                    const std::shared_ptr<const AbstractGrid>& grid,
                    const util::TimeStamp& t0,
                    const strmap_t<std::string>& tag_rename = {});

  // By default, only the GLL grid's topography fields need a rename (the
  // topography file always uses 'ncol_d' for the GLL 'ncol' dimension, so
  // that a single topography file can serve both GLL and PG2 runs).
  // Virtual, since the pg2 model init may need to remap more dimension
  // names to the CGLL grid ones.
  virtual strmap_t<std::string>
  get_tag_rename (const std::string& group_name,
                  const std::string& grid_name) const;

  // Map from eamxx field name to the (differing) name used for that field
  // in the topography file. Virtual, since a derived class may need to
  // change which fields are expected there (e.g., the PG2 grid does not
  // load phis from the topography file, since dynamics computes it).
  virtual strmap_t<std::string> get_topography_file_names () const;

  ekat::ParameterList       m_params;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
