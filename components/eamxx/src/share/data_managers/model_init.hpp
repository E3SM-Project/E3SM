#ifndef EAMXX_MODEL_INIT_HPP
#define EAMXX_MODEL_INIT_HPP

#include "share/data_managers/field_manager.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_parameter_list.hpp>

#include <any>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace scream
{

class IOPDataManager;

// Handles the initialization of the fields in a FieldManager, either for an
// initial (startup) run or for a restarted run.
//
// The input parameter list is (a copy of) the 'initial_conditions' sublist of
// the atm input parameters. The entries relevant to this class are
//  - filename: the file to read fields from. For a startup run, this is the
//    IC file. For a restart run, the caller (who has access to the case's
//    rpointer file) must resolve the actual restart file name, and set it as
//    'filename'.
//  - constant_fields: a sublist with two sublists, 'overrides' and 'defaults'.
//    Each entry is '$field_name: $value', where $value is either a scalar or,
//    for vector fields, an array with one value per component. Fields in
//    'overrides' are ALWAYS set to the given constant, even if they are in
//    the IC file. Fields in 'defaults' are set to the given constant ONLY IF
//    they are NOT in the IC file. If a field is in both lists, 'overrides' wins.
//    Only used for startup runs.
//  - topography_filename: the file to read the topography fields (phis, sgh,
//    sgh30) from. These fields are NOT read from the IC file, nor are they
//    affected by 'defaults' constants (but 'overrides' constants still apply).
//    Only used for startup runs.
//  - perturbed_fields, perturbation_limit, perturbation_minimum_pressure,
//    perturbation_random_seed, generate_perturbation_random_seed: if
//    perturbed_fields is not empty, a random relative perturbation is applied
//    to those (already inited) fields, which must live on the physics_gll grid.
//    Only used for startup runs.
//
// If an IOPDataManager is provided, fields are read from file so that every
// column gets the data of the file column closest to the IOP target lat/lon,
// and the IOP data (read at the run's t0) is used to set some of the fields on
// the physics_gll grid (if any), before any perturbation is applied.
//
// For a startup run, ModelInit inits the fields in the STARTUP group (of each
// grid); for a restart run, it inits those in the RESTART group. As in the
// atmosphere driver, we assume that these groups only contain *leaf* fields
// (e.g., U/V rather than horiz_winds), since that is what IC/restart files
// store. After a field is inited, its time stamp is set to the given t0. If all
// the children of a field's parent have been inited, the parent's time stamp
// is set as well.
//
// This class is meant to be subclassed, for components that need to handle
// part of the initialization on their own (e.g., because some of the model's
// state lives in data structures that are not (or not only) EAMxx fields).
//  - run() is virtual, so that a derived class can wrap the base class
//    behavior with its own pre/post steps;
//  - get_fields() is virtual, so that a derived class can change which fields
//    the base class considers for a given group/grid;
//  - all the other methods are protected helpers that a derived class can
//    use to init arbitrary fields/field managers.
class ModelInit {
public:
  template<typename T>
  using strmap_t = std::map<std::string,T>;

  using restart_extra_data_t = strmap_t<std::shared_ptr<std::any>>;

  ModelInit (const ekat::ParameterList& params);

  virtual ~ModelInit () = default;

  // Data (other than fields) that must be loaded from the restart file, in the
  // form of GLOBAL attributes of the same name. Each any must hold a
  // value of type int, std::int64_t, float, double, or std::string.
  // The data is only loaded for restart runs.
  void set_restart_extra_data (const restart_extra_data_t& data) {
    m_restart_extra_data = data;
  }

  // Provide an IOP data manager, if running an IOP case. Only used for startup runs.
  void set_iop_data_manager (const std::shared_ptr<IOPDataManager>& iop) {
    m_iop_data_manager = iop;
  }

  // Init all the fields in the STARTUP (RunType::Initial) or RESTART
  // (RunType::Restart) group, on all grids of the field manager.
  virtual void run (const std::shared_ptr<FieldManager>& fm,
                    const util::TimeStamp& t0,
                    const RunType run_type);

protected:

  // Init the fields of the STARTUP (RunType::Initial) or RESTART (RunType::Restart) group on
  // all grids of the field manager. For initial runs, also check that no constant is provided
  // for fields with subfields. This is the first step of run().
  void init_all_grids (const std::shared_ptr<FieldManager>& fm,
                       const util::TimeStamp& t0,
                       const RunType run_type);

  // Do what is left to do after the fields of a startup run have been inited: use the IOP
  // data (if any) to set some of the physics_gll fields, and apply the IC perturbation
  // (if requested). This is the last step of run() for startup runs.
  void finalize_startup (const std::shared_ptr<FieldManager>& fm,
                         const util::TimeStamp& t0);

  // Init all the fields in the STARTUP group of grid_name.
  void init_startup_fields (const std::shared_ptr<FieldManager>& fm,
                            const std::string& grid_name,
                            const util::TimeStamp& t0);

  // Init all the fields in the RESTART group of grid_name from the restart file.
  void init_restart_fields (const std::shared_ptr<FieldManager>& fm,
                            const std::string& grid_name,
                            const util::TimeStamp& t0);

  // Load the restart extra data from the restart file
  void load_restart_extra_data () const;

  // The fields (still needing initialization) in group_name on grid_name.
  virtual std::vector<Field>
  get_fields (const std::shared_ptr<FieldManager>& fm,
              const std::string& group_name,
              const std::string& grid_name);

  // Set the field to the constant provided in constants_pl for it.
  // The entry can be '$name: $val', or '$name: [$val]', or '$name: [$val1,...,$valN]',
  // where the values are doubles (write 1.0, not 1, since the ekat yaml parser reads 1
  // as an int). In the first two cases all entries of the field are set to val; in the
  // last case, the field must be a vector field with N components, and each component
  // is set to the corresponding entry. Integer fields can be set from (real) values, but
  // the conversion must be exact (e.g., 2.0 is ok, 2.5 throws).
  // Does NOT update the time stamp.
  void set_constant_field (Field& f, const ekat::ParameterList& constants_pl) const;

  // Read the fields from the file (on the given grid). Optionally, rename some layout
  // tags (e.g., to account for a different dimension name on file). Does NOT update
  // the time stamps. If an IOP data manager is present (and use_iop=true), use it to
  // read the fields (in which case tag_rename is ignored).
  void read_fields (const std::string& filename,
                    const std::vector<Field>& fields,
                    const std::shared_ptr<const AbstractGrid>& grid,
                    const strmap_t<std::string>& tag_rename = {},
                    const bool use_iop = true);

  // Map from eamxx field name to the name of that field in the topography file,
  // for the topography fields on grid_name. These are the STARTUP fields that
  // are read from the topography file, rather than from the IC file.
  // Virtual, since a derived class may need to change which fields are expected
  // there (e.g., because dynamics computes some of them).
  virtual strmap_t<std::string>
  get_topography_file_names (const std::string& grid_name) const;

  // Layout tags to rename when reading the topography fields of grid_name.
  // By default, only the GLL grid needs a rename (the topography file always
  // uses 'ncol_d' for the GLL 'ncol' dimension, so that a single topography
  // file can serve both GLL and PG2 runs).
  virtual strmap_t<std::string>
  get_topography_tag_rename (const std::string& grid_name) const;

  // If m_params has a (non-empty) 'perturbed_fields' entry, apply a random
  // (relative) perturbation to each of those (already inited) physics_gll
  // fields, at levels below 'perturbation_minimum_pressure'. A no-op otherwise.
  void perturb_fields (const std::shared_ptr<FieldManager>& fm);

  // The random seed to use for perturb_fields, resolved from the
  // 'generate_perturbation_random_seed'/'perturbation_random_seed' params.
  int get_perturbation_seed () const;

  // A (LEV-layout, int) field, used by perturb_fields to decide which
  // levels of gll_grid to perturb: level k is perturbed iff the field's
  // value at k is nonzero. A level is included iff its reference pressure
  // (computed from gll_grid's hyam/hybm geometry data) exceeds the
  // 'perturbation_minimum_pressure' param.
  Field build_perturbation_level_mask (const std::shared_ptr<const AbstractGrid>& gll_grid) const;

  // Set the time stamp of the fields to t0. Also set that of their children, and that
  // of their parents, if all the children of the parent have been inited.
  void update_time_stamps (const std::vector<Field>& fields,
                           const util::TimeStamp& t0) const;

  // A field with subfields (e.g., horiz_winds, with U/V as children) is not
  // part of the STARTUP group (only the subfields are), so a constant provided
  // for it would be silently ignored. Throw if that is the case.
  void check_no_constants_for_composite_fields (const std::shared_ptr<FieldManager>& fm);

  ekat::ParameterList   m_params;

  restart_extra_data_t  m_restart_extra_data;

  std::shared_ptr<IOPDataManager> m_iop_data_manager;

  // Grids for which the IOP data manager has already set up its io info
  std::set<std::string> m_iop_io_grids;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
