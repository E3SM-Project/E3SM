#ifndef EAMXX_MODEL_INIT_HPP
#define EAMXX_MODEL_INIT_HPP

#include "share/data_managers/field_manager.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_parameter_list.hpp>

#include <any>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace scream
{

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
//    Each entry is '$field_name: [$value1,...,$valueN]', an array of reals (the XML
//    defaults enforce that, so that is what atmchange and buildnml produce). A single
//    value sets all the entries of the field, while N>1 values are only valid for a
//    vector field with N components, and set each component to the corresponding value.
//    For convenience, the input.yaml of standalone tests can use a plain scalar
//    ('$field_name: $value') in place of a one-element array.
//    Fields in 'overrides' are ALWAYS set to the given constant, even if they are in
//    the IC file. Fields in 'defaults' are set to the given constant ONLY IF
//    they are NOT in the IC file. If a field is in both lists, 'overrides' wins.
//    Integer fields can be set from (real) values, as long as the conversion is exact.
//    Only used for startup runs.
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

  // Init all the fields in the STARTUP (RunType::Initial) or RESTART
  // (RunType::Restart) group, on all grids of the field manager.
  virtual void run (const std::shared_ptr<FieldManager>& fm,
                    const util::TimeStamp& t0,
                    const RunType run_type);

protected:

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
              const std::string& grid_name) const;

  // Get the constant(s) in constants_pl (the 'overrides' or 'defaults' sublist of 'constant_fields')
  // for the field with the given name. The entry can be either an array of doubles, or a single
  // double, which is the same as an array of one element (the ekat yaml parser reads '1' as an
  // int, so write 1.0). Throws if the entry has a different type.
  static std::vector<double> get_constant_values (const ekat::ParameterList& constants_pl,
                                                  const std::string& field_name);

  // Set the field to the given constants: a single value sets all the entries of the field,
  // while N>1 values require a vector field with N components, and set each component to the
  // corresponding value. For integer fields, the values must be integral (e.g., 2.0 is ok,
  // 2.5 throws), since we do not want to silently truncate.
  // Does NOT update the time stamp.
  static void set_constant_field (const Field& f, const std::vector<double>& values);

  // Read the fields from the file (on the given grid). Does NOT update the time stamps.
  void read_fields (const std::string& filename,
                    const std::vector<Field>& fields,
                    const std::shared_ptr<const AbstractGrid>& grid) const;

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
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
