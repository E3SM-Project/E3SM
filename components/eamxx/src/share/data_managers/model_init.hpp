#ifndef EAMXX_MODEL_INIT_HPP
#define EAMXX_MODEL_INIT_HPP

#include "share/data_managers/field_manager.hpp"

#include <ekat_parameter_list.hpp>
#include <ekat_comm.hpp>

namespace scream
{

// Handles initialization of a FieldManager's fields, either from a startup
// (initial run) or restart run. Fields are inited, in order of precedence,
// from
//   1) a constant value (scalar or per-component) specified in the input
//      parameter list,
//   2) a copy of another field, also specified in the input parameter list,
//   3) a startup/restart/topography file.
// The parameter list uses the same schema as the (historical) 'initial_conditions'
// input parameter list: an entry '$field_name: $value' provides a constant
// IC for field $field_name (a scalar, or, for vector fields, either a scalar
// or an array with one value per component), while an entry
// '$field_name: "$other_field_name"' means $field_name should be
// initialized as a copy of $other_field_name.
class ModelInit {
public:
  ModelInit (const ekat::ParameterList& params, const ekat::Comm& comm);

  virtual ~ModelInit () = default;

  virtual void run (const std::shared_ptr<FieldManager>& fm,
                    const util::TimeStamp& t0,
                    const RunType run_type);

protected:

  template<typename T>
  using strmap_t = std::map<std::string,T>;

  using strvec_t = std::vector<std::string>;

  // Fields in group_name (on grid_name) that still need to be initialized.
  // Virtual, since the pg2 model init may change the field names
  virtual std::vector<Field>
  get_fields (const std::shared_ptr<FieldManager>& fm,
              const std::string& group_name,
              const std::string& grid_name);

  // If m_params has an entry '$name: $value' (a number, or an array of
  // numbers, for vector fields), assign $value to f, and stamp its time.
  // Returns true if such an entry was found (and used).
  bool set_constant_field (Field& f, const util::TimeStamp& t0);

  // If m_params has an entry '$name: "$other_name"', return $other_name.
  // Otherwise, return an empty string.
  std::string get_copy_source (const std::string& name) const;

  ekat::ParameterList       m_params;
  ekat::Comm                m_comm;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
