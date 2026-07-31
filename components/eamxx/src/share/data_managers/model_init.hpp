#ifndef EAMXX_MODEL_INIT_HPP
#define EAMXX_MODEL_INIT_HPP

#include "share/data_managers/field_manager.hpp"

#include <ekat_parameter_list.hpp>
#include <ekat_comm.hpp>

namespace scream
{

class ModelInit {
public:
  ModelInit (const ekat::ParameterList& params);

  virtual void run (const std::shared_ptr<FieldManager>& fm,
                    const util::TimeStamp& t0,
                    const RunType run_type);

protected:

  void read_input_fields ();

  template<typename T>
  using strmap_t = std::map<std::string,T>;

  using strvec_t = std::vector<std::string>;

  // Virtual, since the pg2 model init may change the field names
  virtual std::vector<Field>
  get_fields (const std::shared_ptr<FieldManager>& fm,
              const std::string& group_name,
              const std::string& grid_name);

  void read_fields (const std::string& filename,
                    std::vector<Field>& fields,
                    const std::shared_ptr<const AbstractGrid>& grid);

  // By default, no tag renamed
  virtual strmap_t<std::string> get_tag_rename (const std::string& /* grid_name */) const { return {}; }

  // If an entry '$name = $value' is found in the 'constant_fields' entry of the param list,
  // then set field $name to $value.
  void set_constant_fields (std::vector<Field>& fields,
                            const util::TimeStamp& t0);

  ekat::ParameterList       m_params;
  strmap_t<ScalarWrapper>   m_constant_values;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
