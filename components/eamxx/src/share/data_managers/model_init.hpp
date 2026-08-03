#ifndef EAMXX_MODEL_INIT_HPP
#define EAMXX_MODEL_INIT_HPP

#include "share/data_manager/fields_manager.hpp"

#include <ekat_parameter_list.hpp>
#include <ekat_comm.hpp>

namespace scream
{

class ModelInit {
public:
  ModelInit (const ekat::ParameterList& params);

  void run (const std::shared_ptr<FieldsManager>& fm,
            const util::TimeStamp& t0,
            const RunType run_type);

protected:

  void read_input_fields ();

  template<typename T>
  using strmap_t = std::map<std::string,T>;

  // Gather fields to be inited
  std::vector<Field>
  gather_fields (const std::shared_ptr<FieldsManager>& fm,
                 const RunType run_type);
  
  void read_topography_fields (strmap_t<strmap_t<Field>>& fm);

  // If user sets fields to a constant in input param list,
  // init fields and remove them from the map(s)
  void set_constant_fields (strmap_t<strmap_t<Field>>& fields);
  strmap_t<double> parse_constant_fields ();

  ekat::ParameterList             m_params;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HPP
