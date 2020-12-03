#ifndef SCREAM_ZM_INPUTS_INITIALIZER_HPP
#define SCREAM_ZM_INPUTS_INITIALIZER_HPP

#include "share/field/field_initializer.hpp"
#include <unordered_map>
#include <string>

namespace scream {

class ZMInputsInitializer : public FieldInitializer
{
public:

  virtual ~ZMInputsInitializer () = default;

  // The name of the initializer
  std::string name () const override { return "ZMInputsInitializer"; }
  
  // Initialize fields
  void initialize_fields () override;

  const std::set<FieldIdentifier>& get_inited_fields () const override {
    return m_fields_id;
  }

protected:

  void add_field (const field_type& f) override;
  void add_field (const field_type& f, const field_type& f_ref,
                  const remapper_ptr_type& remapper) override;

  std::map<std::string,const field_type>  m_fields;

  std::set<FieldIdentifier> m_fields_id;

  std::shared_ptr<AbstractRemapper<Real>> m_remapper;
};

} // namespace scream

#endif // SCREAM_ZM_INPUTS_INITIALIZER_HPP
