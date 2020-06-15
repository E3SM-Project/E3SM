#ifndef SCREAM_P3_STANDALONE_FIELD_INITIALIZER_HPP
#define SCREAM_P3_STANDALONE_FIELD_INITIALIZER_HPP

#include "share/field/field_initializer.hpp"

namespace scream {

class P3StandAloneInit : public FieldInitializer
{
public:

  virtual ~P3StandAloneInit () = default;

  // The name of the initializer
  std::string name () const { return "P3StandAloneInit"; }
  
  // Initialize fields
  void initialize_fields ();

  const std::set<FieldIdentifier>& get_inited_fields () const {
    return m_fields_id;
  }

protected:

  void add_field (const field_type& f);

  std::map<std::string,const field_type>  m_fields;

  std::set<FieldIdentifier> m_fields_id;
};

} // namespace scream

#endif // SCREAM_P3_STANDALONE_FIELD_INITIALIZER_HPP
