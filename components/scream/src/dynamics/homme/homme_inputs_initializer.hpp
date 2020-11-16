#ifndef SCREAM_HOMME_INPUTS_INITIALIZER_HPP
#define SCREAM_HOMME_INPUTS_INITIALIZER_HPP

#include "share/field/field_initializer.hpp"

namespace scream {

class HommeInputsInitializer : public FieldInitializer
{
public:
  using field_type       = Field<      Real>;

  // Constructor(s) & Destructor
  virtual ~HommeInputsInitializer () = default;

  std::string name () const override { return "HommeInputsInitializer"; }

  void initialize_fields () override;
  const std::set<FieldIdentifier>& get_inited_fields () const override {
    return m_fids;
  }

protected:

  void add_field (const field_type& f) override;

  // Members
  std::set<FieldIdentifier>               m_fids;
  std::set<ekat::CaseInsensitiveString>   m_names;

  bool m_fields_inited = false;
};

} // namespace scream

#endif // SCREAM_HOMME_INPUTS_INITIALIZER_HPP
