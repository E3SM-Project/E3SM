#ifndef SCREAM_HOMME_INPUTS_INITIALIZER_HPP
#define SCREAM_HOMME_INPUTS_INITIALIZER_HPP

#include "share/field/field_initializer.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

namespace scream {

class HommeInputsInitializer : public FieldInitializer
{
public:
  using field_type = Field<Real>;
  using remapper_ptr_type = std::shared_ptr<AbstractRemapper<Real>>;

  HommeInputsInitializer () = default;

  // Constructor(s) & Destructor
  virtual ~HommeInputsInitializer () = default;

  std::string name () const override { return "HommeInputsInitializer"; }

  void initialize_fields () override;
  const std::set<FieldIdentifier>& get_inited_fields () const override {
    return m_fids;
  }

protected:

  void add_field (const field_type& f) override;
  void add_field (const field_type& f, const field_type& f_ref,
                  const remapper_ptr_type& remapper) override;

  // Members
  std::set<FieldIdentifier>                 m_fids;

  std::shared_ptr<AbstractRemapper<Real>>   m_remapper;

  bool m_fields_inited = false;
};

} // namespace scream

#endif // SCREAM_HOMME_INPUTS_INITIALIZER_HPP
