#ifndef SCREAM_FIELD_VALUE_INITIALIZER_HPP
#define SCREAM_FIELD_VALUE_INITIALIZER_HPP

#include "share/field/field_initializer.hpp"

namespace scream {

/*
 * A class responsible to init a field to a given constant
 */

class FieldValueInitializer : public FieldInitializer {
public:

  FieldValueInitializer (const Real v) : m_val(v) {}

  // A label, mostly to be used for debugging reasons
  std::string name () const { return "FieldValueInitializer [" + std::to_string(m_val) + "]"; }

  // The actual method that will be invoked to make this class
  // initialize all the fields it is in charge of.
  void initialize_fields () {
    Kokkos::deep_copy(m_field.get_view(),m_val);
  }

  const std::set<FieldIdentifier>& get_inited_fields () const {
    static std::set<FieldIdentifier> s;
    s.insert(m_field.get_header().get_identifier());
    return s;
  }

protected:

  const Real m_val;

  field_type    m_field;
 
  // Store the field to be init-ed.
  void add_field (const field_type& f) {
    m_field = f;
  }

  void add_field (const field_type& /* f */, const field_type& f_ref,
                  const remapper_ptr_type& /* remapper */) {
    // No need to remap anything, since it's trivial to init f_ref instead of f.
    m_field = f_ref;
  }
};

// Create a field initializer, and correctly set up the (weak) pointer to self.
inline std::shared_ptr<FieldValueInitializer>
create_field_value_initializer (const Real val) {
  auto ptr = std::make_shared<FieldValueInitializer>(val);
  ptr->setSelfPointer(ptr);
  return ptr;
}

} // namespace scream

#endif // SCREAM_FIELD_VALUE_INITIALIZER_HPP
