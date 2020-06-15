#ifndef SCREAM_FIELD_INITIALIZER_HPP
#define SCREAM_FIELD_INITIALIZER_HPP

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_std_enable_shared_from_this.hpp"
#include "share/field/field.hpp"

namespace scream {

/*
 * A class responsible to initialize a field
 *
 * A FieldInitializer object has to be able to initialize
 * one or more field. This capability will be exploited by
 * the AtmosphereDriver (AD), to make sure all atm inputs
 * are initialized at the beginning of the simulation.
 *
 * Most likely, some atm proc will also inherit from this class,
 * since they will claim the role of initializer for some fields.
 */

class FieldInitializer : public util::enable_shared_from_this<FieldInitializer> {
public:
  using device_type = DefaultDevice; // may need to template class on this

  using field_type       = Field<      Real,device_type>;
  using const_field_type = Field<const Real,device_type>;

  virtual ~FieldInitializer () = default;

  // A label, mostly to be used for debugging reasons
  virtual std::string name () const = 0;

  // The actual method that will be invoked to make this class
  // initialize all the fields it is in charge of.
  virtual void initialize_fields () = 0;

  virtual const std::set<FieldIdentifier>& get_inited_fields () const = 0;

  // Note: a const_field_type& is NOT a "field_type const&".
  //       Instead, it is a field with a const value type.
  //       There is an implicit conversion from nonconst to const,
  //       so we can have a single method accepting const_field_type
  void add_me_as_initializer (const const_field_type& f) {
    f.get_header_ptr()->get_tracking().set_initializer(weak_from_this());
  }

  // Store the field to be init-ed.
  virtual void add_field (const field_type& f) = 0;
};

// Create a field initializer, and correctly set up the (weak) pointer to self.
template <typename FieldInitType>
inline std::shared_ptr<FieldInitializer>
create_field_initializer () {
  auto ptr = std::make_shared<FieldInitType>();
  ptr->setSelfPointer(ptr);
  return ptr;
}

} // namespace scream

#endif // SCREAM_FIELD_INITIALIZER_HPP
