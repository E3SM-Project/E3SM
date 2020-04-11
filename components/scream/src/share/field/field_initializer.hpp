#ifndef SCREAM_FIELD_INITIALIZER_HPP
#define SCREAM_FIELD_INITIALIZER_HPP

#include "share/scream_types.hpp"
#include "share/field/field.hpp"
#include "share/util/scream_std_enable_shared_from_this.hpp"

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
  using field_type = Field<Real,device_type>;

  // A label, mostly to be used for debugging reasons
  virtual const std::string& name () const = 0;

  // The actual method that will be invoked to ensure initialization of field f.
  virtual void initialize_field (const field_type& f) const = 0;

  void add_me_as_initializer (const field_type& f) {
    f.get_header_ptr()->get_tracking().set_initializer(weak_from_this());
  }
};

} // namespace scream

#endif // SCREAM_FIELD_INITIALIZER_HPP
