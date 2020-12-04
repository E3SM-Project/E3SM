#ifndef SCREAM_FIELD_INITIALIZER_HPP
#define SCREAM_FIELD_INITIALIZER_HPP

#include "share/grid/grids_manager.hpp"
#include "share/field/field.hpp"
#include "share/scream_types.hpp"

#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"

namespace scream {

/*
 * A class responsible to initialize a field
 *
 * A FieldInitializer object has to be able to initialize one or more fields.
 * This capability will be exploited by the AtmosphereDriver (AD), to make sure
 * all atm inputs are initialized at the beginning of the simulation.
 *
 * If the initializer initializes a field on a grid that is *not* the reference grid,
 * it will also be asked to init a copy of the field on the reference grid.
 * One way to achieve this is to pass to MyInitializer the grids manager, so that 
 * it can build a remapper from its grid(s) to the reference grid.
 * The remapper will then be used to remap the inited fields onto the reference grid.
 */

class FieldInitializer : public ekat::enable_shared_from_this<FieldInitializer> {
public:
  using field_type        = Field<      Real>;
  using const_field_type  = Field<const Real>;
  using remapper_ptr_type = std::shared_ptr<AbstractRemapper<Real>>;

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

  // If f is not on the ref grid, pass both f and f_ref to the initializer.
  // The initializer must initialize BOTH.
  virtual void add_field (const field_type& f, const field_type& f_ref,
                          const remapper_ptr_type& remapper) = 0;
};

// Create a field initializer, and correctly set up the (weak) pointer to self.
template <typename FieldInitType>
inline std::shared_ptr<FieldInitType>
create_field_initializer () {
  auto ptr = std::make_shared<FieldInitType>();
  ptr->setSelfPointer(ptr);
  return ptr;
}

} // namespace scream

#endif // SCREAM_FIELD_INITIALIZER_HPP
