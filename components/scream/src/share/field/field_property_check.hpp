#ifndef SCREAM_FIELD_PROPERTY_CHECK_HPP
#define SCREAM_FIELD_PROPERTY_CHECK_HPP

namespace scream
{

// Forward declaration of Field.
template<typename ScalarType, typename Device> class Field;

// =================== FIELD PROPERTY CHECK ======================== //

// A Field can have zero or more "property check" objects associated with it.
// Each of these objects performs a check on a field to verify it has a certain
// property. If the property is not satisfied, user can requests that the
// property check try to fix it.
//
// A property can be as simple as the field's values being bounded, or can be
// more involved, referencing other fields (e.g., |f1 - f2| < C), and can even
// involve differential operators (e.g., |div(f)| < eps).
//
// FieldPropertyCheck is an abstract base class that provides an interface to
// be implemented by a subclass.
template<typename ScalarType, typename Device>
class FieldPropertyCheck {
public:

  // Constructor(s)
  FieldPropertyCheck () = default;

  // Override this method to perform a property check on a Field. The method
  // returns true if the property check passes, and false if it fails.
  virtual bool check(const Field<ScalarType, DeviceType>& field) const = 0;

  // Override this method to return true if the property check is capable of
  // attempting to fix a field to make it satisfy a property check.
  virtual bool can_repair() const = 0;

  // Override this method to attempt to repair a field that doesn't pass this
  // property check. The field must be checked again to determine whether the
  // repair is successful.
  virtual void repair(Field<ScalarType, DeviceType>& field) = 0;

};

} // namespace scream

#endif // SCREAM_FIELD_PROPERTY_CHECK_HPP
