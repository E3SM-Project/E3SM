#ifndef SCREAM_FIELD_PROPERTY_CHECK_HPP
#define SCREAM_FIELD_PROPERTY_CHECK_HPP

namespace scream
{

// Forward declaration of Field.
template<typename RealType> class Field;

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
template<typename RealType>
class FieldPropertyCheck {
public:

  // Constructor(s)
  FieldPropertyCheck () = default;

  // No copy constructor.
  FieldPropertyCheck(const FieldPropertyCheck&) = delete;

  // No assignment operator, either.
  FieldPropertyCheck& operator=(const FieldPropertyCheck&) = delete;

  virtual ~FieldPropertyCheck() {}

  // Override this method to perform a property check on a Field. The method
  // returns true if the property check passes, and false if it fails.
  virtual bool check(const Field<RealType>& field) const = 0;

  // Override this method to return true if the property check is capable of
  // attempting to fix a field to make it satisfy a property check.
  virtual bool can_repair() const = 0;

  // Override this method to attempt to repair a field that doesn't pass this
  // property check. The field must be checked again to determine whether the
  // repair is successful. NOTE the const in this method!
  virtual void repair(Field<RealType>& field) const = 0;

  // Override this method to provide a more efficient way to check and repair
  // a field in a single pass (applicable only to property checks that can
  // repair a field).
  virtual void check_and_repair(Field<RealType>& field) const {
    if (!check(field) && can_repair()) {
      repair(field);
    }
  }

};

} // namespace scream

#endif // SCREAM_FIELD_PROPERTY_CHECK_HPP
