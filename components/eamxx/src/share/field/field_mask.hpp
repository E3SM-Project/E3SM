#ifndef EAMXX_FIELD_MASK_HPP
#define EAMXX_FIELD_MASK_HPP

#include "field.hpp"

namespace scream {

/*
 * A specialization of a Field that represent a mask
 *
 * It is just like a field, but:
 * - data type is hard-coded to IntType
 * - units are hard-coded to nondimensional
 * - overload the update method to align with the idea of combining masks
 */

class FieldMask : public Field
{
public:
  FieldMask (const std::string& name, const FieldLayout& layout, const std::string& grid_name, bool allocate = false);

  // Calls the one above using layout and grid name from the input field f
  FieldMask (const Field& f, const std::string& mname, bool allocate = false);

  FieldMask (const FieldMask& src) = default;
  FieldMask& operator= (const FieldMask& src) = default;

  // If and_op=true, perform bitwise AND, otheriwse perform bitwise OR
  // If negate_lhs=true, negate values of *this before combining.
  // If negate_rhs=true, negate values of x before combining.
  void update(const FieldMask& x, bool and_op, bool negate_lhs, bool negate_rhs);

  // Short version of logical operators
  FieldMask& operator&= (const FieldMask& x);
  FieldMask& operator|= (const FieldMask& x);

  // Note: there is no "in-place" version of operator!, so we must use some name for this.
  // Using operator! is misleading, since for builtin types !a does not change the value stored in a
  FieldMask& flip ();

  // Create a copy of this, then flip it
  FieldMask operator! () const;

  // Clone the mask, returning a completely separate allocation, but storing same values
  FieldMask clone_mask(const std::string& name = "");
};

} // namespace scream

#endif // EAMXX_FIELD_MASK_HPP
