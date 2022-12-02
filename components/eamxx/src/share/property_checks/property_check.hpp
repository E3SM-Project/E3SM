#ifndef SCREAM_PROPERTY_CHECK_HPP
#define SCREAM_PROPERTY_CHECK_HPP

#include "share/field/field.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>
#include <list>

namespace scream
{

/*
 * Abstract interface for property checks
 * 
 * A property check (PC) object is responsible to check that
 * a certain property holds. The class can (but does not have to)
 * also implement a way to 'repair' the simulation if
 * such property is not satisfied.
 *
 * PC's are stored in an AtmosphereProcess (AP), and can be
 * run before and/or after the process is executed.
 * 
 * The typical class deriving from this interface will implement
 * some checks on one or more Field objects, verifying that they
 * satisfy the properties. For instance, we can check that a
 * Field does not contain NaN values, or that it is always within
 * certain bounds.
 * 
 * More complicate checks can verify that 2+ fields together verify
 * a certain property. For instance, we might want to verify that
 * the sum of certain quantities stay the same (like an energy or
 * mass check).
 */

// A 'Fail' result means the property cannot be fixed,
// while 'Repairable' means the property is "not too far off",
// and can be fixed.
enum class CheckResult {
  Pass,
  Fail,
  Repairable
};

enum class PropertyType {
  PointWise,    // The property is computed pointwise at each entry of the field(s)
  ColumnWise,   // The property is computed over each column of the field(s)
  Global        // The property is computed globally
};

class PropertyCheck {
public:

  virtual ~PropertyCheck()  = default;

  // Name of the property being checked
  virtual std::string name () const = 0;

  struct ResultAndMsg {
    CheckResult       result;
    std::string       msg;

    // For pointwise/columnwise checks, if failing, can return tags/indices
    // for the location where the test failed.
    std::vector<int>        fail_loc_indices;
    std::vector<FieldTag>   fail_loc_tags;
  };

  virtual PropertyType type () const = 0;

  // Check if the property is satisfied, and return true if it is
  virtual ResultAndMsg check () const = 0;

  // Set fields, and whether they can be repaired.
  void set_fields (const std::list<Field>& fields,
                   const std::list<bool>& repairable);

  // Whether this PC is capable of fixing things if the check fails.
  // Defaults to false.
  bool can_repair() const {
    return m_repairable_fields.size()>0;
  }

  // Return the list of fields involved in this property check
  const std::list<Field>& fields () const {
    return m_fields;
  }

  // If can_repair()=true, return a list of the fields that would be
  // repaired. Note that this may be a subset of the list of fields
  // involved in the property check. E.g., for a check f1+f2=C, with
  // C fixed, the repair might be to set f2=C-f1. So this method
  // would return only f2. If repair is not allowed, this method
  // returns an empty list.
  const std::list<Field*>& repairable_fields () const {
    return m_repairable_fields;
  }

  // If a check fails, attempt to repair things. Default is to throw.
  void repair () const;

  // Check the property, and if not satisfied, proceed to repair.
  // The default impl is to run check() and repair() in sequence. If you can
  // perform both in a single call, for performance reason, you should
  // override this method.
  virtual void check_and_repair () const;

  // Whether the input check is the same as this class
  virtual bool same_as (const PropertyCheck& pc) const;

protected:
  virtual void repair_impl () const {
    EKAT_ERROR_MSG ("Error! The method 'repair_impl' has not been overridden.\n"
        "  PropertyCheck name: " + name() + "\n");
  }

  std::list<Field>    m_fields;
  std::list<Field*>   m_repairable_fields;
};

} // namespace scream

#endif // SCREAM_PROPERTY_CHECK_HPP
