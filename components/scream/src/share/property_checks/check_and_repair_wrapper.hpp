#ifndef SCREAM_CHECK_AND_REPAIR_WRAPPER_HPP
#define SCREAM_CHECK_AND_REPAIR_WRAPPER_HPP

#include "share/property_checks/property_check.hpp"

#include <memory>

namespace scream
{

/*
 * Use two separate instances for check and repair phases
 * 
 * This class is useful if we want to use two separate
 * checks for the two phases. For instace, we might want
 * to do a non-negativity check and repair, but allow tiny
 * negative values (due to roundoff) in the check phase.
 * We could then use two separate lower bound checks:
 * one with LB<0 small for the check phase, and one with
 * LB=0 for the repair phase.
 */

class CheckAndRepairWrapper : public PropertyCheck {
public:
  using base_type     = PropertyCheck;
  using base_ptr_type = std::shared_ptr<base_type>;

  // Constructor(s)
  explicit CheckAndRepairWrapper (base_ptr_type check,
                                  base_ptr_type repair = nullptr)
   : m_check  (check)
   , m_repair (repair)
  {
    // repair can be null, but check cannot
    EKAT_REQUIRE_MSG (m_check,
      "Error in CheckAndRepairWrapper constructor!\n"
      " Input PropertyCheck for check phase is null.\n");

    // Build list of fields as the union of the check and repair PC lists
    std::list<Field> fields (m_check->fields());
    std::list<bool> repairable (fields.size(),false);

    if (m_repair) {
      for (const auto& f : m_repair->fields()) {
        if (not ekat::contains(fields,f)) {
          fields.push_back(f);
        }
      }
      repairable.resize(fields.size(),false);

      // Mark repairable fields
      for (auto ptr : m_repair->repairable_fields()) {
        auto it = ekat::find(fields,*ptr);
        auto pos = std::distance(fields.begin(),it);
        auto rep_it = std::next(repairable.begin(),pos);
        *rep_it = true;
      }
    }

    set_fields (fields,repairable);
  }

  // Name of property check - override this method to give a name to the check.
  std::string name () const override {
    std::string n;
    n += "CheckAndRepairWrapper\n";
    n += "  - Check : " + m_check->name() + "\n";
    n += "  - Repair: " + (m_repair ? m_repair->name() : "NONE") + "\n";
    return n;
  }

  CheckResult check() const override {
    return m_check->check();
  }

protected:

  void repair_impl() const override {
    EKAT_REQUIRE_MSG (m_repair, "Error! No FieldPropertyCheck stored for repair phase.\n");
    m_repair->repair();
  }

private:
  base_ptr_type     m_check;
  base_ptr_type     m_repair;
};

} // namespace scream

#endif // SCREAM_CHECK_AND_REPAIR_WRAPPER_HPP
