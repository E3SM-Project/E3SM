#ifndef SCREAM_CHECK_AND_REPAIR_WRAPPER_HPP
#define SCREAM_CHECK_AND_REPAIR_WRAPPER_HPP

#include "share/field/field_property_check.hpp"

namespace scream
{

class CheckAndRepairWrapper : public FieldPropertyCheck {
public:
  using base_type     = FieldPropertyCheck;
  using base_ptr_type = std::shared_ptr<base_type>;

  // Constructor(s)
  explicit CheckAndRepairWrapper (base_ptr_type check,
                                  base_ptr_type repair = nullptr)
   : m_check  (check)
   , m_repair (repair)
  {
    // repair can be null, but check cannot
    EKAT_REQUIRE_MSG (m_check,
      "Error in CheckAndRepair constructor!\n"
      " Input FieldPropertyCheck for check phase is null.\n");
  }

  // Name of property check - override this method to give a name to the check.
  std::string name () const {
    std::string n;
    n += "CheckAndRepairWrapper\n";
    n += "  - Check : " + m_check->name() + "\n";
    n += "  - Repair: " + (m_repair ? m_repair->name() : "NONE") + "\n";
    return n;
  }

  bool check(const Field& field) const {
    return m_check->check(field);
  }

  bool can_repair() const {
    return m_repair and m_repair->can_repair();
  }

  void repair(Field& field) const {
    EKAT_REQUIRE_MSG (m_repair, "Error! No FieldPropertyCheck stored for repair phase.\n");
    m_repair->repair(field);
  }

private:
  base_ptr_type     m_check;
  base_ptr_type     m_repair;
};

} // namespace scream

#endif // SCREAM_CHECK_AND_REPAIR_WRAPPER_HPP
