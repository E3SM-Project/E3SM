#include "share/property_checks/property_check.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>
#include <list>

namespace scream
{

void PropertyCheck::
set_fields (const std::list<Field>& fields,
            const std::list<bool>& repairable)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (m_fields.size()==0,
      "Error! Cannot reset fields once set.\n"
        "  - PropertyCheck name: " + name() + "\n");

  // Set fields now, since derived classes might need them
  // to create the string returned by name().
  m_fields = fields;

  EKAT_REQUIRE_MSG (fields.size()>0,
      "Error! Input fields list is empty.\n"
        "  - PropertyCheck name: " + name() + "\n");
  EKAT_REQUIRE_MSG (repairable.size()==fields.size(),
      "Error! The method 'set_fields' requires lists of same size.\n"
      "  - Fields list size: " + std::to_string(fields.size()) + "\n"
      "  - Repairable list size: " + std::to_string(repairable.size()) + "\n");
  for (const auto& f : fields) {
    EKAT_REQUIRE_MSG (f.is_allocated(),
        "Error! Fields must be allocated *before* being set in a PropertyCheck.\n"
        "  - PropertyCheck name: " + name() + "\n"
        "  - Field name: " + f.name() + "\n");
  }

  // Do an additional sanity check: the repairable fields must be
  //  - not read-only
  //  - a subset of m_fields
  auto it_f = m_fields.begin();
  auto it_b = repairable.begin();
  for (; it_b!=repairable.end(); ++it_b, ++it_f) {
    if (*it_b) {
      EKAT_REQUIRE_MSG (not it_f->is_read_only(),
          "Error! One of the repairable fields is read only.\n"
          "  - PropertyCheck name: " + name() + "\n"
          "  - Field name: " + it_f->name() + "\n");

      m_repairable_fields.push_back(&(*it_f));
    }
  }
}

// If a check fails, attempt to repair things. Default is to throw.
void PropertyCheck::repair () const {
  EKAT_REQUIRE_MSG (can_repair(),
      "Error! The method 'repair' was called despite can_repair() returns false.\n"
      "  PropertyCheck name: " + name() + "\n");

  repair_impl ();
}

// Check the property, and if not satisfied, proceed to repair.
// The default impl is to run check() and repair() in sequence. If you can
// perform both in a single call, for performance reason, you should
// override this method.
void PropertyCheck::check_and_repair () const {
  EKAT_REQUIRE_MSG(can_repair(),
      "Error! This property check cannot repair.\n"
      "  - PropertyCheck name: " + name() + "\n");
  auto check_result = check();
  if (check_result.result != CheckResult::Fail) {
    repair();
  }
}

bool PropertyCheck::same_as (const PropertyCheck& pc) const
{
  if (this->name()!=pc.name()) {
    return false;
  }

  for (const auto& f : m_fields) {
    bool found = false;
    for (const auto& pc_f : pc.fields()) {
      if (f.get_header().get_identifier()==pc_f.get_header().get_identifier()) {
        found = true;
      }
    }
    if (not found) {
      return false;
    }
  }

  for (const auto& f : m_repairable_fields) {
    bool found = false;
    for (const auto& pc_f : pc.repairable_fields()) {
      if (f->get_header().get_identifier()==pc_f->get_header().get_identifier()) {
        found = true;
      }
    }
    if (not found) {
      return false;
    }
  }
  return true;
}

} // namespace scream
