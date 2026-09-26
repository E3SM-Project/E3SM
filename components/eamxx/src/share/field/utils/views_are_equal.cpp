#include "share/field/field_utils.hpp"

#include <cmath>
#include <cstdio>

namespace scream {

bool views_are_equal(const Field& f1, const Field& f2, const Real tol)
{
  const auto& l1 = f1.get_header().get_identifier().get_layout();
  const auto& l2 = f2.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (l1==l2,
      "Error! views_are_equal - the two fields don't have matching layouts.\n");

  // Take advantage of field utils update, min and max to assess the max
  // difference between the two fields simply.
  auto ft = f1.clone(CloneFlags::CopyData);
  ft.update(f2,1,-1);

  // Zero-out invalid entries (where f1 or f2 are invalid),
  // so they cannot affect the min/max check below.
  if (f1.has_valid_mask()) {
    ft.deep_copy(0,f1.get_valid_mask(),true);
  }
  if (f2.has_valid_mask()) {
    ft.deep_copy(0,f2.get_valid_mask(),true);
  }

  auto d_min = field_min(ft).as<Real>();
  auto d_max = field_max(ft).as<Real>();

  return (std::abs(d_min) <= tol) and (std::abs(d_max) <= tol);
}

} // namespace scream
