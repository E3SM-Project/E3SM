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

  const bool masked1 = f1.has_valid_mask();
  const bool masked2 = f2.has_valid_mask();
  if (masked1!=masked2) {
    // Only one field is masked: the comparison is only well posed if that
    // mask is trivial (all entries valid), since the other field has no
    // matching notion of which entries to ignore.
    auto mask = masked1 ? f1.get_valid_mask() : f2.get_valid_mask();
    auto mask_min = field_min(mask).as<int>();
    if (mask_min<=0) {
      return false;
    }
  } else if (masked1) {
    // Both masked: the masks themselves must agree.
    if (not views_are_equal(f1.get_valid_mask(),f2.get_valid_mask())) {
      return false;
    }
  }

  // Take advantage of field utils update, min and max to assess the max
  // difference between the two fields simply.
  auto ft = f1.clone(CloneFlags::CopyData);
  ft.update(f2,1,-1);
  if (masked1 or masked2) {
    Field mask = masked1 ? f1.get_valid_mask() : f2.get_valid_mask();

    // Zero-out invalid entries (where the two, now known to agree or to be
    // trivial, masks are 0), so they cannot affect the min/max check below,
    // regardless of what garbage (e.g. a fill value) they hold.
    ft.deep_copy(0,mask,true);
  }

  auto d_min = field_min(ft).as<Real>();
  auto d_max = field_max(ft).as<Real>();
  return std::abs(d_min) <= tol and std::abs(d_max) <= tol;
}

} // namespace scream
