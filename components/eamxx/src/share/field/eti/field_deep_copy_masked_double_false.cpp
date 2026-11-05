#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<false,double> (const double, const Field&) const;

} // namespace scream
