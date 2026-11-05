#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<false,int> (const int, const Field&) const;

} // namespace scream
