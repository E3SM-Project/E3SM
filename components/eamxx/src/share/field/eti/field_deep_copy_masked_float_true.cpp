#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<true,float> (const float, const Field&) const;

} // namespace scream
