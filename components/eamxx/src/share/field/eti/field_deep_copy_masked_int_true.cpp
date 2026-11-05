#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_masked<true,int> (const int, const Field&) const;

} // namespace scream
