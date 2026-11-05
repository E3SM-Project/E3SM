#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Replace, int, int>(const Field&, const int, const int, const int) const;

} // namespace scream
