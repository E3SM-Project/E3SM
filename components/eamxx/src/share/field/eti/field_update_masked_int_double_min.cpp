#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Min, int, double>(const Field&, const int, const int, const int, const Field&) const;

} // namespace scream
