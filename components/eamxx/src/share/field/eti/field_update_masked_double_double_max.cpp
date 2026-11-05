#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Max, double, double>(const Field&, const double, const double, const double, const Field&) const;

} // namespace scream
