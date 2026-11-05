#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Max, double, double>(const Field&, const double, const double, const double) const;

} // namespace scream
