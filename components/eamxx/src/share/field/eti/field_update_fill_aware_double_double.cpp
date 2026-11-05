#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<double, double>(const CombineMode, const Field&, const double, const double, const double) const;

} // namespace scream
