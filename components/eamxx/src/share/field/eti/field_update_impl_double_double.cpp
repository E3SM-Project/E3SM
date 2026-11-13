#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<double, double>(const CombineMode, const Field&, const double, const double, const double) const;

} // namespace scream
