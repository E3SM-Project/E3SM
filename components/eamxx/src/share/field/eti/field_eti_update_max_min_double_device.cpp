#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Max, Device, double>(const Field&, const double, const double);
template void Field::update<CombineMode::Min, Device, double>(const Field&, const double, const double);

} // namespace scream
