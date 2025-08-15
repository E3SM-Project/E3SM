#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,  Device, double>(const Field&, const double, const double);
template void Field::update<CombineMode::Multiply,Device, double>(const Field&, const double, const double);
template void Field::update<CombineMode::Divide,  Device, double>(const Field&, const double, const double);

} // namespace scream
