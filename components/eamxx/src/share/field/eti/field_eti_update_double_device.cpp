#include "share/field/field.hpp"

namespace scream {

template void Field::update<Device, CombineMode::Update,  double>(const Field&, const double, const double);
template void Field::update<Device, CombineMode::Multiply,double>(const Field&, const double, const double);
template void Field::update<Device, CombineMode::Divide,  double>(const Field&, const double, const double);

} // namespace scream
