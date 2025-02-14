#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Device, false, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, false, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, false, double>(const Field&, const double, const double);

} // namespace scream
