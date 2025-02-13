#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Device, true, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, true, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, true, double>(const Field&, const double, const double);

} // namespace scream
