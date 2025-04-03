#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Device, false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, false, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,  Device, false, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, false, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, false, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,  Device, false, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, false, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, false, double, int>(const Field&, const double, const double);

} // namespace scream
