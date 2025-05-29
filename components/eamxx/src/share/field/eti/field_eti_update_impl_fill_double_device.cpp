#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Device, true, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, true, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, true, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,  Device, true, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, true, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, true, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,  Device, true, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Device, true, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Device, true, double, int>(const Field&, const double, const double);

} // namespace scream
