#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Device, false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Device, false, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, Device, false, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Device, false, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, Device, false, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Device, false, double, int>(const Field&, const double, const double);

} // namespace scream
