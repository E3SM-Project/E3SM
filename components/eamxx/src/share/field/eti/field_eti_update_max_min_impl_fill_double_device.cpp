#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Device, true, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Device, true, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, Device, true, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Device, true, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, Device, true, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Device, true, double, int>(const Field&, const double, const double);

} // namespace scream
