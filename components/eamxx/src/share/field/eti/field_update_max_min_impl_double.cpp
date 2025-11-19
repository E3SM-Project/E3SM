#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, false, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, false, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, false, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, false, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, false, double, int>(const Field&, const double, const double);

} // namespace scream
