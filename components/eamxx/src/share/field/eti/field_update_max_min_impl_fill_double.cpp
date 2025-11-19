#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, true, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, true, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, true, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, true, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, true, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, true, double, int>(const Field&, const double, const double);

} // namespace scream
