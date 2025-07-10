#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Host, false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Host, false, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, Host, false, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Host, false, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Max, Host, false, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min, Host, false, double, int>(const Field&, const double, const double);

} // namespace scream
