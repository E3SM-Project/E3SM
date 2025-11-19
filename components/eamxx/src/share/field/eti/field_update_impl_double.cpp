#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   false, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,   false, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, false, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   false, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,   false, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, false, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   false, double, int>(const Field&, const double, const double);

} // namespace scream
