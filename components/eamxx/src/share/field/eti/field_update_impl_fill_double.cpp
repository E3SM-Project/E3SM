#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   true, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, true, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   true, double, double>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,   true, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, true, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   true, double, float>(const Field&, const double, const double);

template void Field::update_impl<CombineMode::Update,   true, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, true, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   true, double, int>(const Field&, const double, const double);

} // namespace scream
