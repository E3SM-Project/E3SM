#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Replace,  double, float>(const Field&, const double, const double, const double);
template void Field::update_fill_aware<CombineMode::Update,   double, float>(const Field&, const double, const double, const double);
template void Field::update_fill_aware<CombineMode::Multiply, double, float>(const Field&, const double, const double, const double);
template void Field::update_fill_aware<CombineMode::Divide,   double, float>(const Field&, const double, const double, const double);
template void Field::update_fill_aware<CombineMode::Max,      double, float>(const Field&, const double, const double, const double);
template void Field::update_fill_aware<CombineMode::Min,      double, float>(const Field&, const double, const double, const double);

} // namespace scream
