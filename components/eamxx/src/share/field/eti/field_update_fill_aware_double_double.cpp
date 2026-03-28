#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Replace,  double, double>(const Field&, const double, const double);
template void Field::update_fill_aware<CombineMode::Update,   double, double>(const Field&, const double, const double);
template void Field::update_fill_aware<CombineMode::Multiply, double, double>(const Field&, const double, const double);
template void Field::update_fill_aware<CombineMode::Divide,   double, double>(const Field&, const double, const double);
template void Field::update_fill_aware<CombineMode::Max,      double, double>(const Field&, const double, const double);
template void Field::update_fill_aware<CombineMode::Min,      double, double>(const Field&, const double, const double);

} // namespace scream
