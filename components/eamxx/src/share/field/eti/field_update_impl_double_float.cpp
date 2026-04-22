#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Update,   double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Max,      double, float>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min,      double, float>(const Field&, const double, const double);

} // namespace scream
