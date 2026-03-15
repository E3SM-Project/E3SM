#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Update,   double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Max,      double, int>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min,      double, int>(const Field&, const double, const double);

} // namespace scream
