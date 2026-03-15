#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<false, CombineMode::Replace,  double, float>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Update,   double, float>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Multiply, double, float>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Divide,   double, float>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Max,      double, float>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Min,      double, float>(const Field&, const double, const double, const Field*);

} // namespace scream
