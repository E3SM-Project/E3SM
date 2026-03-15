#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<false, CombineMode::Replace,  double, int>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Update,   double, int>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Multiply, double, int>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Divide,   double, int>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Max,      double, int>(const Field&, const double, const double, const Field*);
template void Field::update_impl<false, CombineMode::Min,      double, int>(const Field&, const double, const double, const Field*);

} // namespace scream
