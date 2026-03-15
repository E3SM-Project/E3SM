#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  double, double>(const Field&, const double, const double, const Field&);
template void Field::update_impl<CombineMode::Update,   double, double>(const Field&, const double, const double, const Field&);
template void Field::update_impl<CombineMode::Multiply, double, double>(const Field&, const double, const double, const Field&);
template void Field::update_impl<CombineMode::Divide,   double, double>(const Field&, const double, const double, const Field&);
template void Field::update_impl<CombineMode::Max,      double, double>(const Field&, const double, const double, const Field&);
template void Field::update_impl<CombineMode::Min,      double, double>(const Field&, const double, const double, const Field&);

} // namespace scream
