#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Update,   false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply, false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,   false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Max,      false, double, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Min,      false, double, double>(const Field&, const double, const double);

} // namespace scream
