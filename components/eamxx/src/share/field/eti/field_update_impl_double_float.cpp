#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  false, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Update,   false, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Multiply, false, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Divide,   false, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Max,      false, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Min,      false, double, float>(const Field&, const double, const double) const;

} // namespace scream
