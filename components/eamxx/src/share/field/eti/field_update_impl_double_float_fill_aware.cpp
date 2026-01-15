#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   true, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Multiply, true, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Divide,   true, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Max,      true, double, float>(const Field&, const double, const double) const;
template void Field::update_impl<CombineMode::Min,      true, double, float>(const Field&, const double, const double) const;

} // namespace scream
