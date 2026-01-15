#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   true, float, float>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Multiply, true, float, float>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Divide,   true, float, float>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Max,      true, float, float>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Min,      true, float, float>(const Field&, const float, const float) const;

} // namespace scream
