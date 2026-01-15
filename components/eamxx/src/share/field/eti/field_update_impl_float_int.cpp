#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  false, float, int>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Update,   false, float, int>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Multiply, false, float, int>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Divide,   false, float, int>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Max,      false, float, int>(const Field&, const float, const float) const;
template void Field::update_impl<CombineMode::Min,      false, float, int>(const Field&, const float, const float) const;

} // namespace scream
