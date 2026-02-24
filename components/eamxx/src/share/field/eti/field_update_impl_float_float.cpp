#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Update,   false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply, false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,   false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Max,      false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min,      false, float, float>(const Field&, const float, const float);

} // namespace scream
