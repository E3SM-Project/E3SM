#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  float, float>(const Field&, const float, const float, const float);
template void Field::update_impl<CombineMode::Update,   float, float>(const Field&, const float, const float, const float);
template void Field::update_impl<CombineMode::Multiply, float, float>(const Field&, const float, const float, const float);
template void Field::update_impl<CombineMode::Divide,   float, float>(const Field&, const float, const float, const float);
template void Field::update_impl<CombineMode::Max,      float, float>(const Field&, const float, const float, const float);
template void Field::update_impl<CombineMode::Min,      float, float>(const Field&, const float, const float, const float);

} // namespace scream
