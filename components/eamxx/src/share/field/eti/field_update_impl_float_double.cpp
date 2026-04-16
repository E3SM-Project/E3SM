#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  float, double>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Update,   float, double>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply, float, double>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,   float, double>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Max,      float, double>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min,      float, double>(const Field&, const float, const float);

} // namespace scream
