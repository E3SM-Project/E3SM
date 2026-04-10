#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Replace,  float, float>(const Field&, const float, const float);
template void Field::update_fill_aware<CombineMode::Update,   float, float>(const Field&, const float, const float);
template void Field::update_fill_aware<CombineMode::Multiply, float, float>(const Field&, const float, const float);
template void Field::update_fill_aware<CombineMode::Divide,   float, float>(const Field&, const float, const float);
template void Field::update_fill_aware<CombineMode::Max,      float, float>(const Field&, const float, const float);
template void Field::update_fill_aware<CombineMode::Min,      float, float>(const Field&, const float, const float);

} // namespace scream
