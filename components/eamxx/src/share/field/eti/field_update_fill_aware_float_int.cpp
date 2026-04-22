#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Replace,  float, int>(const Field&, const float, const float, const float);
template void Field::update_fill_aware<CombineMode::Update,   float, int>(const Field&, const float, const float, const float);
template void Field::update_fill_aware<CombineMode::Multiply, float, int>(const Field&, const float, const float, const float);
template void Field::update_fill_aware<CombineMode::Divide,   float, int>(const Field&, const float, const float, const float);
template void Field::update_fill_aware<CombineMode::Max,      float, int>(const Field&, const float, const float, const float);
template void Field::update_fill_aware<CombineMode::Min,      float, int>(const Field&, const float, const float, const float);

} // namespace scream
