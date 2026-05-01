#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Replace,  float, float>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Update,   float, float>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Multiply, float, float>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Divide,   float, float>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Max,      float, float>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Min,      float, float>(const Field&, const float, const float, const float, const Field&);

} // namespace scream
