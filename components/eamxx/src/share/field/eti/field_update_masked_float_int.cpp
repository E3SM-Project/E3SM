#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Replace,  float, int>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Update,   float, int>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Multiply, float, int>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Divide,   float, int>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Max,      float, int>(const Field&, const float, const float, const float, const Field&);
template void Field::update_masked<CombineMode::Min,      float, int>(const Field&, const float, const float, const float, const Field&);

} // namespace scream
