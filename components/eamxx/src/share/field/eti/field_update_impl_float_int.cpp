#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<false, CombineMode::Replace,  float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Update,   float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Multiply, float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Divide,   float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Max,      float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Min,      float, int>(const Field&, const float, const float, const Field*);

} // namespace scream
