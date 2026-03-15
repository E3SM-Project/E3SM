#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<false, CombineMode::Replace,  float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Update,   float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Multiply, float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Divide,   float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Max,      float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<false, CombineMode::Min,      float, float>(const Field&, const float, const float, const Field*);

} // namespace scream
