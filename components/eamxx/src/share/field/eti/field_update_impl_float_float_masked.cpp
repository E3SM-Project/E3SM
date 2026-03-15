#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<true, CombineMode::Replace,  float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Update,   float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Multiply, float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Divide,   float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Max,      float, float>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Min,      float, float>(const Field&, const float, const float, const Field*);

} // namespace scream
