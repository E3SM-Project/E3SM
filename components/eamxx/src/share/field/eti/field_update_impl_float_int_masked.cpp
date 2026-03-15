#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<true, CombineMode::Replace,  float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Update,   float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Multiply, float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Divide,   float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Max,      float, int>(const Field&, const float, const float, const Field*);
template void Field::update_impl<true, CombineMode::Min,      float, int>(const Field&, const float, const float, const Field*);

} // namespace scream
