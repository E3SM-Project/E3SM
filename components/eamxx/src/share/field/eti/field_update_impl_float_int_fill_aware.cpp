#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply, true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,   true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Max,      true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min,      true, float, int>(const Field&, const float, const float);

} // namespace scream
