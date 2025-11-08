#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply, false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,   false, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Update,   false, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply, false, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,   false, float, int>(const Field&, const float, const float);

} // namespace scream
