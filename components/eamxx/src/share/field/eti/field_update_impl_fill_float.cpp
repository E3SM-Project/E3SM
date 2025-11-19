#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   true, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply, true, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,   true, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Update,   true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply, true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,   true, float, int>(const Field&, const float, const float);

} // namespace scream
