#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,   float>(const Field&, const float, const float);
template void Field::update<CombineMode::Multiply, float>(const Field&, const float, const float);
template void Field::update<CombineMode::Divide,   float>(const Field&, const float, const float);

} // namespace scream
