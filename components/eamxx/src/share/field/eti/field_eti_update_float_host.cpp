#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,  Host, float>(const Field&, const float, const float);
template void Field::update<CombineMode::Multiply,Host, float>(const Field&, const float, const float);
template void Field::update<CombineMode::Divide,  Host, float>(const Field&, const float, const float);

} // namespace scream
