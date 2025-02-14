#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,  Device, float>(const Field&, const float, const float);
template void Field::update<CombineMode::Multiply,Device, float>(const Field&, const float, const float);
template void Field::update<CombineMode::Divide,  Device, float>(const Field&, const float, const float);

} // namespace scream
