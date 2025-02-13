#include "share/field/field.hpp"

namespace scream {

template void Field::update<Device, CombineMode::Update,  float>(const Field&, const float, const float);
template void Field::update<Device, CombineMode::Multiply,float>(const Field&, const float, const float);
template void Field::update<Device, CombineMode::Divide,  float>(const Field&, const float, const float);

} // namespace scream
