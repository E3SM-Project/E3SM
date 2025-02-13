#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Device, true, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply,Device, true, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,  Device, true, float>(const Field&, const float, const float);

} // namespace scream
