#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Device, false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Device, false, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Max, Device, false, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Device, false, float, int>(const Field&, const float, const float);

} // namespace scream
