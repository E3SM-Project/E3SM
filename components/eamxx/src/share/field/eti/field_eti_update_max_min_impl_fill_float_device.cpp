#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Device, true, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Device, true, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Max, Device, true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Device, true, float, int>(const Field&, const float, const float);

} // namespace scream
