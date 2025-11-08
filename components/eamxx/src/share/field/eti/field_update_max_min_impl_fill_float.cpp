#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, true, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, true, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Max, true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, true, float, int>(const Field&, const float, const float);

} // namespace scream
