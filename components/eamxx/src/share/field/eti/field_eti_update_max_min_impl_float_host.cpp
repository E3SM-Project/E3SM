#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Host, false, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Host, false, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Max, Host, false, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Host, false, float, int>(const Field&, const float, const float);

} // namespace scream
