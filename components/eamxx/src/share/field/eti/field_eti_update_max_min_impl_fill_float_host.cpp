#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Host, true, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Host, true, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Max, Host, true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Min, Host, true, float, int>(const Field&, const float, const float);

} // namespace scream
