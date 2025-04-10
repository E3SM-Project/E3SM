#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Host, true, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply,Host, true, float, float>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,  Host, true, float, float>(const Field&, const float, const float);

template void Field::update_impl<CombineMode::Update,  Host, true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Multiply,Host, true, float, int>(const Field&, const float, const float);
template void Field::update_impl<CombineMode::Divide,  Host, true, float, int>(const Field&, const float, const float);

} // namespace scream
