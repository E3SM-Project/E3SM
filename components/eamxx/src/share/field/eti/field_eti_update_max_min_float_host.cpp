#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Max, Host, float>(const Field&, const float, const float);
template void Field::update<CombineMode::Min, Host, float>(const Field&, const float, const float);

} // namespace scream
