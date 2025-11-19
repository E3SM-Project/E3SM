#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Max, float>(const Field&, const float, const float);
template void Field::update<CombineMode::Min, float>(const Field&, const float, const float);

} // namespace scream
