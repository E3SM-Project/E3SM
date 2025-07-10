#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Max, Device, int>(const Field&, const int, const int);
template void Field::update<CombineMode::Min, Device, int>(const Field&, const int, const int);

} // namespace scream
